/*
Copyright (c) 2020 Richard King

The stressRefine analysis executable "SRwithMkl" is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"SRwithMkl" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The terms of the GNU General Public License are explained in the file COPYING.txt,
also available at <https://www.gnu.org/licenses/>

Note that "SRwithMkl" makes use of the pardiso sparse solver
in the Intel MKL library, with which it must be linked.
Copyright (c) 2018 Intel Corporation.
You may use and redistribute the Intel MKL library, without modification, provided the conditions of
the Intel Simplified Software license are met:
https://software.intel.com/en-us/license/intel-simplified-software-license

It is perfectly permissable to replace the use of the pardiso software from the MKL library
with an equivalent open-source solver
*/

//////////////////////////////////////////////////////////////////////
//
// SRpardiso.cpp: implementation of the SRpardiso class
//for the intel mkl solver
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "SRmodel.h"
#include "SRinput.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

extern SRmodel model;

#ifdef _DEBUG
static bool solverecho = true;
#endif

#define MAXINT32 2E9

#define BKPDIFF false

SRpardiso::SRpardiso()
{
	needInt64 = false;//this should be false except to force testing of 64 branch of mkl bkp and solver
	maxInt32 = MAXINT32;
	bkpSimpleDiff = false;//set this to true to diff fast vs simple pardiso bookkeeping
}


void SRpardiso::bookkeep()
{
	if (model.UseSimplePardiso())
	{
		bookkeepSimple();

		if (!bkpSimpleDiff)
			return;
	}
	int ndof = 3;
	nzTotal = 0;
	int nfunelmax = model.GetmaxNumElementFunctions();
	int nel = model.GetNumElements();
	AllElEquationNumbers.Allocate(nel);//global equation numbers for each local equation in element, includes "-1" for constrained dofs
	AllElSortedFunNumbers.Allocate(nel); //global function numbers for each local function in element

	int neq = model.GetnumEquations();
	int ngfun = model.GetNumFunctions();
	rowNumels.Allocate(ngfun);
	rowEls.Allocate(ngfun);

	int rowfun, roweq;

	//first pass through elements: count number of elements that own each equation (rowNumels),
	// and fill up ElEquationNumbers and ElSortedNumbers
	for (int e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		int nfun = elem->GetNumFunctions();
		SRmklIntVector* ElSortedFunNumbers = AllElSortedFunNumbers.GetPointer(e);
		ElSortedFunNumbers->Allocate(nfun);
		int neqel = ndof * nfun;
		int leq = 0;
		SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(e);
		ElEquationNumbers->Allocate(neqel);

		for (int lfun = 0; lfun < nfun; lfun++)
		{
			int gfun = elem->GetFunctionNumber(lfun);
			ElSortedFunNumbers->d[lfun] = gfun;
			rowNumels.PlusPlus(gfun);
			for (int dof = 0; dof < 3; dof++)
			{
				int geq = model.GetFunctionEquation(gfun, dof);
				ElEquationNumbers->Put(leq, geq);
				leq++;
			}
		}
		ElSortedFunNumbers->Sort();
	}


	SRmklIntVector* oneRowEls;

	for (int row = 0; row < ngfun; row++)
	{
		oneRowEls = rowEls.GetPointer(row);
		int rownel = rowNumels.Get(row);
		oneRowEls->Allocate(rownel);
	}
	rowNumels.Zero();

	//second pass through elements:fill list of elements that own each function
	for (int e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);

		int nfun = elem->GetNumFunctions();
		for (int lfun = 0; lfun < nfun; lfun++)
		{
			int gfun = elem->GetFunctionNumber(lfun);
			oneRowEls = rowEls.GetPointer(gfun);
			int c = rowNumels.Get(gfun);
			oneRowEls->Put(c, e);
			rowNumels.PlusPlus(gfun);
		}
	}


	if (!fillColind(neq, ngfun))
	{
		needInt64 = true;
		//need int 64 for colind bookkeeping:
		fillColind64(neq, ngfun);
	}

	AllElSortedFunNumbers.Free();

	rowNumels.Free();
	for (int row = 0; row < ngfun; row++)
	{
		oneRowEls = rowEls.GetPointer(row);
		oneRowEls->Free();
	}
	rowEls.Free();

	if (model.UseSimplePardiso() && bkpSimpleDiff)
	{
		int firstrowdiff = -1;
		int firstcoldiff = -1;
		for (int i = 0; i <= neq; i++)
		{
			if (rowIndex.Get(i) != rowIndexSimple.Get(i))
			{
				firstrowdiff = i;
				break;
			}
		}
		for (int i = 0; i < nzTotal; i++)
		{
			if (colIndex.Get(i) != colIndexSimple.Get(i))
			{
				firstcoldiff = i;
				break;
			}
		}
		if (firstcoldiff != -1 || firstrowdiff != -1)
		{
			OUTPRINT(" firstcoldiff: %d firstrowdiff %d", firstcoldiff, firstrowdiff);
			ERROREXIT;
		}
	}
}

bool SRpardiso::fillColind(int neq, int ngfun)
{
	if (needInt64)
		return false;
	int numThreads = model.GetMaxNumThreads();
	SRmklIntVector* oneRowEls;
	rowColinds.Allocate(neq);
	rowNZ.Allocate(neq);
	colIndexNonZeroforOneRowParallel.Allocate(numThreads);
	for (int i = 0; i < numThreads; i++)
	{
		SRmklIntVector* colind1 = colIndexNonZeroforOneRowParallel.GetPointer(i);
		colind1->Allocate(ngfun);
	}
	nzStore.Allocate(numThreads * ngfun);

	double sInit = dsecnd();

	if (numThreads == 1)
	{
		for (int rowfun = 0; rowfun < ngfun; rowfun++)
			fillRowNZByFun(rowfun, ngfun, 0);
	}
	else
	{
#pragma omp parallel
		{
#pragma omp for nowait
			for (int rowfun = 0; rowfun < ngfun; rowfun++)
			{
				int pnum = omp_get_thread_num();
				fillRowNZByFun(rowfun, ngfun, pnum);
			}
		}
	}

	nzStore.Free();
	double sElapsed = (dsecnd() - sInit);
	//OUTPRINT("fillrownz elapsed: %lg\n", sElapsed);
	model.bkpSecs += sElapsed;

	for (int i = 0; i < numThreads; i++)
	{
		SRmklIntVector* colind1 = colIndexNonZeroforOneRowParallel.GetPointer(i);
		colind1->Free();
	}
	colIndexNonZeroforOneRowParallel.Free();

	for (int row = 0; row < neq; row++)
		nzTotal += rowNZ.Get(row);

	if (nzTotal < maxInt32)
	{
		globalStiff.Allocate((int)nzTotal);
		rowIndex.Allocate(neq + 1);
		rowIndex.Put(0, 1);
		for (int row = 1; row <= neq; row++)
		{
			int rit0 = rowIndex.Get(row - 1);
			int nzt = rowNZ.Get(row - 1);
			int rit = rit0 + nzt;
			rowIndex.Put(row, rit);
		}
		colIndex.Allocate(nzTotal);
		int colg = 0;
		int colind;
		for (int row = 0; row < neq; row++)
		{
			SRmklIntVector* colindr = rowColinds.GetPointer(row);
			for (int cnz = 0; cnz < colindr->GetNum(); cnz++)
			{
				colind = (int)colindr->Get(cnz);
				colIndex.Put(colg, colind + 1);//colIndex is passed into parDiso solver, which expects 1-based
				colg++;
			}
		}
	}
	else
	{
		needInt64 = true;
		rowNZ.Free();
	}

	return !needInt64;

}

void SRpardiso::fillColind64(int neq, int ngfun)
{
	int numThreads = model.GetMaxNumThreads();
	SRmklIntVector* oneRowEls;
	rowColinds64.Allocate(neq);
	rowNZ.Allocate(neq);
	colIndexNonZeroforOneRowParallel64.Allocate(numThreads);
	for (MKL_INT64 i = 0; i < numThreads; i++)
	{
		SRmklIntVector64* colind1 = colIndexNonZeroforOneRowParallel64.GetPointer(i);
		colind1->Allocate(ngfun);
	}
	nzStore.Allocate(numThreads * ngfun);

	double sInit = dsecnd();

	if (numThreads == 1)
	{
		for (MKL_INT64 rowfun = 0; rowfun < ngfun; rowfun++)
			fillRowNZByFun64(rowfun, ngfun, 0);
	}
	else
	{
#pragma omp parallel
		{
#pragma omp for nowait
			for (MKL_INT64 rowfun = 0; rowfun < ngfun; rowfun++)
			{
				int pnum = omp_get_thread_num();
				fillRowNZByFun64(rowfun, ngfun, pnum);
			}
		}
	}

	nzStore.Free();
	double sElapsed = (dsecnd() - sInit);
	//OUTPRINT("fillrownz64 elapsed: %lg\n", sElapsed);
	model.bkpSecs += sElapsed;

	for (MKL_INT64 i = 0; i < numThreads; i++)
	{
		SRmklIntVector64* colind1 = colIndexNonZeroforOneRowParallel64.GetPointer(i);
		colind1->Free();
	}
	colIndexNonZeroforOneRowParallel64.Free();

	for (MKL_INT64 row = 0; row < neq; row++)
		nzTotal += rowNZ.Get(row);

	globalStiff64.Allocate(nzTotal);
	rowIndex64.Allocate(neq + 1);
	rowIndex64.Put(0, 1);
	for (MKL_INT64 row = 1; row <= neq; row++)
	{
		MKL_INT64 rit0 = rowIndex64.Get(row - 1);
		MKL_INT64 nzt = rowNZ.Get(row - 1);
		MKL_INT64 rit = rit0 + nzt;
		rowIndex64.Put(row, rit);
	}
	colIndex64.Allocate(nzTotal);
	MKL_INT64 colg = 0;
	MKL_INT64 colind;
	for (MKL_INT64 row = 0; row < neq; row++)
	{
		SRmklIntVector64* colindr = rowColinds64.GetPointer(row);
		for (MKL_INT64 cnz = 0; cnz < colindr->GetNum(); cnz++)
		{
			colind = colindr->Get(cnz);
			colIndex64.Put(colg, colind + 1);//colIndex is passed into parDiso solver, which expects 1-based
			colg++;
		}
	}
}


void SRpardiso::smoothBookkeep()
{
	//bookkeeping for the Intel pardiso solver.
	//version for solving for strain smoothing 

	nzTotal = 0;
	int nfunelmax = model.GetmaxNumElementFunctions();
	SRmklIntVector packedtmp;
	packedtmp.Allocate(nfunelmax);
	int nel = model.GetNumElements();
	AllElEquationNumbers.Allocate(nel);//global equation numbers for each local equation in element, includes "-1" for constrained dofs
	AllElSortedFunNumbers.Allocate(nel); //global equation numbers for each local equation in element, excludes constrained dofs

	int neq = model.GetNumSmoothEquations();
	rowNumels.Allocate(neq);
	rowEls.Allocate(neq);

	int rowfun, roweq;

	//first pass through elements: count number of elements that own each equation (rowNumels),
	// and fill up ElEquationNumbers and ElSortedNumbers
	for (int e = 0; e < nel; e++)
	{
		if (!model.post.getElSmooth(e))
			continue;
		SRelement* elem = model.GetElement(e);
		int nfun = elem->GetNumFunctions();
		int leq = 0;
		SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(e);
		ElEquationNumbers->Allocate(nfun);

		for (int lfun = 0; lfun < nfun; lfun++)
		{
			int gfun = elem->GetFunctionNumber(lfun);
			int geq = model.GetSmoothFunctionEquation(gfun);
			ElEquationNumbers->Put(leq, geq);
			packedtmp.Put(lfun, geq);
			rowNumels.PlusPlus(geq);
			leq++;
		}

		SRmklIntVector* ElSortedFunNumbers = AllElSortedFunNumbers.GetPointer(e);
		ElSortedFunNumbers->Allocate(nfun);
		ElSortedFunNumbers->Copy(packedtmp, nfun);
		ElSortedFunNumbers->Sort();
	}

	packedtmp.Free();

	SRmklIntVector* oneRowEls;

	for (int row = 0; row < neq; row++)
	{
		oneRowEls = rowEls.GetPointer(row);
		int rownel = rowNumels.Get(row);
		oneRowEls->Allocate(rownel);
	}
	rowNumels.Zero();

	//second pass through elements:fill list of elements that own each equation (smoothing) or function (!smoothing) and store as rowEls
	for (int e = 0; e < nel; e++)
	{
		if (!model.post.getElSmooth(e))
			continue;
		SRelement* elem = model.GetElement(e);

		int nfun = elem->GetNumFunctions();
		for (int lfun = 0; lfun < nfun; lfun++)
		{
			int gfun = elem->GetFunctionNumber(lfun);
			int geq = model.GetSmoothFunctionEquation(gfun);
			oneRowEls = rowEls.GetPointer(geq);
			int c = rowNumels.Get(geq);
			oneRowEls->Put(c, e);
			rowNumels.PlusPlus(geq);
		}
	}

	if (!needInt64)
		smoothfillColind(neq);
	else
		smoothfillColind64(neq);

	AllElSortedFunNumbers.Free();

	rowNumels.Free();
	for (int row = 0; row < neq; row++)
	{
		oneRowEls = rowEls.GetPointer(row);
		oneRowEls->Free();
	}
	rowEls.Free();
}

void SRpardiso::smoothfillColind(int neq)
{
	int numThreads = model.GetMaxNumThreads();
	SRmklIntVector* oneRowEls;
	rowColinds.Allocate(neq);
	rowNZ.Allocate(neq);
	colIndexNonZeroforOneRowParallel.Allocate(numThreads);
	for (int i = 0; i < numThreads; i++)
	{
		SRmklIntVector* colind1 = colIndexNonZeroforOneRowParallel.GetPointer(i);
		colind1->Allocate(neq);
	}
	nzStore.Allocate(numThreads * neq);

	double sInit = dsecnd();
	if (numThreads == 1)
	{
		for (int row = 0; row < neq; row++)
		{
			SRmklIntVector* colindr = rowColinds.GetPointer(row);
			fillRowNZForSmoothing(row, neq, colindr, 0);
		}
	}
	else
	{
#pragma omp parallel
		{
#pragma omp for nowait
			for (int row = 0; row < neq; row++)
			{
				int pnum = omp_get_thread_num();
				SRmklIntVector* colindr = rowColinds.GetPointer(row);
				fillRowNZForSmoothing(row, neq, colindr, pnum);
			}
		}
	}

	nzStore.Free();
	double sElapsed = (dsecnd() - sInit);
	//OUTPRINT("fillrownz elapsed: %lg\n", sElapsed);
	model.bkpSecs += sElapsed;

	for (int i = 0; i < numThreads; i++)
	{
		SRmklIntVector* colind1 = colIndexNonZeroforOneRowParallel.GetPointer(i);
		colind1->Free();
	}
	colIndexNonZeroforOneRowParallel.Free();

	for (int row = 0; row < neq; row++)
		nzTotal += rowNZ.Get(row);

	globalStiff.Allocate((int)nzTotal);
	rowIndex.Allocate(neq + 1);
	rowIndex.Put(0, 1);
	for (int row = 1; row <= neq; row++)
	{
		int rit0 = rowIndex.Get(row - 1);
		int nzt = rowNZ.Get(row - 1);
		int rit = rit0 + nzt;
		rowIndex.Put(row, rit);
	}
	colIndex.Allocate(nzTotal);
	int colg = 0;
	int colind;
	for (int row = 0; row < neq; row++)
	{
		SRmklIntVector* colindr = rowColinds.GetPointer(row);
		for (int cnz = 0; cnz < colindr->GetNum(); cnz++)
		{
			colind = (int)colindr->Get(cnz);
			colIndex.Put(colg, colind + 1);//colIndex is passed into parDiso solver, which expects 1-based
			colg++;
		}
	}

}

void SRpardiso::smoothfillColind64(int neq)
{
	int numThreads = model.GetMaxNumThreads();
	SRmklIntVector* oneRowEls;
	rowColinds64.Allocate(neq);
	rowNZ.Allocate(neq);
	colIndexNonZeroforOneRowParallel64.Allocate(numThreads);
	for (MKL_INT64 i = 0; i < numThreads; i++)
	{
		SRmklIntVector64* colind1 = colIndexNonZeroforOneRowParallel64.GetPointer(i);
		colind1->Allocate(neq);
	}
	nzStore.Allocate(numThreads * neq);

	double sInit = dsecnd();
	if (numThreads == 1)
	{
		for (MKL_INT64 row = 0; row < neq; row++)
		{
			SRmklIntVector64* colindr = rowColinds64.GetPointer(row);
			fillRowNZ64ForSmoothing(row, neq, colindr, 0);
		}
	}
	else
	{
#pragma omp parallel
		{
#pragma omp for nowait
			for (MKL_INT64 row = 0; row < neq; row++)
			{
				int pnum = omp_get_thread_num();
				SRmklIntVector64* colindr = rowColinds64.GetPointer(row);
				fillRowNZ64ForSmoothing(row, neq, colindr, pnum);
			}
		}
	}

	nzStore.Free();
	double sElapsed = (dsecnd() - sInit);
	//OUTPRINT("fillrownz elapsed: %lg\n", sElapsed);
	model.bkpSecs += sElapsed;

	for (MKL_INT64 i = 0; i < numThreads; i++)
	{
		SRmklIntVector64* colind1 = colIndexNonZeroforOneRowParallel64.GetPointer(i);
		colind1->Free();
	}
	colIndexNonZeroforOneRowParallel64.Free();

	for (MKL_INT64 row = 0; row < neq; row++)
		nzTotal += rowNZ.Get(row);

	globalStiff64.Allocate((int)nzTotal);
	rowIndex64.Allocate(neq + 1);
	rowIndex64.Put(0, 1);
	for (MKL_INT64 row = 1; row <= neq; row++)
	{
		MKL_INT64 rit0 = rowIndex64.Get(row - 1);
		MKL_INT64 nzt = rowNZ.Get(row - 1);
		MKL_INT64 rit = rit0 + nzt;
		rowIndex64.Put(row, rit);
	}
	colIndex64.Allocate((MKL_INT64)nzTotal);
	MKL_INT64 colg = 0;
	MKL_INT64 colind;
	for (MKL_INT64 row = 0; row < neq; row++)
	{
		SRmklIntVector64* colindr = rowColinds64.GetPointer(row);
		for (MKL_INT64 cnz = 0; cnz < colindr->GetNum(); cnz++)
		{
			colind = colindr->Get(cnz);
			colIndex64.Put(colg, colind + 1);//colIndex is passed into parDiso solver, which expects 1-based
			colg++;
		}
	}

}

void SRpardiso::clear()
{
	rowColinds.Free();
	rowIndex.Free();
	colIndex.Free();
	rowNZ.Free();
	AllElEquationNumbers.Free();
	rowColinds64.Free();
	rowIndex64.Free();
	colIndex64.Free();
}

void SRpardiso::fillRowNZByFun(int rowfun, int ngfun, int pnum)
{
	//fill number of nonzero elements for each row (1 row per global equation)
	//also fill the nonzero column indices for each row
	//input:
	//rowfun = global function number
	//ngfun = number global function
	//pnum = processor number in multi-threaded loop
	//notes:
	//fills class variable rownz with number of nonzero elements for the rows associated with rowfun
	//also fills column indices in class variable rowColinds for the same rows

	SRmklIntVector* colind1 = colIndexNonZeroforOneRowParallel.GetPointer(pnum);
	colind1->Zero();
	int *colindOneData = colind1->d;
	int nz = 0;
	SRmklIntVector* oneRowEls = rowEls.GetPointer(rowfun);
	int*oneRowelData = oneRowEls->d;
	int* nzstoreData = nzStore.d + pnum*ngfun;
	int roweq;
	int colfun;
	int nel = rowNumels.d[rowfun];
	for (int et = 0; et < nel; et++)
	{
		int e = oneRowEls->Get(et);
		SRelement* elem = model.GetElement(e);
		SRmklIntVector* elsortedFunNumbers = AllElSortedFunNumbers.GetPointer(e);
		int elnfun = elsortedFunNumbers->GetNum();
		int elrowfun = elsortedFunNumbers->Find(rowfun);
		if (elrowfun == -1)
			ERROREXIT;
		for (int elcolfun = elrowfun; elcolfun < elnfun; elcolfun++)
		{
			colfun = elsortedFunNumbers->Get(elcolfun);
			if (colindOneData[colfun] == 0)
			{
				colindOneData[colfun] = 1;
				nzstoreData[nz] = colfun;
				nz++;
			}
		}
	}

	bool uncon = model.GetFunUncon(rowfun);
	SRmklIntVector* colindr;
	int nz3 = nz * 3;
	int rowdof;
	int coleq;
	if (uncon)
	{
		roweq = model.GetFunctionEquation(rowfun, 0);
		int nzv[3];
		nzv[0] = nz3;
		nzv[1] = nzv[0] - 1;
		nzv[2] = nzv[0] - 2;
		rowNZ.Put(roweq, nzv[0]);
		rowNZ.Put(roweq + 1, nzv[1]);
		rowNZ.Put(roweq + 2, nzv[2]);
		for (rowdof = 0; rowdof < 3; rowdof++)
		{
			SRmklIntVector* colindr = rowColinds.GetPointer(roweq + rowdof);
			colindr->Allocate(nzv[rowdof]);
			int* colindrData = colindr->d;
			colfun = nzstoreData[0];
			int col = 0;
			if (uncon)
			{
				for (int coldof = rowdof; coldof < 3; coldof++)
				{
					coleq = model.GetFunctionEquation(colfun, coldof);
					colindrData[col] = coleq;
					col++;
				}
				for (int cfun = 1; cfun < nz; cfun++)
				{
					colfun = nzstoreData[cfun];
					coleq = model.GetFunctionEquation(colfun, 0);
					colindrData[col] = coleq;
					colindrData[col + 1] = coleq + 1;
					colindrData[col + 2] = coleq + 2;
					col += 3;
				}
			}
			colindr->Sort();
		}
	}
	else
	{
		for (rowdof = 0; rowdof < 3; rowdof++)
		{
			roweq = model.GetFunctionEquation(rowfun, rowdof);
			if (roweq < 0)
				continue;
			int neqz = nz3 - rowdof;
			SRmklIntVector* colindr = rowColinds.GetPointer(roweq);
			colindr->Allocate(neqz);//this is conservative allocation if the function is owned by elements with constraints but that is minor. Will fix length below.
			int* colindrData = colindr->d;
			colfun = nzstoreData[0];
			int col = 0;
			for (int coldof = rowdof; coldof < 3; coldof++)
			{
				coleq = model.GetFunctionEquation(colfun, coldof);
				if (coleq >= 0)
				{
					colindrData[col] = coleq;
					col++;
				}
			}
			for (int cfun = 1; cfun < nz; cfun++)
			{
				colfun = nzstoreData[cfun];
				for (int coldof = 0; coldof < 3; coldof++)
				{
					coleq = model.GetFunctionEquation(colfun, coldof);
					if (coleq >= 0)
					{
						colindrData[col] = coleq;
						col++;
					}
				}
			}
			colindr->num = col; //could save a little memory by realloc, but it's not worth the time for copying the data
			rowNZ.Put(roweq, col);
			colindr->Sort();
		}
	}
}

void SRpardiso::fillRowNZForSmoothing(int row, int neq, SRmklIntVector* colindr, int pnum)
{
	//fill number of nonzero elements for each row (1 row per global equation)
	//version for solving for smoothing 
	//also fill the nonzero column indices for each row
	//input:
	//row = global equation number
	//neq = number global function
	//pnum = processor number in multi-threaded loop
	//output:
	//colindr = column indices for this row
	//notes:
	//fills class variable rownz with number of nonzero elements for the rows associated with rowfun

	SRmklIntVector* colind1 = colIndexNonZeroforOneRowParallel.GetPointer(pnum);
	colind1->Zero();
	int *colindOneData = colind1->d;
	int nz = 0;
	SRmklIntVector* oneRowEls = rowEls.GetPointer(row);
	int* nzstoreData = nzStore.d + pnum*neq;
	for (int et = 0; et < rowNumels.d[row]; et++)
	{
		int e = oneRowEls->Get(et);
		SRelement* elem = model.GetElement(e);
		SRmklIntVector* ElSortedNumbers = AllElSortedFunNumbers.GetPointer(e);
		int elnz = ElSortedNumbers->GetNum();
		int elroweq = ElSortedNumbers->Find(row);
		int coleq;
		for (int elcoleq = elroweq; elcoleq < elnz; elcoleq++)
		{
			coleq = ElSortedNumbers->Get(elcoleq);
			if (colindOneData[coleq] == 0)
			{
				colindOneData[coleq] = 1;
				nzstoreData[nz] = coleq;
				nz++;
			}
		}
	}
	rowNZ.Put(row, nz);
	colindr->Allocate(nz);
	int* colindrData = colindr->d;
	int coleq;
	for (int c = 0; c < nz; c++)
	{
		coleq = nzstoreData[c];
		colindrData[c] = coleq;
	}
	colindr->Sort();
}

void SRpardiso::fillRowNZByFun64(MKL_INT64 rowfun, MKL_INT64 ngfun, int pnum)
{
	//fill number of nonzero elements for each row (1 row per global equation)
	//64 bit integer version
	//also fill the nonzero column indices for each row
	//input:
	//rowfun = global function number
	//ngfun = number global function
	//pnum = processor number in multi-threaded loop
	//notes:
	//fills class variable rownz with number of nonzero elements for the rows associated with rowfun
	//also fills column indices in class variable rowColinds for the same rows


	SRmklIntVector64* colind1 = colIndexNonZeroforOneRowParallel64.GetPointer(pnum);
	colind1->Zero();
	MKL_INT64* colindOneData = colind1->d;
	MKL_INT64 nz = 0;
	SRmklIntVector* oneRowEls = rowEls.GetPointer(rowfun);
	int* oneRowelData = oneRowEls->d;
	int* nzstoreData = nzStore.d + pnum*ngfun;
	int roweq;
	int colfun;
	int nel = rowNumels.d[rowfun];
	for (MKL_INT64 et = 0; et < nel; et++)
	{
		MKL_INT64 e = oneRowEls->Get(et);
		SRelement* elem = model.GetElement(e);
		SRmklIntVector* elsortedFunNumbers = AllElSortedFunNumbers.GetPointer(e);
		MKL_INT64 elnfun = elsortedFunNumbers->GetNum();
		MKL_INT64 elrowfun = elsortedFunNumbers->Find(rowfun);
		if (elrowfun == -1)
			ERROREXIT;
		for (MKL_INT64 elcolfun = elrowfun; elcolfun < elnfun; elcolfun++)
		{
			colfun = elsortedFunNumbers->Get(elcolfun);
			if (colindOneData[colfun] == 0)
			{
				colindOneData[colfun] = 1;
				nzstoreData[nz] = colfun;
				nz++;
			}
		}
	}

	bool uncon = model.GetFunUncon(rowfun);
	SRmklIntVector64* colindr;
	MKL_INT64 nz3 = nz * 3;
	MKL_INT64 rowdof;
	MKL_INT64 coleq;
	if (uncon)
	{
		roweq = model.GetFunctionEquation(rowfun, 0);
		MKL_INT64 nzv[3];
		nzv[0] = nz3;
		nzv[1] = nzv[0] - 1;
		nzv[2] = nzv[0] - 2;
		rowNZ.Put(roweq, nzv[0]);
		rowNZ.Put(roweq + 1, nzv[1]);
		rowNZ.Put(roweq + 2, nzv[2]);
		for (rowdof = 0; rowdof < 3; rowdof++)
		{
			colindr = rowColinds64.GetPointer(roweq + rowdof);
			colindr->Allocate(nzv[rowdof]);
			MKL_INT64* colindrData = colindr->d;
			colfun = nzstoreData[0];
			MKL_INT64 col = 0;
			if (uncon)
			{
				for (MKL_INT64 coldof = rowdof; coldof < 3; coldof++)
				{
					coleq = model.GetFunctionEquation(colfun, coldof);
					colindrData[col] = coleq;
					col++;
				}
				for (MKL_INT64 cfun = 1; cfun < nz; cfun++)
				{
					colfun = nzstoreData[cfun];
					coleq = model.GetFunctionEquation(colfun, 0);
					colindrData[col] = coleq;
					colindrData[col + 1] = coleq + 1;
					colindrData[col + 2] = coleq + 2;
					col += 3;
				}
			}
			colindr->Sort();
		}
	}
	else
	{
		for (rowdof = 0; rowdof < 3; rowdof++)
		{
			roweq = model.GetFunctionEquation(rowfun, rowdof);
			if (roweq < 0)
				continue;
			MKL_INT64 neqz = nz3 - rowdof;
			colindr = rowColinds64.GetPointer(roweq);
			colindr->Allocate(neqz);//this is conservative allocation if the function is owned by elements with constraints but that is minor. Will fix length below.
			MKL_INT64* colindrData = colindr->d;
			colfun = nzstoreData[0];
			MKL_INT64 col = 0;
			for (MKL_INT64 coldof = rowdof; coldof < 3; coldof++)
			{
				coleq = model.GetFunctionEquation(colfun, coldof);
				if (coleq >= 0)
				{
					colindrData[col] = coleq;
					col++;
				}
			}
			for (MKL_INT64 cfun = 1; cfun < nz; cfun++)
			{
				colfun = nzstoreData[cfun];
				for (MKL_INT64 coldof = 0; coldof < 3; coldof++)
				{
					coleq = model.GetFunctionEquation(colfun, coldof);
					if (coleq >= 0)
					{
						colindrData[col] = coleq;
						col++;
					}
				}
			}
			colindr->num = col; //could save a little memory by realloc, but it's not worth the time for copying the data
			rowNZ.Put(roweq, col);
			colindr->Sort();
		}
	}
}


void SRpardiso::fillRowNZ64ForSmoothing(MKL_INT64 row, MKL_INT64 neq, SRmklIntVector64* colindr, int pnum)
{
	//fill number of nonzero elements for each row (1 row per global equation)
	//64 bit integer version for solving for smoothing 
	//also fill the nonzero column indices for each row
	//input:
	//row = global equation number
	//neq = number global function
	//pnum = processor number in multi-threaded loop
	//output:
	//colindr = column indices for this row
	//notes:
	//fills class variable rownz with number of nonzero elements for the rows associated with rowfun

	SRmklIntVector64* colind1 = colIndexNonZeroforOneRowParallel64.GetPointer(pnum);
	colind1->Zero();
	MKL_INT64* colindOneData = colind1->d;
	MKL_INT64 nz = 0;
	SRmklIntVector* oneRowEls = rowEls.GetPointer(row);
	int* nzstoreData = nzStore.d + pnum*neq;
	for (MKL_INT64 et = 0; et < rowNumels.d[row]; et++)
	{
		MKL_INT64 e = oneRowEls->Get(et);
		SRelement* elem = model.GetElement(e);
		SRmklIntVector* ElSortedNumbers = AllElSortedFunNumbers.GetPointer(e);
		MKL_INT64 elnz = ElSortedNumbers->GetNum();
		MKL_INT64 elroweq = ElSortedNumbers->Find(row);
		MKL_INT64 coleq;
		for (MKL_INT64 elcoleq = elroweq; elcoleq < elnz; elcoleq++)
		{
			coleq = ElSortedNumbers->Get(elcoleq);
			if (colindOneData[coleq] == 0)
			{
				colindOneData[coleq] = 1;
				nzstoreData[nz] = coleq;
				nz++;
			}
		}
	}
	rowNZ.Put(row, nz);
	colindr->Allocate(nz);
	MKL_INT64* colindrData = colindr->d;
	MKL_INT64 coleq;
	for (MKL_INT64 c = 0; c < nz; c++)
	{
		coleq = nzstoreData[c];
		colindrData[c] = coleq;
	}
	colindr->Sort();
}

void SRpardiso::smoothAssemble(int elnum, double *elstiff)
{
	if (needInt64)
	{
		smoothAssemble64(elnum, elstiff);
		return;
	}
	SRelement* elem = model.GetElement(elnum);
	SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(elnum);

	int neqel = elem->GetNumFunctions();
	elem->FillStiffDiag(neqel);
	int rowColLocal;
	int rowColGlobal;
	for (int row = 0; row < neqel; row++)
	{
		int roweq = ElEquationNumbers->Get(row);
		if (roweq < 0)
			continue;
		int col0 = rowIndex.Get(roweq) - 1;//row index is one-based
		SRmklIntVector* colindr = rowColinds.GetPointer(roweq);
		for (int col = 0; col < neqel; col++)
		{
			int coleq = ElEquationNumbers->Get(col);
			if (coleq < roweq)
				continue;
			rowColLocal = elem->GetStiffnessLocation(row, col);
			double kelrowcol = elstiff[rowColLocal];
			rowColGlobal = colindr->Find(coleq) + col0;
			globalStiff.PlusAssign(rowColGlobal, kelrowcol);
		}
	}
}


void SRpardiso::assemble()
{
	if (needInt64)
	{
		assemble64();
		return;
	}

	double* elstiff = NULL;
	SRfile* elfile = model.GetElementStiffnessFile();
	if (!model.areElementsInMemory())
		elfile->Open(SRinbinaryMode);
	int nel = model.GetNumElements();
	for (int e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(e);
		elstiff = model.ReadElementStiffness(elem);
		int neqel = ElEquationNumbers->GetNum();
		int rowColLocal;
		int rowColGlobal;
		for (int row = 0; row < neqel; row++)
		{
			int roweq = ElEquationNumbers->Get(row);
			if (roweq < 0)
				continue;
			int col0 = rowIndex.Get(roweq) - 1;//row index is one-based
			SRmklIntVector* colindr = rowColinds.GetPointer(roweq);
			for (int col = 0; col < neqel; col++)
			{
				int coleq = ElEquationNumbers->Get(col);
				if (coleq < roweq)
					continue;
				rowColLocal = elem->GetStiffnessLocation(row, col);
				double kelrowcol = elstiff[rowColLocal];
				rowColGlobal = colindr->Find(coleq) + col0;
				globalStiff.PlusAssign(rowColGlobal, kelrowcol);
			}
		}
	}

	if (!model.areElementsInMemory())
		elfile->Close();
	for (int row = 0; row < model.GetNumEquations(); row++)
	{
		SRmklIntVector* colindr = rowColinds.GetPointer(row);
		colindr->Free();
	}
	rowColinds.Free();

}

void SRpardiso::assemble64()
{
	double* elstiff = NULL;
	SRfile* elfile = model.GetElementStiffnessFile();
	if (!model.areElementsInMemory())
		elfile->Open(SRinbinaryMode);
	MKL_INT64 nel = model.GetNumElements();
	for (MKL_INT64 e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(e);
		elstiff = model.ReadElementStiffness(elem);
		MKL_INT64 neqel = ElEquationNumbers->GetNum();
		MKL_INT64 rowColLocal;
		MKL_INT64 rowColGlobal;
		for (MKL_INT64 row = 0; row < neqel; row++)
		{
			MKL_INT64 roweq = ElEquationNumbers->Get(row);
			if (roweq < 0)
				continue;
			MKL_INT64 col0 = rowIndex64.Get(roweq) - 1;//row index is one-based
			SRmklIntVector64* colindr = rowColinds64.GetPointer(roweq);
			for (MKL_INT64 col = 0; col < neqel; col++)
			{
				int coleq = ElEquationNumbers->Get(col);
				if (coleq < roweq)
					continue;
				rowColLocal = elem->GetStiffnessLocation(row, col);
				double kelrowcol = elstiff[rowColLocal];
				rowColGlobal = colindr->Find(coleq) + col0;
				globalStiff64.PlusAssign(rowColGlobal, kelrowcol);
			}
		}
	}

	if (!model.areElementsInMemory())
		elfile->Close();

	for (MKL_INT64 row = 0; row < model.GetNumEquations(); row++)
	{
		SRmklIntVector64* colindr = rowColinds64.GetPointer(row);
		colindr->Free();
	}
	rowColinds64.Free();
}

void SRpardiso::solve(bool smoothing)
{
	if (needInt64)
	{
		solve64(smoothing);
		return;
	}
	_MKL_DSS_HANDLE_t parSolveHandle;
	MKL_INT iparm[64];//solver control parameters. iparm[0] = 0 means use defaults
	SRmklIntVector perm;
	MKL_INT neq;
	if (smoothing)
		neq = model.GetNumSmoothEquations();
	else
		neq = model.GetNumEquations();
	SRdoubleVector solutionTmp;
	solutionTmp.Allocate(neq);
	SRdoubleVector solutionTmp2;
	solutionTmp2.Allocate(neq);
	perm.Allocate(neq);
	MKL_INT *permp = perm.GetVector();
	MKL_INT *ja = colIndex.GetVector();
	MKL_INT *ia = rowIndex.GetVector();
	int ptv[64];
	for (int i = 0; i < 64; i++)
		ptv[i] = 0;
	parSolveHandle = (_MKL_DSS_HANDLE_t)ptv;
	MKL_INT maxfct = 1;
	MKL_INT mnum = 1;
	MKL_INT mtype = 2; // sym posdef
	MKL_INT solvphase;
	MKL_INT nrhs = 1;

	double sInit = dsecnd();
	iparm[0] = 0;
	MKL_INT msglvl = 0;//message level, 0 for no stats, 1 for stats to screen;
	MKL_INT error = 0;
	double *a = globalStiff.GetVector();
	double *x = solutionTmp.GetVector();
	double* b = solutionTmp2.GetVector();

	bool oocNeeded = false;

	solvphase = 11; //analyze the matrix;
	pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	int memNeeded = iparm[14]; //see mkl documentation
	int mem2 = iparm[15] + iparm[16];
	if (mem2 > memNeeded)
		memNeeded = mem2;

	double MbytesNeeded = ((double)memNeeded) / 1024.0;

	if (!smoothing)
		OUTPRINT("Memory needed for equation solving (Mbytes): %lg\n", MbytesNeeded);

	double availmem = SRmachDep::availMemCheck();

	LOGPRINT("Pardiso sparse equation solving ");

	if (memNeeded > availmem)
		oocNeeded = true;

	if (oocNeeded)
	{
		LOGPRINT(" not enough memory for incore solution, switching to out of core");
		OUTPRINT(" not enough memory for incore solution, switching to out of core");
		//not enough memeory. try switching to OOC:
		iparm[59] = 2;//flag for ooc
		solvphase = 11; //reanalyze the matrix in OOC mode:
		pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	}

	solvphase = 22; //factor the matrix;
	pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	if (error == -2 && !oocNeeded)
	{
		//try again with ooc on. redo analysis and factor in one pass with solvphase 12
		oocNeeded = true;
		iparm[59] = 2;//flag for ooc
		solvphase = 12; //analysize and factor the matrix;
		pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	}
	if (error != 0)
	{
		LOGPRINT("pardiso factor phase. error = %d", error);
		if (!smoothing)
		{
			OUTPRINT("Probable poorly constrained model. Please check constraints and material properties", error);
			LOGPRINT("Probable poorly constrained model. Please check constraints and material properties", error);
			REPPRINT("Poorly constrained model. Please check constraints and material properties.");
			REPPRINT("If loads are in equilibrium, soft spring stabilization may fix the problem.");
		}
		ERROREXIT;
	}

	solvphase = 33; //backsolve;
	if (smoothing)
	{
		for (int comp = 0; comp < 6; comp++)
		{
			b = model.post.getSmoothedStrainVec(comp);
			pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
			for (int i = 0; i < neq; i++)
				b[i] = x[i];
			if (error != 0)
			{
				LOGPRINT("pardiso backsolve phase. error = %d", error);
				ERROREXIT;
			}
		}

	}
	else
	{
		b = model.GetSolutionVector();
		pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
		if (error != 0)
		{
			LOGPRINT("pardiso backsolve phase. error = %d", error);
			LOGPRINT("Probable poorly constrained model. Please check constraints and material properties", error);
			REPPRINT("Poorly constrained model. Please check constraints and material properties.");
			REPPRINT("If loads are in equilibrium, soft spring stabilization may fix the problem.");
			ERROREXIT;
		}
		model.CopyToSolutionVector(solutionTmp);
	}
	colIndex.Free();
	double sElapsed = (dsecnd() - sInit);
	model.solvSecs += sElapsed;

#ifndef _DEBUG
	//release pardiso memory:
	solvphase = -1;
	pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
#endif
}

void SRpardiso::smoothAssemble64(int elnum, double *elstiff)
{
	SRelement* elem = model.GetElement(elnum);
	SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(elnum);

	int neqel = elem->GetNumFunctions();
	elem->FillStiffDiag(neqel);
	int rowColLocal;
	MKL_INT64 rowColGlobal;
	for (MKL_INT64 row = 0; row < neqel; row++)
	{
		MKL_INT64 roweq = ElEquationNumbers->Get(row);
		if (roweq < 0)
			continue;
		MKL_INT64 col0 = rowIndex64.Get(roweq) - 1;//row index is one-based
		SRmklIntVector64* colindr = rowColinds64.GetPointer(roweq);
		for (MKL_INT64 col = 0; col < neqel; col++)
		{
			MKL_INT64 coleq = ElEquationNumbers->Get(col);
			if (coleq < roweq)
				continue;
			rowColLocal = elem->GetStiffnessLocation(row, col);
			double kelrowcol = elstiff[rowColLocal];
			rowColGlobal = colindr->Find(coleq) + col0;
			globalStiff64.PlusAssign(rowColGlobal, kelrowcol);
		}
	}
}

void SRpardiso::solve64(bool smoothing)
{
	OUTPRINT("large problem size, solving with 64 bit integers\n");
	LOGPRINT("large problem size, solving with 64 bit integers\n");

	_MKL_DSS_HANDLE_t parSolveHandle;
	MKL_INT64 iparm[64];//solver control parameters. iparm[0] = 0 means use defaults
	SRmklIntVector64 perm;

	MKL_INT64 neq;
	if (smoothing)
		neq = model.GetNumSmoothEquations();
	else
		neq = model.GetNumEquations();
	SRdoubleVector solutionTmp;
	solutionTmp.Allocate(neq);
	SRdoubleVector solutionTmp2;
	solutionTmp2.Allocate(neq);
	perm.Allocate(neq);
	MKL_INT64 *permp = perm.GetVector();
	MKL_INT64 *ja = colIndex64.GetVector();
	MKL_INT64 *ia = rowIndex64.GetVector();
	MKL_INT64 ptv[64];
	for (MKL_INT64 i = 0; i < 64; i++)
		ptv[i] = 0;
	parSolveHandle = (_MKL_DSS_HANDLE_t)ptv;
	MKL_INT64 maxfct = 1;
	MKL_INT64 mnum = 1;
	MKL_INT64 mtype = 2; // sym posdef
	MKL_INT64 solvphase;
	MKL_INT64 nrhs = 1;

	double sInit = dsecnd();
	iparm[0] = 0;
	MKL_INT64 msglvl = 0;//message level, 0 for no stats, 1 for stats to screen;
	MKL_INT64 error = 0;
	double *a = globalStiff64.GetVector();
	double *x = solutionTmp.GetVector();
	double* b = solutionTmp2.GetVector();

	bool oocNeeded = false;

	solvphase = 11; //analyze the matrix;
	pardiso_64(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);

	int memNeeded = iparm[14]; //see mkl documentation
	int mem2 = iparm[15] + iparm[16];
	if (mem2 > memNeeded)
		memNeeded = mem2;

	double availmem = SRmachDep::availMemCheck();

	if (memNeeded > availmem)
		oocNeeded = true;

	if (oocNeeded)
	{
		LOGPRINT(" not enough memory for incore solution, switching to out of core");
		OUTPRINT(" not enough memory for incore solution, switching to out of core");
		//not enough memeory. try switching to OOC:
		iparm[59] = 2;//flag for ooc
		solvphase = 11; //reanalyze the matrix in OOC mode:
		pardiso_64(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	}

	solvphase = 22; //factor the matrix;
	pardiso_64(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	if (error == -2 && !oocNeeded)
	{
		//try again with ooc on. redo analysis and factor in one pass with solvphase 12
		oocNeeded = true;
		iparm[59] = 2;//flag for ooc
		solvphase = 12; //analysize and factor the matrix;
		pardiso_64(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	}

	if (error != 0)
	{
		LOGPRINT("pardiso factor phase. error = %d", error);
		if (!smoothing)
		{
			LOGPRINT("Probable poorly constrained model. Please check constraints and material properties", error);
		}
		ERROREXIT;
	}
	solvphase = 33; //backsolve;
	if (smoothing)
	{
		for (int comp = 0; comp < 6; comp++)
		{
			b = model.post.getSmoothedStrainVec(comp);
			pardiso_64(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
			for (int i = 0; i < neq; i++)
				b[i] = x[i];
			if (error != 0)
			{
				LOGPRINT("pardiso backsolve phase. error = %d", error);
				ERROREXIT;
			}
		}

	}
	else
	{
		b = model.GetSolutionVector();
		pardiso_64(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
		if (error != 0)
		{
			LOGPRINT("pardiso backsolve phase. error = %d", error);
			LOGPRINT("Probable poorly constrained model. Please check constraints and material properties", error);
			ERROREXIT;
		}
		model.CopyToSolutionVector(solutionTmp);
	}
	colIndex.Free();
	double sElapsed = (dsecnd() - sInit);
	model.solvSecs += sElapsed;

#ifndef _DEBUG
	//release pardiso memory:
	solvphase = -1;
	pardiso_64(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
#endif
}
