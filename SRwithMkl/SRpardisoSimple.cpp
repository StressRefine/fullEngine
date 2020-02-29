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
// SRpardisoSimple.cpp: implementation of the SRpardisoSimple class
// for the intel mkl solver, simler version for debugging
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "SRmodel.h"
#include "SRinput.h"
#ifndef NOSOLVER

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

extern SRmodel model;

void SRpardiso::bookkeepSimple(bool smoothing)
{
	nzTotal = 0;
	int neqelmax = model.GetmaxNumElementFunctions();
	int nelnzmax = neqelmax*(neqelmax + 1) / 2;
	SRmklIntVector packedEqtmp;
	packedEqtmp.Allocate(nelnzmax);
	int nel = model.GetNumElements();
	AllElEquationNumbers.Allocate(nel);//global equation numbers for each local equation in element, includes "-1" for constrained dofs
	AllElpackedEquationNumbers.Allocate(nel); //global equation numbers for each local equation in element, excludes constrained dofs
	elminrow.Allocate(nel);
	elmaxrow.Allocate(nel);

	int ndof = 3;
	if (smoothing)
		ndof = 1;
	for (int e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		int nfun = elem->GetNumFunctions();
		int neqel = ndof * nfun;
		int leq = 0;
		int nelNz = 0;
		SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(e);
		ElEquationNumbers->Allocate(neqel);
		if (smoothing)
		{
			for (int lfun = 0; lfun < nfun; lfun++)
			{
				int gfun = elem->GetFunctionNumber(lfun);
				int geq = model.GetSmoothFunctionEquation(gfun);
				ElEquationNumbers->Put(leq, geq);
				if (geq >= 0)
				{
					packedEqtmp.Put(nelNz, geq);
					nelNz++;
				}
				leq++;
			}
		}
		else
		{
			for (int lfun = 0; lfun < nfun; lfun++)
			{
				int gfun = elem->GetFunctionNumber(lfun);
				for (int dof = 0; dof < 3; dof++)
				{
					int geq = model.GetFunctionEquation(gfun, dof);
					ElEquationNumbers->Put(leq, geq);
					if (geq >= 0)
					{
						packedEqtmp.Put(nelNz, geq);
						nelNz++;
					}
					leq++;
				}
			}
		}
		SRmklIntVector* ElpackedEquationNumbers = AllElpackedEquationNumbers.GetPointer(e);
		ElpackedEquationNumbers->Allocate(nelNz);
		ElpackedEquationNumbers->Copy(packedEqtmp, nelNz);
		ElpackedEquationNumbers->Sort();
		int elid = elem->GetId();
		elminrow.Put(elid, ElpackedEquationNumbers->Get(0));
		elmaxrow.Put(elid, ElpackedEquationNumbers->Get(nelNz - 1));
	}

	packedEqtmp.Free();
	int neq;
	if (smoothing)
		neq = model.GetNumSmoothEquations();
	else
		neq = model.GetnumEquations();
	colIndexNonZeroforOneRow.Allocate(neq);
	rowColinds.Allocate(neq);
	rowNZ.Allocate(neq);

	if (!fillColindSimple(neq))
		fillColindSimple64(neq);

	elmaxrow.Free();
	elminrow.Free();

	if (bkpSimpleDiff)
	{
		rowIndexSimple.Allocate(neq + 1);
		rowIndexSimple.Copy(rowIndex);
		rowIndex.Free();
		colIndexSimple.Allocate(nzTotal);
		colIndexSimple.Copy(colIndex);
		colIndex.Free();
	}
}

void SRpardiso::fillRowNZSimple(int row, int neq, SRmklIntVector* colindr)
{
	int nz = 0;
	colIndexNonZeroforOneRow.Zero();
	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		int elid = elem->GetId();
		int minrow = elminrow.Get(elid);
		int maxrow = elmaxrow.Get(elid);
		SRmklIntVector* ElpackedEquationNumbers = AllElpackedEquationNumbers.GetPointer(e);
		if (row < minrow || row > maxrow)
			continue;
		int elnz = ElpackedEquationNumbers->GetNum();
		for (int elroweq = 0; elroweq <elnz; elroweq++)
		{
			int roweq = ElpackedEquationNumbers->Get(elroweq);
			if (roweq == row)
			{
				int coleq;
				for (int elcoleq = 0; elcoleq < elnz; elcoleq++)
				{
					coleq = ElpackedEquationNumbers->Get(elcoleq);
					if (coleq < roweq)
						continue;
					if (colIndexNonZeroforOneRow.Get(coleq) == 0)
					{
						colIndexNonZeroforOneRow.Put(coleq, 1);
						nz++;
					}
				}
			}
		}
	}
	rowNZ.Put(row, nz);
	nzTotal += nz;

	colindr->Allocate(nz);
	int cnz = 0;
	for (int c = row; c < neq; c++)
	{
		if (colIndexNonZeroforOneRow.Get(c) != 0)
		{
			colindr->Put(cnz, c);
			cnz++;
		}
	}
	colindr->Sort();
}

bool SRpardiso::fillColindSimple(int neq)
{
	if (needInt64)
		return false;
	for (int row = 0; row < neq; row++)
	{
		SRmklIntVector* colindr = rowColinds.GetPointer(row);
		fillRowNZSimple(row, neq, colindr);
	}
	if (nzTotal < maxInt32)
	{
		AllElpackedEquationNumbers.Free();
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
		needInt64 = true;

	return !needInt64;
}

void SRpardiso::fillColindSimple64(int neq)
{
	//free int 32 data, allocate corresponding MKL_INT64 data:
	colIndexNonZeroforOneRow.Free();
	rowIndex.Free();
	colIndex.Free();
	for (int row = 0; row < neq; row++)
	{
		SRmklIntVector* colindr = rowColinds.GetPointer(row);
		colindr->Free();
	}
	rowColinds.Free();

	colIndexNonZeroforOneRow64.Allocate(neq);

	rowIndex64.Allocate(neq + 1);
	colIndex64.Allocate(nzTotal);
	rowColinds64.Allocate(neq);

	for (int row = 0; row < neq; row++)
	{
		SRmklIntVector64* colindr = rowColinds64.GetPointer(row);
		fillRowNZSimple64(row, neq, colindr);
	}
	AllElpackedEquationNumbers.Free();
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

void SRpardiso::fillRowNZSimple64(MKL_INT64 row, MKL_INT64 neq, SRmklIntVector64* colindr)
{
	MKL_INT64 nz = 0;
	colIndexNonZeroforOneRow64.Zero();
	for (MKL_INT64 e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		MKL_INT64 elid = elem->GetId();
		MKL_INT64 minrow = elminrow.Get(elid);
		MKL_INT64 maxrow = elmaxrow.Get(elid);
		SRmklIntVector* ElpackedEquationNumbers = AllElpackedEquationNumbers.GetPointer(e);
		if (row < minrow || row > maxrow)
			continue;
		MKL_INT64 elnz = ElpackedEquationNumbers->GetNum();
		for (MKL_INT64 elroweq = 0; elroweq <elnz; elroweq++)
		{
			MKL_INT64 roweq = ElpackedEquationNumbers->Get(elroweq);
			if (roweq == row)
			{
				MKL_INT64 coleq;
				for (MKL_INT64 elcoleq = 0; elcoleq < elnz; elcoleq++)
				{
					coleq = ElpackedEquationNumbers->Get(elcoleq);
					if (coleq < roweq)
						continue;
					if (colIndexNonZeroforOneRow64.Get(coleq) == 0)
					{
						colIndexNonZeroforOneRow64.Put(coleq, 1);
						nz++;
					}
				}
			}
		}
	}
	rowNZ.Put(row, nz);
	nzTotal += nz;

	colindr->Allocate(nz);
	MKL_INT64 cnz = 0;
	for (MKL_INT64 c = row; c < neq; c++)
	{
		if (colIndexNonZeroforOneRow64.Get(c) != 0)
		{
			colindr->Put(cnz, c);
			cnz++;
		}
	}
	colindr->Sort();
}

//routines for pardiso solver for faceforcegroup traction smoothing

void SRpardiso::solveFFG(SRFaceForceGroup *ffg)
{
	bookkeepFFG(ffg);
	for (int f = 0; f < ffg->faceIds.GetNum(); f++)
		assembleFFG(ffg, f);

	for (int row = 0; row < ffg->nodeIds.GetNum(); row++)
	{
		SRmklIntVector* colindr = rowColinds.GetPointer(row);
		colindr->Free();
	}
	rowColinds.Free();

	_MKL_DSS_HANDLE_t parSolveHandle;
	MKL_INT iparm[64];//solver control parameters. iparm[0] = 0 means use defaults
	SRmklIntVector perm;
	MKL_INT neq;
	neq = ffg->nodeIds.GetNum();
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

	iparm[0] = 0;
	MKL_INT msglvl = 0;//message level, 0 for no stats, 1 for stats to screen;
	MKL_INT error = 0;
	double *a = ffg->smoothStiff.GetVector();
	SRdoubleVector xv;
	SRdoubleVector bv;
	xv.Allocate(neq);
	bv.Allocate(neq);
	double* x = xv.GetVector();
	double* b = bv.GetVector();;

	bool oocNeeded = false;

	solvphase = 11; //analyze the matrix;
	pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	//not enough memory. handled in full SRpardiso branch
	if (error == -2)
		ERROREXIT;

	LOGPRINT("Pardiso sparse equation solving ");

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
		ERROREXIT;
	}

	solvphase = 33; //backsolve;
	for (int dof = 0; dof < 3; dof++)
	{
		for (int i = 0; i < neq; i++)
			bv.Put(i, ffg->ForceDof.Get(dof, i));
		b = bv.d;
		pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
		if (error != 0)
		{
			LOGPRINT("pardiso backsolve phase. error = %d", error);
			ERROREXIT;
		}
		for (int i = 0; i < neq; i++)
			ffg->ForceDof.Put(dof, i, x[i]);
	}
	xv.Free();
}

void SRpardiso::assembleFFG(SRFaceForceGroup* ffg, int elnum)
{
	int fid = ffg->faceIds.Get(elnum);
	SRface* face = model.GetFace(fid);
	int forceid = face->forceIds.Get(0);
	SRforce* faf = model.GetForce(forceid);
	SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(elnum);
	double *elstiff = faf->StiffMat.d;

	int neqel = face->GetNumNodesTotal();
	int rowColLocal;
	int rowColGlobal;
	for (int row = 0; row < neqel; row++)
	{
		int roweq = ElEquationNumbers->Get(row);
		int col0 = rowIndex.Get(roweq) - 1;//row index is one-based
		SRmklIntVector* colindr = rowColinds.GetPointer(roweq);
		for (int col = 0; col < neqel; col++)
		{
			int coleq = ElEquationNumbers->Get(col);
			if (coleq < roweq)
				continue;
			rowColLocal = getFFGElemStiffLocation(faf, row, col);
			double kelrowcol = elstiff[rowColLocal];
			rowColGlobal = colindr->Find(coleq) + col0;
			ffg->smoothStiff.PlusAssign(rowColGlobal, kelrowcol);
		}
	}
}



int SRpardiso::getFFGElemStiffLocation(SRforce* faf, int row, int col)
{
	if (col >= row)
		return faf->stiffDiag.Get(row) + col - row;
	else
		return faf->stiffDiag.Get(col) + row - col;
}

void SRpardiso::fillRowNZFFG(SRFaceForceGroup *ffg, int row, int neq, SRmklIntVector* colindr)
{
	int nz = 0;
	colIndexNonZeroforOneRow.Zero();
	for (int e = 0; e < ffg->faceIds.GetNum(); e++)
	{
		int minrow = elminrow.Get(e);
		int maxrow = elmaxrow.Get(e);
		SRmklIntVector* ElpackedEquationNumbers = AllElpackedEquationNumbers.GetPointer(e);
		if (row < minrow || row > maxrow)
			continue;
		int elnz = ElpackedEquationNumbers->GetNum();
		for (int elroweq = 0; elroweq <elnz; elroweq++)
		{
			int roweq = ElpackedEquationNumbers->Get(elroweq);
			if (roweq == row)
			{
				int coleq;
				for (int elcoleq = 0; elcoleq < elnz; elcoleq++)
				{
					coleq = ElpackedEquationNumbers->Get(elcoleq);
					if (coleq < roweq)
						continue;
					if (colIndexNonZeroforOneRow.Get(coleq) == 0)
					{
						colIndexNonZeroforOneRow.Put(coleq, 1);
						nz++;
					}
				}
			}
		}
	}
	rowNZ.Put(row, nz);
	nzTotal += nz;

	colindr->Allocate(nz);
	int cnz = 0;
	for (int c = row; c < neq; c++)
	{
		if (colIndexNonZeroforOneRow.Get(c) != 0)
		{
			colindr->Put(cnz, c);
			cnz++;
		}
	}
	colindr->Sort();
}

void SRpardiso::bookkeepFFG(SRFaceForceGroup *ffg)
{
	nzTotal = 0;
	int neqelmax = ffg->nfaceFunMax;
	int nelnzmax = neqelmax*(neqelmax + 1) / 2;
	SRmklIntVector packedEqtmp;
	packedEqtmp.Allocate(nelnzmax);
	int nel = ffg->faceIds.GetNum();
	AllElEquationNumbers.Allocate(nel);
	AllElpackedEquationNumbers.Allocate(nel);
	elminrow.Allocate(nel);
	elmaxrow.Allocate(nel);

	for (int e = 0; e < nel; e++)
	{
		int fid = ffg->faceIds.Get(e);
		SRface* face = model.GetFace(fid);
		int nfun = face->GetNumNodesTotal();
		int nelNz = 0;
		SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(e);
		ElEquationNumbers->Allocate(nfun);
		for (int lfun = 0; lfun < nfun; lfun++)
		{
			int gfun = ffg->faceFunLoc.Get(lfun,e);
			ElEquationNumbers->Put(lfun, gfun);
			packedEqtmp.Put(nelNz, gfun);
			nelNz++;
		}
		SRmklIntVector* ElpackedEquationNumbers = AllElpackedEquationNumbers.GetPointer(e);
		ElpackedEquationNumbers->Allocate(nelNz);
		ElpackedEquationNumbers->Copy(packedEqtmp, nelNz);
		ElpackedEquationNumbers->Sort();
		elminrow.Put(e, ElpackedEquationNumbers->Get(0));
		elmaxrow.Put(e, ElpackedEquationNumbers->Get(nelNz - 1));
	}

	packedEqtmp.Free();
	int neq = ffg->nodeIds.GetNum();
	colIndexNonZeroforOneRow.Allocate(neq);
	rowColinds.Allocate(neq);
	rowNZ.Allocate(neq);

	fillColindFFG(ffg, neq);

	elmaxrow.Free();
	elminrow.Free();
}

void SRpardiso::fillColindFFG(SRFaceForceGroup *ffg, int neq)
{
	for (int row = 0; row < neq; row++)
	{
		SRmklIntVector* colindr = rowColinds.GetPointer(row);
		fillRowNZFFG(ffg, row, neq, colindr);
	}
	AllElpackedEquationNumbers.Free();
	ffg->smoothStiff.Allocate((int)nzTotal);
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

#endif

