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
// SRsolver.h: interface for the SRsolver class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRSOLVER_INCLUDED)
#define SRSOLVER_INCLUDED

#include "SRmklUtil.h"
#include "SRutil.h"

class SRFaceForceGroup;
class SRforce;

class SRpardiso
{
	friend class SRsolver;
public:
	SRpardiso();
	void fillRowNZForSmoothing(int row, int neq, SRmklIntVector* colindr, int pnum);
	void fillRowNZByFun(int rowfun, int ngfun, int pnum);
	void fillRowNZ64ForSmoothing(MKL_INT64 row, MKL_INT64 neq, SRmklIntVector64* colindr, int pnum);
	void fillRowNZByFun64(MKL_INT64 rowfun, MKL_INT64 ngfun, int pnum);
	bool fillColind(int neq, int ngfun);
	void fillColind64(int neq, int ngfun);
	void smoothfillColind(int neq);
	void smoothfillColind64(int neq);
	void bookkeep();
	void smoothBookkeep();
	void clear();
	void assemble();
	void assemble64();
	void smoothAssemble(int elnum, double *elstiff);
	void smoothAssemble64(int elnum, double *elstiff);
	void solve(bool smoothing = false);
	void solve64(bool smoothing = false);
	//simpler routines for debugging:
	void fillRowNZSimple(int row, int neq, SRmklIntVector* colindr);
	void fillRowNZSimple64(MKL_INT64 row, MKL_INT64 neq, SRmklIntVector64* colindr);
	void bookkeepSimple(bool smoothing = false);
	bool fillColindSimple(int neq);
	void fillColindSimple64(int neq);
	//for multifaceforce smoothing:
	void fillColindFFG(SRFaceForceGroup *ffg, int neq);
	void fillRowNZFFG(SRFaceForceGroup *ffg, int row, int neq, SRmklIntVector* colindr);
	void bookkeepFFG(SRFaceForceGroup *ffg);
	void assembleFFG(SRFaceForceGroup* ffg, int elnum);
	void solveFFG(SRFaceForceGroup *ffg);
	int getFFGElemStiffLocation(SRforce* faf, int row, int col);

private:
	bool needInt64;
	MKL_INT64 nzTotal;
	SRmklIntVector rowNumels;
	SRmklIntVector rowNZ;
	SRmklIntVector colIndexNonZeroforOneRow;
	SRmklIntVector64 colIndexNonZeroforOneRow64;
	SRmklIntVector rowIndex;
	SRmklIntVector colIndex;
	SRmklDoubleVector globalStiff;
	SRmklIntVector64 rowIndex64;
	SRmklIntVector64 colIndex64;
	SRmklDouble64Vector globalStiff64;
	SRvector <SRmklIntVector> AllElEquationNumbers;
	SRvector <SRmklIntVector> AllElSortedFunNumbers;
	SRvector <SRmklIntVector> rowColinds;
	SRvector <SRmklIntVector64> rowColinds64;
	SRvector <SRmklIntVector> rowEls;
	SRvector <SRmklIntVector> colIndexNonZeroforOneRowParallel;
	SRvector <SRmklIntVector64> colIndexNonZeroforOneRowParallel64;
	SRmklIntVector nzStore;
	//data for simpler bookkeeping routines:
	SRvector <SRmklIntVector> AllElpackedEquationNumbers;
	SRmklIntVector elmaxrow;
	SRmklIntVector elminrow;
	double maxInt32;
	//for diffing, back up copies:
	SRmklIntVector rowIndexSimple;
	SRmklIntVector colIndexSimple;
	bool bkpSimpleDiff;
};

class SRsolver
{
	friend class SRmodel;

public:
	SRsolver();
	void Cleanup();
	void DoSolution();
	SRpardiso* parDisoPtr(){ return &parDisoSolver; };

private:
	SRpardiso parDisoSolver;
};

class SRskyline
{
public:
	SRskyline(int n);
	void FillSymDiag();
	bool BackSolve(double x[]);
	bool Decomp();
	bool BackSolve(int neq, int diag[], double stiff[], double x[]);
	bool Decomp(int neq, int diag[], double stiff[]);
	int Location(int row, int col);

	int numEq;
	SRintVector stiffDiag;
	SRdoubleVector stiffMat;
};


#endif //!defined(SRSOLVER_INCLUDED)
