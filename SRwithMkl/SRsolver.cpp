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
// SRsolver.cpp: implementation of the SRsolver class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "SRmodel.h"
#include "SRinput.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

extern SRmodel model;

#ifdef _DEBUG
static bool solverecho=true;
#endif

SRsolver::SRsolver()
{
}

void SRsolver::Cleanup()
{
#ifndef NOSOLVER
	parDisoSolver.clear();
#endif
}

void SRsolver::DoSolution()
{
#ifndef NOSOLVER
	LOGPRINT("Solution Setup: bookkeep, assemble, preanalyze stiffness matrix\n");
	parDisoSolver.bookkeep();
	parDisoSolver.assemble();
	parDisoSolver.solve();
#endif
}


//skyline solver

SRskyline::SRskyline(int n)
{
	int matlen = n*(n + 1) / 2 + 1;//+1 because decomp and backsolve are 1-based
	stiffMat.Allocate(matlen);
	stiffDiag.Allocate(n + 2);//+2 because decomp and backsolve are 1-based, and need diag of n+1 
	numEq = n;
}


// 1-based diagonal for full symmetric matrix:
void SRskyline::FillSymDiag()
{
	int n = numEq;
	int diag;
	stiffDiag.Put(1, 1);
	for (int i = 2; i <= n + 1; i++)
	{
		diag = stiffDiag.Get(i - 1) + i -1;
		stiffDiag.Put(i, diag);
	}
}

int SRskyline::Location(int row, int col)
{
	//determine location in stiffness matrix, stored as skyline vector, corresponding to
	//row,col
	//input:
	//row = row of stiffness matrix
	//col = column of stiffness matrix
	//return:
	//location in stiffness matrix
	//+1 is because diag is 1-based for compatibility with Decomp and BackSolve

	return stiffDiag.Get(col + 1) + col - row;
}

bool SRskyline::Decomp(int neq, int diag[], double stiff[])
{
	//LU-decomposition of skyline matrix
	//input:
	//neq=no. of equations
	//diag=vector of pointers to diagonal elements of stiffness matrix
	//stiff=stiffness matrix stored as vector in skyline format
	//output:
	//stiff is overwritten with LU-decomposition

	int n, kn, kl, ku, kh, ic, klt, j, ki, nd, l, k, kk;
	double c, b;
	for (n = 1; n <= neq; n++)
	{
		kn = diag[n];
		kl = kn + 1;
		ku = diag[n + 1] - 1;
		kh = ku - kl;
		if (kh < 0)
		{
			if (stiff[kn] < TINY)
				return false;
			else
				continue;
		}
		if (kh > 0)
		{
			k = n - kh;
			ic = 0;
			klt = ku;
			for (j = 1; j <= kh; j++)
			{
				ic++;
				klt--;
				ki = diag[k];
				nd = diag[k + 1] - ki - 1;
				if (nd > 0)
				{
					kk = MATHMIN(ic, nd);
					c = 0.0;
					for (l = 1; l <= kk; l++)
						c += (stiff[ki + l] * stiff[klt + l]);
					stiff[klt] -= c;
				}
				k++;
			}
		}
		k = n;
		b = 0.0;
		for (kk = kl; kk <= ku; kk++)
		{
			k--;
			ki = diag[k];
			if (stiff[ki] < TINY)
				return false;
			c = stiff[kk] / stiff[ki];
			b += (c*stiff[kk]);
			stiff[kk] = c;
		}
		stiff[kn] -= b;
		if (stiff[kn] < TINY)
			return false;
	}
	return true;
}

bool SRskyline::BackSolve(int neq, int diag[], double stiff[], double x[])
{
	//Back-Solve of skyline matrix
	//input:
	//neq = no. of equations
	//diag = vector of pointers to diagonal elements of stiffness matrix
	//stiff = LU-decomposed stiffness matrix stored as vector in skyline format			
	//x = right-hand side vector
	//output:
	//x = overwritten with solution vector

	int n, k, kl, ku, l, kk;
	double c;
	for (n = 1; n <= neq; n++)
	{
		kl = diag[n] + 1;
		ku = diag[n + 1] - 1;
		if ((ku - kl) >= 0)
		{
			k = n;
			c = 0.0;
			for (kk = kl; kk <= ku; kk++)
			{
				k--;
				c += (stiff[kk] * x[k]);
			}
			x[n] -= c;
		}
	}
	for (n = 1; n <= neq; n++)
	{
		k = diag[n];
		if (stiff[k] < TINY)
			return false;
		x[n] /= stiff[k];
	}
	if (neq == 1)
		return true;
	n = neq;
	for (l = 2; l <= neq; l++)
	{
		kl = diag[n] + 1;
		ku = diag[n + 1] - 1;
		if ((ku - kl) >= 0)
		{
			k = n;
			for (kk = kl; kk <= ku; kk++)
			{
				k--;
				x[k] -= (stiff[kk] * x[n]);
			}
		}
		n--;
	}
	return true;
}

bool SRskyline::Decomp()
{
	//note: stiffDiag, stiffMat already stored 1-based, so don't have to
	//correct by passing in, e.g. stiffDiag-1:

	return Decomp(numEq, stiffDiag.GetVector(), stiffMat.GetVector());
}

bool SRskyline::BackSolve(double x[])
{
	//Back-Solve of skyline matrix stored in class variable stiffMat
	//input:
	//x = right-hand side vector
	//output:
	//x = overwritten with solution vector

	//note: X is stored 0 based, so have to correct by passing in x -1::

	return BackSolve(numEq, stiffDiag.GetVector(), stiffMat.GetVector(), x - 1);
}
