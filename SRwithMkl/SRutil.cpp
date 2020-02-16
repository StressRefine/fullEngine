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
// SRutil.cpp: implementation of the SRutil class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <omp.h>
#include "SRmodel.h"

extern SRmodel model;

static bool errorExitCalled = false;

void SRutil::ErrorExit(char *file, int line)
{
	//print error messages and shut down
	//when fatal error occurs
	//input:
	//file = filename where error occurred
	//line = line number where error occurred

	if (errorExitCalled)
		exit(0); //prevent recursion
	errorExitCalled = true;
	REPPRINT("StressRefine abnormal termination");

	SRstring s = file;
	char *t = s.LastChar('\\');
	OUTPRINT("*******************************************************\nFatal Error: \nPlease Contact Customer Support at www.StressRefine.com");
	OUTPRINT("File: %s", file);
	OUTPRINT("line: %d", line);
	OUTPRINT("*******************************************************\n\n");
	LOGPRINT("*******************************************************\nFatal Error: \nPlease Contact Customer Support at www.StressRefine.com");
	LOGPRINT("File: %s", file);
	LOGPRINT("line: %d", line);
	LOGPRINT("*******************************************************\n\n");

	exit(0);
}

void SRutil::SRAssert(char *file, int line, bool expn)
{
	//error exit if expn is not true
	if (!expn)
		SRutil::ErrorExit(file, line);
}

void SRutil::TimeStamp()
{
	//put a time stamp in model output file
	SRstring line;
	SRmachDep::GetTime(line);
	OUTPRINT(line.str);
	OUTPRINTRET;
}

void SRintVector::PushBack(int v)
{
	//append integer "v" to the end of this vector
	SRintVector tmp;
	tmp.Copy(*this);
	num++;
	Allocate(num);
	for (int i = 0; i < tmp.GetNum(); i++)
		d[i] = tmp.d[i];
	d[num - 1] = v;
};

void SRdoubleVector::PushBack(double v)
{
	//append double "v" to the end of this vector
	SRdoubleVector tmp;
	tmp.Copy(*this);
	num++;
	Allocate(num);
	for (int i = 0; i < tmp.GetNum(); i++)
		d[i] = tmp.d[i];
	d[num - 1] = v;
};

void SRdoubleMatrix::PlusAssign(SRdoubleMatrix& that)
{
	//add contains of matrix "that" to this matrix
	//if this matrix has not yet been allocated, 
	//zero it and assign the contents of "that" to it
	if (n==0)
		Allocate(that.n, that.m);
	else
	{
		if (n != that.n || m != that.m)
			ERROREXIT;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			d[i][j] += that.d[i][j];
	}
}

void SRdoubleMatrix::Copy(SRdoubleMatrix& that)
{
	//copy contains of matrix "that" into this matrix
	Allocate(that.n, that.m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			d[i][j] = that.d[i][j];
	}
}

int SRIntCompareFunc(const void *v1, const void *v2)
{
	//compare function for sorting SRintVector
	return *((int *)v1) - *((int *)v2);
}

int SRintVector::Find(int intIn)
{
	//find an integer in this vector using binary search
	//input:
		//intIn = integer value
	//return:
		//location where intIn resides in this vector, -1 if not found
	//note:
		//Sort must be called before using this routine

	int* intOutPtr = (int *)bsearch(&intIn, d, num, sizeof(int), SRIntCompareFunc);
	if (intOutPtr == NULL)
		return -1;
	else
		return intOutPtr - d;
}

void SRintVector::Sort()
{
	//binary sort this vector
	qsort((void *)d, num, sizeof(int), SRIntCompareFunc);
}
