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
// SRmachDep.cpp: implementation of the SRmachDep class.
//                Machine Dependent Functions
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <direct.h>
#include <time.h>
#include "SRString.h"
#include "SRmodel.h"
#include <windows.h>
#ifdef _WINDOWS
#endif

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

static char buf[256];

bool SRmachDep::CreateDir(char *name)
{
	//create a directory with full path "name"
	int i = _mkdir(name);
	if (i == 0)
		return true;
	else
		return false;
	return false;
}

void SRmachDep::GetTime(SRstring &t, bool timeOnly)
{
	//return date and time in a string
	//note: this should be portable- time and ctime are standard c functions
	//put it in ifdef if problems on other platforms

	//example:
	//Fri Apr 29 12:25:12 2001

	time_t ltime;
	time( &ltime );
	char buf[256];
	ctime_s(buf, 256, &ltime);
	t = buf;
	if (!timeOnly)
	{
		//strip trailing \n:
		int len = t.len;
		t.str[len - 1] = '\0';
	}
	else
	{
		SRstring line;
		line.Copy(t);
		line.Token(); //skip day of week;
		line.Token(); //skip month;
		line.Token(); //skip day;
		t = line.Token();
	}
}

double SRmachDep::availMemCheck()
{
	//look up the available memory on this machine
	//leave a little mem for other processes so don't mess up machine

	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof (statex);
	GlobalMemoryStatusEx(&statex);
	double bytes = (double) statex.ullAvailPhys;
	bytes *= 0.8;//leave a little mem for other processes so don't mess up machine
	double kb = bytes*1000.0;
	return kb;
}

