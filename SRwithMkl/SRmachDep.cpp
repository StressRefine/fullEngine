/*
Copyright (c) 2020 Richard King

The stressRefine analysis executable "SRwithMkl" is free software:
you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SRwithMkl is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The terms of the GNU General Public License are explained in the file COPYING.txt,
also available at <https://www.gnu.org/licenses/>
*/

//////////////////////////////////////////////////////////////////////
//
// SRmachDep.cpp: implementation of the SRmachDep class.
//                Machine Dependent Functions
//
//////////////////////////////////////////////////////////////////////

#include "SRstring.h"
#include "SRmodel.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

void SRmachDep::Delay(int msec)
{
#ifdef linux
	usleep(msec * 1000); // takes microseconds
#else
	Sleep(msec);
#endif


}

void SRmachDep::fileOpen(FILE*& fptr, const char* name, const char* mode)
{
#ifdef linux
	fptr = fopen(name, mode);
#else
	fopen_s(&fptr, name, mode);
#endif
}

void SRmachDep::stringCopy(char* dest, int destLen, const char* src)
{
#ifdef linux
	strcpy(dest, src);
#else
	strcpy_s(dest, destLen, src);
#endif
}

void SRmachDep::stringNCopy(char* dest, int destLen, const char* src, int n)
{
#ifdef linux
	strncpy(dest, src, n);
    //make sure dest is null-terminated:
	dest[n] = '\0';
#else
	strncpy_s(dest, destLen, src, n);
#endif
}

void SRmachDep::stringCat(char* dest, int destLen, const char* src)
{
#ifdef linux
	strcat(dest, src);
#else
	strcat_s(dest, destLen, src);
#endif
}

void SRmachDep::stringNCat(char* dest, int destLen, const char* src, int n)
{
#ifdef linux
	strncat(dest, src, n);
#else
	strncat_s(dest, destLen, src, n);
#endif
}

char* SRmachDep::stringToken(char* str, char* sep, char** nextToken)
{
#ifdef linux
	return strtok(str, sep);
#else
	return strtok_s(str, sep, nextToken);
#endif
}

int SRmachDep::stringCmp(const char* str,const char* str2)
{
	return strcmp(str, str2);
}

int SRmachDep::stringICmp(const char* str,const char* str2)
{
#ifdef linux
	return strcasecmp(str, str2);
#else
	return _stricmp(str, str2);
#endif
}

int SRmachDep::stringNCmp(const char* str,const char* str2, int n)
{
	return strncmp(str, str2, n);
}

int SRmachDep::stringNICmp(const char* str,const char* str2,int n)
{
#ifdef linux
	return strncasecmp(str, str2,n);
#else
	return _strnicmp(str, str2,n);
#endif
}

bool SRmachDep::CreateDir(const char* name)
{
	//create a directory with full path "name"
	int i;
#ifdef linux
	i = mkdir(name, 0777);
#else
	i = _mkdir(name);
#endif

	if (i == 0)
		return true;
	else
		return false;
	return false;
}





