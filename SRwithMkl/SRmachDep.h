/*
Copyright (c) 2020 Richard King

The stressRefine analysis executable "SRwithMkl" is free software: you can redistribute it and/or modify
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
// SRmachDep.h: interface for the SRmachDep class.
// MACHINE DEPENDENT FUNCTIONS
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRMACHDEP_INCLUDED)
#define SRMACHDEP_INCLUDED

//MACHINE Dependency flags:

enum dlgfilemode{newmode,openmode,saveasmode};

#ifdef linux
#define NOSOLVER
#include <unistd.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#define SPRINTF sprintf
#define SSCANF sscanf
#define slashStr "/"
#define slashChar '/'
#define MDGETCWD getcwd
#define UNLINK unlink
#define ACCESS access
#else
#include <direct.h>
//#include <tchar.h>
#include <stdio.h>
#include "windows.h"
#include  <io.h>
#define SPRINTF sprintf_s
#define SSCANF sscanf_s
#define UNLINK _unlink
#define ACCESS _access
#define slashStr "\\"
#define slashChar '\\'
#define MDGETCWD _getcwd
#endif

#define FOPEN SRmachDep::fileOpen
#define STRCPY SRmachDep::stringCopy
#define STRNCPY SRmachDep::stringNCopy
#define STRCAT SRmachDep::stringCat
#define STRNCAT SRmachDep::stringNCat
#define STRTOK SRmachDep::stringToken
#define STRCMP SRmachDep::stringCmp
#define STRNCMP SRmachDep::stringNCmp
#define STRNICMP SRmachDep::stringNCmp
#define STRICMP SRmachDep::stringICmp

#define SRBUFSIZE 8192

class SRstring;
class SRmachDep
{
public:
	static bool CreateDir(const char* name);
	static void fileOpen(FILE*& fptr, const char* name, const char* mode);
	static void stringCopy(char* dest, int destLen,const char* src);
	static void stringNCopy(char* dest, int destLen,const char* src, int n);
	static void stringCat(char* dest, int destLen,const char* src);
	static void stringNCat(char* dest, int destLen,const char* src, int n);
	static char* stringToken(char* str, char* sep, char** nextToken);
	static int stringCmp(const char* str,const char* str2);
	static int stringICmp(const char* str,const char* str2);
	static int stringNCmp(const char* str,const char* str2, int n);
	static int stringNICmp(const char* str,const char* str2, int n);
	static void Delay(int msec);
};



#endif // !defined(SRMACHDEP_INCLUDED)
