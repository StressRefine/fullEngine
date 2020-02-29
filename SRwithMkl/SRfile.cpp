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
// SRfile.cpp: implementation of the SRfile class
//
//////////////////////////////////////////////////////////////////////

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include "SRmachDep.h"
#include "SRfile.h"
#include "SRmodel.h"


extern SRmodel model;

SRfile::SRfile()
{
	fileptr = NULL;
	opened = false;
	filename = "";
}

int SRfile::GetFilePos()
{
	return filePos;
};


bool SRfile::Open(FileOpenMode mode,const char *name)
{
    //open file "name"
    //input:
        //mode = SRinputMode, SRoutputMode, SRappendMode, SRoutbinaryMode, SRinbinaryMode, or SRinoutbinaryMode
    //return:
		//true if file was already opened else flalse
	if(opened)
		return false;

	filePos = 0; //position for binary i/o is at beginning of file for newly opened file

	if (name != NULL)
		filename = name;
	if (filename.getLength() == 0)
		return false;

	if (mode == SRinputMode)
	{
		//unsuccessful trying to open a non-existent file for reading
		if (!Existcheck(filename.getStr()))
			return false;
		FOPEN(fileptr,filename.getStr(), "r");
	}
	else if (mode == SRoutputMode)
	{
		//overwrite existing file for output mode:
		if (Existcheck(filename.getStr()))
			UNLINK(filename.getStr());
		FOPEN(fileptr,filename.getStr(), "w");
	}
	else if (mode == SRappendMode)
		FOPEN(fileptr,filename.getStr(), "a");
	else if (mode == SRoutbinaryMode)
		FOPEN(fileptr,filename.getStr(), "wb");
	else if (mode == SRinbinaryMode)
	{
		//unsuccessful trying to open a non-existent file for reading
		if (!Existcheck(filename.getStr()))
			return false;
		FOPEN(fileptr,filename.getStr(), "rb");
	}
	else if (mode == SRinoutbinaryMode)
		FOPEN(fileptr,filename.getStr(), "rb+");
	else
		return false;

	if(fileptr == NULL)
	{
		return false;
	}
	else
	{
		opened = true;
		return true;
	}
}

void SRfile::ToTop()
{
	rewind(fileptr);
}

bool SRfile::Open(SRstring& fn, FileOpenMode mode)
{
	return Open(mode, fn.getStr());
}

void SRfile::SetFileName(SRstring& name)
{
	filename = name;
}

bool SRfile::GetLine(SRstring &line,bool noSlashN)
{
    //get a line from a file
    //input:
        //noSlashN = true to not return the "\n" character at end of line else false
    //output:
        //line = the fetched line stored as SRstring
    //return
        //true if successful else false (e.g. EOF)
	line.Clear();
	char *tmp,c;
	int len;
	tmp = fgets(linebuf, MAXLINELENGTH, fileptr);
	if (tmp == NULL)
		return false;
	if(noSlashN)
	{
		len = strlen(tmp);;
		c = tmp[len-1];
		if(c == '\n')
			tmp[len-1] = '\0';
	}
	else
	{
		len = strlen(tmp);
		c = tmp[len-1];
		if(c != '\n')
		{
			//append \n if not already there:
			tmp[len] = '\n';
		}
	}
	line = tmp;
	return true;
}

bool SRfile::Close()
{
    //close a file
    //return:
		//false if file was not opened or close is unsuccessful, else true

	if(!opened)
		return false;
	opened = false;
	if(fclose(fileptr) != 0)
		return false;
	fileptr = NULL;
	return true;
}


bool SRfile::Print(const char *fmt,...)
{
//print to file; layer over vfprintf
//input:
    //fmt = format string
    //... = variable input

	va_list arglist;
	va_start(arglist, fmt);
	int ret = vfprintf(fileptr, fmt, arglist);
	va_end(arglist);
	if(ret < 0)
		return false;
	else
		return true;
}

bool SRfile::VPrint(const char *fmt, va_list arglist)
{
	//print to file; layer over vfprintf
	//input:
	//fmt = format string
	//va_list = variable argument list

	int ret = 0;
	vprintf(fmt, arglist);
	ret = vfprintf(fileptr, fmt, arglist);
	va_end(arglist);
	if (ret < 0)
		return false;
	else
		return true;
}

bool SRfile::PrintLine(const char *fmt,...)
{
//print to file; layer over vfprintf. append \n
//input:
    //fmt = format string
    //... = variable input

	va_list arglist;
	va_start(arglist, fmt);
	bool ret = VPrintLine(fmt, arglist);
	va_end(arglist);
	return ret;
}

bool SRfile::VPrintLine(const char *fmt, va_list arglist)
{
//print to file; layer over vfprintf. append \n
//input:
    //fmt = format string
    //va_list = variable argument list

	int len = strlen(fmt);
	int ret = 0;
	if(fmt[len-1] != '\n')
	{
		STRCPY(linebuf, MAXLINELENGTH, fmt);
		STRCAT(linebuf, MAXLINELENGTH, "\n");
		ret = vfprintf(fileptr, linebuf, arglist);
	}
	else
	{
		ret = vfprintf(fileptr, fmt, arglist);
	}
	va_end(arglist);
	if(ret < 0)
		return false;
	else
		return true;
}

bool SRfile::PrintReturn()
{
    //print "\n" to a file
	int ret = fprintf(fileptr,"\n");
	if(ret < 0)
		return false;
	else
		return true;
}

bool SRfile::WriteBinary(int n, void* v, bool intArg)
{
	//binary output to this file
	//input:
		//n = number of words
		//v = data to write
		//intArg = true if v is int vector, false for double
	//return:
		//true is successful else false
	//note:
		//class variable filePos contains current position of file pointer,
		//which is updated after the write
	if (!opened)
		return false;
	int nrem, pos, size;
	if (intArg)
		size = sizeof(int);
	else
		size = sizeof(double);
	int nwords = SRBUFSIZE / size;
	double* dv = (double*)v;
	int* iv = (int*)v;
	nrem = n;
	pos = 0;
	while (1)
	{
		if (nrem >= nwords)
		{
			if (intArg)
			{
				if (fwrite(iv + pos, size, nwords, fileptr) < 0)
					return false;
			}
			else
			{
				if (fwrite(dv + pos, size, nwords, fileptr) < 0)
					return false;
			}
			pos += nwords;
			nrem -= nwords;
		}
		else
		{
			if (intArg)
			{
				if (fwrite(iv + pos, size, nrem, fileptr) < 0)
					return false;
			}
			else
			{
				if (fwrite(dv + pos, size, nrem, fileptr) < 0)
					return false;
			}
			break;
		}
	}
	filePos += n;
	return true;
}

bool SRfile::ReadBinary(int n, void* v, bool intArg)
{
	//binary input from this file
	//input:
		//n = number of words
		//v = data to read
		//intArg = true if v is int vector, false for double
	//return:
		//true is successful else false
	//note:
		//class variable filePos contains current position of file pointer,
		//which is updated after the read

	if (!opened)
		return false;
	int size;
	if (intArg)
		size = sizeof(int);
	else
		size = sizeof(double);
	int nwords = SRBUFSIZE / size;
	int nrem, pos;
	nrem = n;
	pos = 0;
	int nread;
	double* dv = (double*)v;
	int* iv = (int*)v;
	while (1)
	{
		if (nrem >= nwords)
		{
			if (intArg)
			{
				if ((int)fread(iv + pos, size, nwords, fileptr) < nwords)
					return false;
			}
			else
			{
				if ((int)fread(dv + pos, size, nwords, fileptr) < nwords)
					return false;
			}
			pos += nwords;
			nrem -= nwords;
		}
		else
		{
			if (intArg)
			{
				nread = (int)fread(iv + pos, size, nrem, fileptr);
				if (nread < nrem)
					return false;
			}
			else
			{
				if ((int)fread(dv + pos, size, nrem, fileptr) < nrem)
					return false;
			}
			break;
		}
	}
	filePos += n;
	return true;
}

bool SRfile::SeekBinary(int pos, bool intArg)
{
	//reposition file pointer for this file
	//input:
		//pos = new file position
		//intArg = true if file has int data, false for double
	//return:
		//true is successful else false
	//note:
		//class variable filePos is updated to "pos"

	if (!opened)
		return false;
	int size;
	if (intArg)
		size = sizeof(int);
	else
		size = sizeof(double);
	if (pos == filePos)
		return true;
	filePos = pos;
	if (pos == 0)
	{
		rewind(fileptr);
		return true;
	}
	else
	{
		int i, nbytes = pos * size;
		i = fseek(fileptr, nbytes, SEEK_SET);
		if (i)
			return false;
		else
			return true;
	}
}

bool SRfile::OutOpenNoFail()
{
	//Open a file whose name was previously assigned for output
	//return true if open successful else false
	//overwrite existing file for output mode:
	if (Existcheck(filename.getStr()))
		UNLINK(filename.getStr());
	FOPEN(fileptr, filename.getStr(), "w");
	if (fileptr == NULL)
		return false;
	else
	{
		opened = true;
		return true;
	}
}

//static:

bool SRfile::Existcheck(const char* name)
{
	//determine if file "name" exists
	//return true if it exists else false
#ifdef linux
	if (access(name, 0) != -1)
		return true;
	else
		return false;
#else
	if (_access(name, 0) != -1)
		return true;
	else
		return false;
#endif
}



bool SRfile::Existcheck(SRstring& name)
{
	return Existcheck(name.getStr());
}

void SRfile::Delete(const char *name)
{
    //delete file "name"
	if(Existcheck(name))
		UNLINK(name);
}

void SRfile::Delete()
{
    //delete this file
	if(opened)
		Close();
	SRfile::Delete(filename.getStr());
}

bool SRfile::PrintOutFile()
{
	//print blank line to out file

	SRfile* f = &model.outputFile;
	bool opened = f->opened;
	if (!opened)
	{
		if (!f->Open(SRappendMode))
			return false;
	}
	f->PrintReturn();
	if (!opened)
		f->Close();
	return true;
}

bool SRfile::PrintOutFile(const char *fmt, ...)
{
	//print to file using format fmt; append \n
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false

	SRfile* f = &model.outputFile;
	bool opened = f->opened;
	if (!opened)
	{
		if (!f->Open(SRappendMode))
			return false;
	}
	va_list arglist;
	va_start(arglist, fmt);
	bool ret = f->VPrintLine(fmt, arglist);
	if (!opened)
		f->Close();
	return ret;
}

bool SRfile::PrintRepFile(const char* fmt, ...)
{
	//print to file using format fmt; append \n
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false

	SRfile* f = &model.repFile;
	bool opened = f->opened;
	if (!opened)
	{
		if (!f->Open(SRappendMode))
			return false;
	}
	va_list arglist;
	va_start(arglist, fmt);
	bool ret = f->VPrintLine(fmt, arglist);
	if (!opened)
		f->Close();
	return ret;
}

bool SRfile::PrintLogFile(const char* fmt, ...)
{
	//print to file using format fmt; append \n
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false

	SRfile* f = &model.logFile;
	bool opened = f->opened;
	if (!opened)
	{
		if (!f->Open(SRappendMode))
			return false;
	}
	va_list arglist;
	va_start(arglist, fmt);
	bool ret = f->VPrintLine(fmt, arglist);
	if (!opened)
		f->Close();
	return ret;
}

bool SRfile::PrintOutFileNoReturn(const char *fmt, ...)
{
	//print to file using format fmt;
	//input:
	//fmt = format string
	//... variable data to print
	//return:
	//true if successful else false

	SRfile* f = &model.outputFile;
	bool opened = f->opened;
	if (!opened)
	{
		if (!f->Open(SRappendMode))
			return false;
	}
	va_list arglist;
	va_start(arglist, fmt);
	bool ret = f->VPrint(fmt, arglist);
	if (!opened)
		f->Close();
	return ret;
}

bool SRfile::Existcheck()
{
	return SRfile::Existcheck(filename.getStr());
}

void SRfile::OpenOutFile()
{
	model.outputFile.Open(SRappendMode);
}
void SRfile::CloseOutFile()
{
	model.outputFile.Close();
}


bool SRfile::Screenprint(const char *fmt, ...)
{
	//print to cmd screen
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false
	va_list arglist;
	va_start(arglist, fmt);
	int ret = 0;
	ret = vprintf(fmt, arglist);
	return ret;
}

