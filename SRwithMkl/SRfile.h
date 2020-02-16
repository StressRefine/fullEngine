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
// SRfile.h: interface for the SRfile class.
//
//////////////////////////////////////////////////////////////////////

#if !(defined SRFILE_INCLUDED)
#define SRFILE_INCLUDED

#include "SRstring.h"
#include "SRutil.h"

#define MAXLINELENGTH 256

#define OUTOPEN SRfile::OpenOutFile
#define OUTCLOSE SRfile::CloseOutFile
#define OUTPRINT SRfile::PrintOutFile
#define OUTPRINTNORET SRfile::PrintOutFileNoReturn
#define OUTPRINTRET SRfile::PrintOutFile("\n")
#define LOGPRINT SRfile::LogPrint
#define FRDPRINT model.cgxFrdFile.PrintLine
#define FRDPRINTNORET model.cgxFrdFile.Print
#define FRDPRINTRET model.cgxFrdFile.PrintReturn()
#define SCREENPRINT SRfile::Screenprint
#define REPPRINT SRfile::PrintRepFile


enum FileOpenMode{ SRinputMode, SRoutputMode, SRappendMode, SRoutbinaryMode, SRinbinaryMode, SRinoutbinaryMode };

class SRfile
{
	friend class SRinput;
public:
	static bool CreateDir(char* name);
	static void GetCurrentDir(SRstring& dir);
	static void Delete(char* name);
	static bool Existcheck(char* name);
	static bool Existcheck(SRstring& name){ return Existcheck(name.str); };
	static bool PrintOutFileNoReturn(char *fmt, ...);
	static bool PrintOutFile(char *fmt, ...);
	static void PrintRepFile(char *fmt, ...);
	static bool LogPrint(char *fmt, ...);
	static bool LogPrintLine(char *fmt, va_list arglist);
	static void OpenOutFile();
	static void CloseOutFile();
	static bool PrintLogFileNoReturn(char *fmt, ...);
	static bool PrintLogFileWithTimeStamp(char *fmt, ...);
	static bool PrintLogFile(char *fmt, ...);
	bool Existcheck();
	bool VPrintLine(char* fmt, va_list arglist);
	bool VPrint(char* fmt, va_list arglist);
	bool SeekBinary(int pos, bool intArg = false);
	bool ReadBinary(int n, void* v, bool intArg = false);
	bool WriteBinary(int n, void* v, bool intArg = false);
	void Delete();
	bool PrintReturn();
	bool Close();
	void ToTop(){ rewind(fileptr); };
	bool GetLine(SRstring& line, bool noSlashN = true);
	bool OutOpenNoFail();
	bool Open(FileOpenMode mode, char* name = NULL);
	bool Open(SRstring& fn, FileOpenMode mode){ return Open(mode, fn.str); };
	bool Print(char* s, ...);
	bool PrintLine(char* s, ...);
	void SetFileName(SRstring& name){ filename = name; };
	int GetFilePos(){ return filePos; };
	SRstring& GetFileName(){ return filename; };
	bool Rename(SRstring& newName, bool overWriteOk = false);
	bool Rename(char* newName, bool overWriteOk = false);
	bool GetBdfLine(SRstring& line, bool& isComment);

	static bool SRfile::Screenprint(char *fmt, ...);
	static bool StatPrintLine(char *fmt, va_list arglist);

	SRfile(){ fileptr = NULL; opened = false; filename = ""; bdfLineSaved = false; };
	~SRfile(){Close();};

	SRstring tmpstr;
	FILE* fileptr;
	bool opened;
	SRstring filename;
	char linebuf[MAXLINELENGTH];
	int filePos;
	bool bdfLineSaved;
	SRstring bdfLineSave;
};
#endif //if !(defined SRFILE_INCLUDED)