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
// SRfile.h: interface for the SRfile class.
//
//////////////////////////////////////////////////////////////////////

#if !(defined SRFILE_INCLUDED)
#define SRFILE_INCLUDED

#include "SRstring.h"
#include "SRutil.h"

#define MAXLINELENGTH 1024

#define OUTOPEN SRfile::OpenOutFile
#define OUTCLOSE SRfile::CloseOutFile
#define OUTPRINT SRfile::PrintOutFile
#define LOGPRINT SRfile::PrintLogFile
#define REPPRINT SRfile::PrintRepFile
#define OUTPRINTNORET SRfile::PrintOutFileNoReturn
#define SCREENPRINT SRfile::Screenprint
#define FRDPRINT model.cgxFrdFile.PrintLine
#define FRDPRINTNORET model.cgxFrdFile.Print
#define FRDPRINTRET model.cgxFrdFile.PrintReturn()


enum FileOpenMode{ SRinputMode, SRoutputMode, SRappendMode, SRoutbinaryMode, SRinbinaryMode, SRinoutbinaryMode };

class SRfile
{
	friend class SRoutput;
	friend class SRinput;
public:
	SRfile();
	~SRfile(){Close();};
	static void Delete(const char* name);
	static bool Existcheck(const char* name);
	static bool Existcheck(SRstring& name);
	static bool PrintOutFileNoReturn(const char *fmt, ...);
	static bool PrintOutFile(const char* fmt, ...);
	static bool PrintOutFile();
	static bool PrintLogFile(const char* fmt, ...);
	static bool PrintRepFile(const char* fmt, ...);
	static bool Screenprint(const char *fmt, ...);
	static void OpenOutFile();
	static void CloseOutFile();

	bool OutOpenNoFail();
	bool Existcheck();
	bool VPrintLine(const char* fmt, va_list arglist);
	bool VPrint(const char* fmt, va_list arglist);
	bool SeekBinary(int pos, bool intArg = false);
	bool ReadBinary(int n, void* v, bool intArg = false);
	bool WriteBinary(int n, void* v, bool intArg = false);
	int GetFilePos();
	void Delete();
	bool PrintReturn();
	bool Close();
	void ToTop();
	bool GetLine(SRstring& line, bool noSlashN = true);
	bool Open(FileOpenMode mode, const char* name = NULL);
	bool Open(SRstring& fn, FileOpenMode mode);
	bool Print(const char* s, ...);
	bool PrintLine(const char* s, ...);
	void SetFileName(SRstring& name);

	SRstring tmpstr;
	FILE* fileptr;
	bool opened;
	SRstring filename;
	char linebuf[MAXLINELENGTH];
	int filePos;
};
#endif //if !(defined SRFILE_INCLUDED)
