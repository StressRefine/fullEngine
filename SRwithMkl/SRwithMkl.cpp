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


//
// SRwithMkl.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "mkl.h"

#include <vector>
#include "SRmodel.h"
#include <direct.h>
#include <stdlib.h>
#include <process.h>


using namespace std;


#define LOOP_COUNT 10


//One global instance of "model";
SRmodel model;

void modelSettings(int minp, int maxp, int freezep, int nits,
	double errtol, bool sacr);

int main(int argc, char *argv[])
{

	SRstring line;
	bool batchMode = false;

	if (argc > 1)
	{
		model.wkdir.Copy(argv[1]);
		batchMode = true;
	}
	else
	{

		char buf[256];
		_getcwd(buf, sizeof(buf));
		SCREENPRINT(" _getcwd: %s", buf);
		model.wkdir = buf;
		model.wkdir += "\\";
	}
	//run from the SRUI folder:
	if (model.wkdir.CompareUseLength("C:\\Users\\rich\\Documents\\Visual Studio"))
		model.wkdir.Copy("C:\\Users\\rich\\Documents\\stressrefineWorking\\");
	SCREENPRINT(" wkdir: %s", model.wkdir.str);
	line = model.wkdir;
	line += "engine.log";
	model.logFile.filename = line;
	model.logFile.Delete();
	SRmachDep::GetTime(line);
	LOGPRINT("%s", line.str);
	LOGPRINT("wkdir: %s", model.wkdir.str);

	//options: max p order, error tolerance:
	int minP = 2, maxP = 8, freezeP = 3, nits = 3;
	double errTol = 5;
	bool sacr = true;
	SRstring tok;

	SRstring tail, filename, foldername;

	double sInitial;
	//initial call to dsecond to eliminiate overhead:
	sInitial = dsecnd();

#ifdef _DEBUG
	while (1)
	{
		if (!SRfile::Existcheck("C:\\Users\\rich\\Desktop\\srsleep.txt"))
			break;
	}
#endif

	line = model.wkdir;
	line += "ModelFileName.txt";
	SRfile modelF;
	if (!modelF.Open(line, SRinputMode))
	{
		LOGPRINT("ModelFileName file %s not found\n", line.str);
		ERROREXIT;
	}

	model.Initialize();
	modelF.GetLine(foldername);
	modelF.Close();

	foldername.Right('\\', tail);
	if (tail.len == 0)
	{
		LOGPRINT("Error in ModelFileName file ");
		ERROREXIT;
	}
	model.fileNameTail = tail;
	filename = foldername;
	filename += "\\";
	filename += tail;
	filename += ".msh";
	if (!model.inputFile.Existcheck(filename.str))
	{
		LOGPRINT("input file %s not found\n", filename.str);
		exit(0);
	}

	model.outdir = foldername;
	model.inputFile.SetFileName(filename);
	LOGPRINT("\nStressRefine\n model: %s\n", tail.str);
	SCREENPRINT("\nStressRefine\n model: %s\n", tail.str);

	model.repFile.filename = model.outdir;
	model.repFile.filename += "\\report.txt";
	if (model.repFile.Existcheck())
		model.repFile.Delete();
	model.repFile.Open(SRoutputMode);
	SRmachDep::GetTime(line);
	model.repFile.PrintLine("StressRefine %s   Model: %s\n", line.str, model.fileNameTail.str);
	model.repFile.Close();

	modelSettings(minP, maxP, freezeP, nits, errTol, sacr);

	filename.Left('.',model.settingsFile.filename);
	model.settingsFile.filename += ".srs";

	sInitial = dsecnd();

	model.Create();

	model.Run();

	model.CleanUp();

	double sElapsed = (dsecnd() - sInitial);
	model.repFile.Open(SRappendMode);
	model.repFile.PrintReturn();
	model.repFile.PrintLine("Adaptive Solution Elapsed sec: %lg\n", sElapsed);
	model.repFile.Close();
	OUTPRINT("\nRun Completed\nAdaptive Solution Elapsed sec: %lg\n", sElapsed);
	LOGPRINT("\nStress Refine Run Completed\nAdaptive Solution Elapsed sec: %lg\n", sElapsed);
	SRutil::TimeStamp();

	LOGPRINT("SUCCESSFUL COMPLETION");

	if (!batchMode)
	{
		SCREENPRINT("hit any char to exit");
		int c = getchar();
	}

	return 0;
}

void modelSettings(int minp, int maxp, int freezep, int nits,
	double errtol, bool sacr)
{
	model.setMinPorder(minp);
	model.setMaxPorder(maxp);
	model.pOrderUniform = false;
	if (freezep > 0)
		model.SetFreezePorder(freezep);
	model.setAdaptLoopMax(nits);
	model.setErrorTolerance(errtol / 100.0);
	model.SetDetectSacr(sacr);
}

