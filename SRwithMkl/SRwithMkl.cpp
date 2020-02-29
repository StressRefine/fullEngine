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

#include <stdlib.h>
#include <vector>
#include "SRmodel.h"
#ifndef NOSOLVER
#include "mkl.h"
#endif
#ifndef linux
#include <direct.h>
#include <process.h>
#endif
#include <stdarg.h> 


using namespace std;

//One global instance of "model";
SRmodel model;

void modelSettings(int minp, int maxp, int freezep, int nits,
	double errtol, bool sacr);

void checkForSleep();

int main(int argc, char *argv[])
{
	SRstring line;

#ifdef _DEBUG
	checkForSleep();
#endif

	if (argc > 1)
		model.wkdir.Copy(argv[1]);
	else
	{

		char buf[256];
		MDGETCWD(buf, sizeof(buf));
		SCREENPRINT(" _getcwd: %s", buf);
		model.wkdir = buf;
		model.wkdir += slashStr;
	}

	SCREENPRINT(" wkdir: %s", model.wkdir.getStr());
	line = model.wkdir;
	line += "engine.log";

	model.logFile.filename = line;
	model.logFile.Delete();

	//options: max p order, error tolerance:
	int minP = 2, maxP = 8, freezeP = 3, nits = 3;
	double errTol = 5;
	bool sacr = true;
	SRstring tok;

	SRstring tail, filename, foldername;


	line = model.wkdir;
	line += "ModelFileName.txt";
	SRfile modelF;
	if (!modelF.Open(line, SRinputMode))
	{
		LOGPRINT("ModelFileName file %s not found\n", line.getStr());
		ERROREXIT;
	}

	model.Initialize();
	modelF.GetLine(foldername);
	modelF.Close();

	foldername.Right(slashChar, tail);
	if (tail.getLength() == 0)
	{
		LOGPRINT("Error in ModelFileName file ");
		ERROREXIT;
	}
	model.fileNameTail = tail;
	filename = foldername;
	filename += slashStr;
	filename += tail;
	filename += ".msh";
	if (!model.inputFile.Existcheck(filename.getStr()))
	{
		LOGPRINT("input file %s not found\n", filename.getStr());
		exit(0);
	}

	model.outdir = foldername;
	model.inputFile.SetFileName(filename);
	LOGPRINT("\nStressRefine\n model: %s\n", tail.getStr());
	SCREENPRINT("\nStressRefine\n model: %s\n", tail.getStr());

	model.repFile.filename = model.outdir;
	model.repFile.filename += slashStr;
	model.repFile.filename += "report.txt";
	if (model.repFile.Existcheck())
		model.repFile.Delete();
	model.repFile.Open(SRoutputMode);
	model.repFile.PrintLine("StressRefine  Model: %s\n", model.fileNameTail.getStr());
	model.repFile.Close();

	modelSettings(minP, maxP, freezeP, nits, errTol, sacr);

	filename.Left('.',model.settingsFile.filename);
	model.settingsFile.filename += ".srs";

	model.Create();

	model.Run();

	model.CleanUp();

	model.repFile.Open(SRappendMode);
	model.repFile.PrintReturn();
	model.repFile.PrintLine("Adaptive Solution Complete\n");
	model.repFile.Close();
	OUTPRINT("\nRun Completed\nAdaptive Solution Complete\n");
	LOGPRINT("\nStress Refine Run Completed\n");

	LOGPRINT("SUCCESSFUL COMPLETION");

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

void checkForSleep()
{
	if (SRfile::Existcheck("srsleep.txt"))
		printf(" sleeping while srsleep exists\n");

	while (1)
	{
		if (!SRfile::Existcheck("srsleep.txt"))
			return;
		SRmachDep::Delay(100);
	}
}

