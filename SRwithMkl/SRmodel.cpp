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
// SRmodel.cpp: implementation of the SRmodel class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <stdlib.h>
#include <omp.h>
#include "SRmodel.h"
#include "SRinput.h"
#include "mkl.h"

//////////////////////////////////////////////////////////////////////

SRmodel::SRmodel()
{
	Initialize();
};

void SRmodel::Initialize()
{
	size = 0.0;
	pOrderUniform = false;
	anybricks = false;
	anytets = false;
	anywedges = false;
	minPorder = MINP;
	maxPorder = 5;
	maxPorderFinalAdapt = MAXP;
	maxPorderLowStress = 6;
	thinEdgeMaxPorder = THINEDGEMAXP;
	maxPinModel = MINP;
	adaptLoopMax = MAXADAPTLOOPS;
	ErrorTolerance = ERRTOL;
	anyEnforcedDisplacement = false;
	thermalForce = NULL;
	detectSacrificialElements = true;
	freezeSacrificalP = 2;
	detectThinElements = false;
	freezeSacrEdges = true;
	elemSecs = 0.0;
	bkpSecs = 0.0;
	solvSecs = 0.0;
	maxNodeUid = 0;
	maxElemUid = 0;
	adaptIt = 0;
	elUidAtMaxError = -1;
	simpleElements = false;
	simplePardiso = false;
	breakout = false;
	breakoutByMat = false;
	breakoutAtNode = false;
	saveBreakout = false;
	allConstraintsAsPenalty = false; //change to true to use penalty constraints for all constraints
	allMatsHaveAllowable = true;
	maxAllowableAnyActiveMat = 0.0;
	strainMax = 0.0;
	numFlattenedElHighStress = 0;
	allMatsEquivElast = true;
	anyMatYielded = false;
	prevStressMax = 0.0;
	maxsp1 = 0.0;
	minsp2 = BIG;
	maxFlattened = 0.0;
	maxFlattenedElUId = -1;
	nodeidAtMaxStressPos = -1;
	errorChecker.Initialize();
	maxPJump = 10;
	feasvmmax = 0.0;
	numShellOrBeamNodes = 0;
	partialFlatten = true;
	maxStressElid = -1;
	needSoftSprings = false;
	stressUnitConversion = 1.0;
	stressUnitstr = "Pa (N/m^2)";
	lengthUnitConversion = 1.0;
	lengthUnitstr = "m";
	anyUnsupportedFace = false;
	ReadDispStressSRR = false;
	useUnits = true;
	checkForHotSpotatMax = true;
	anyHotSpotElFound = false;
	breakoutRadSacr = BIG;
	localAdapt = false;
	outputf06 = false;
	allFacesFlat = false;
	doEnergySmooth = true;
	anyUnsupNode = false;
	anyBsurf = false;
	maxNumThreads = 1;

	math.Setup();
}


static char *stressComp[6] = { "xx", "yy", "zz", "xy", "xz", "yz"};

void SRmodel::Run()
{
	//run analysis on this model

	//1. initial operations outside p-loop
	//2. p-loop:
		//function numbering
		//process constraints
		//equation numbering
		//calculate element stiffnesse
		//process forces
		//process enforced displacements
		//solve:
			//assemble and decomp global matrix
			//backsolve force to calculate displacements
		//adapt:
			//calculate raw stresses and perform smoothing
			//calculate errors:
				//1. smoothed vs raw stresses
				//2. traction jumps across shared faces
				//3. traction jumps vs applied loads
			//use errors to determine next p level
		//loop until error within tolerance, all edges at max p, or max iterations reached
	//3. post-process
		//output mesh, displacement results file, and stress results file. output summary quantities

	double sInitial, sElapsed;
	double errForOutput = 0.0;
	if (pOrderUniform)
		passData.Allocate(8);
	else
		passData.Allocate(adaptLoopMax);

	sumFile.Open(SRoutputMode);

	int nits;
	int lastIt = 0;
	SRvec3 dmax;
	int nodeUidAtMaxDisp = 0;

	if (saveBreakout)
	{
		sumFile.PrintLine("saveBreakout run");
		post.saveBreakoutModel();
		return;
	}

	double availmem = SRmachDep::availMemCheck();
	MaxElementMem = 0.5*availmem;
	int mbytes = (int)(availmem / 1.049E6);
	OUTPRINT("available memory for run: %d MBytes\n", mbytes);

	int its, neq = 0;

	if (freezeSacrificalP > maxPorder)
		freezeSacrificalP = maxPorder;
	if (freezeSacrificalP < minPorder)
		freezeSacrificalP = minPorder;

	bool anyLcsEnfd = PreProcessPenaltyConstraints();

	if (detectThinElements)
	{
		int numthin = 0;
		for (int e = 0; e < GetNumElements(); e++)
		{
			SRelement* elem = GetElement(e);
			if (elem->DetectThinEdges())
				numthin++;
		}
		OUTPRINT("number of thin elements: %d", numthin);
	}

	if (pOrderUniform)
	{
		maxPorder = maxPorderFinalAdapt;
		nits = maxPorder - minPorder + 1;
	}
	else
	{
		nits = adaptLoopMax;
		if (maxPorder == minPorder)
			nits = 1;
	}

	bool adapted = true;

	sInitial = dsecnd();
	double sPrev = sInitial;

	nodeUidAtMaxDisp = 0;


	for (its = 0; its < nits; its++)
	{
		adaptIt = its;

		OUTPRINT("\nAdaptive Solution. Iteration: %d\n", its + 1);
		OUTPRINT("Maximum Polynomial Order in Model: %d\n", maxPinModel);
		sumFile.PrintLine("maxP %d", maxPinModel);
		volume = 0.0;

		NumberGlobalFunctions();

		ProcessConstraints();

		NumberEquations();

		solution.Allocate(numEquations);

		allocateSmallElementData(anyLcsEnfd);

		checkElementMapping();

		LOGPRINT("Calculating Element Stiffness\n");

		sPrev = dsecnd();

		CheckElementsFitInMemory();
		SCREENPRINT("Adaptive Iteration: %d. Calculating element stiffnesses\n", adaptIt + 1);
		CalculateElementStiffnesses(anyLcsEnfd);

		sElapsed = (dsecnd() - sPrev);
		elemSecs += sElapsed;

		if (needSoftSprings)
			AddSoftSprings();

		ProcessForces();

		EnforcedDisplacementAssemble();

		//solution:
		SCREENPRINT("Adaptive Iteration: %d. Solving Equations\n", adaptIt + 1);
		LOGPRINT("Solving %d Equations\n", numEquations);
		OUTPRINT("\nNumber of Equations: %d\n", numEquations);
		sumFile.PrintLine("Equations: %d", numEquations);

		solver.DoSolution();

		NodalMaxDisp(nodeUidAtMaxDisp, dmax);
		OUTPRINT("\nMaximum nodal displacement (at node %d)\n", nodeUidAtMaxDisp);
		OUTPRINT("ux: %lg\nuy: %lg\nuz: %lg\n\n", dmax.d[0], dmax.d[1], dmax.d[2]);
		sumFile.PrintLine("max_nodal_disp: %lg %lg %lg", dmax.d[0], dmax.d[1], dmax.d[2]);

		int pbeforeAdapt = maxPinModel;

		SCREENPRINT("Adaptive Iteration: %d. Stress Smoothing and Calculating Stress Error\n", adaptIt + 1);
		LOGPRINT("Stress Smoothing and Calculating Stress Error\n");
		if (its != (nits - 1))
			adapted = Adapt(its);
		else
			adapted = Adapt(its, CHECKERRORONLY);//final pass. CheckErroronly is true

		errForOutput = errorChecker.CalculateMaxErrorForOutput();
		OUTPRINT("Estimated Error in Stress Calculation: %6.2lg percent", errForOutput);
		OUTPRINT("Maximum von Mises Stress in Model: %lg\n", stressUnitConversion*stressMax);
		sumFile.PrintLine("MaxError: %lg", errorMax);
		sumFile.PrintLine("MaxVm: %lg", stressUnitConversion*stressMax);

		SRPassData& pd = passData.Get(its);
		pd.err = errForOutput;
		pd.maxp = pbeforeAdapt;
		pd.neq = numEquations;
		pd.maxsvm = stressMax;
		pd.maxsp1 = maxsp1;
		pd.maxCustom = maxCustom;


		prevStressMax = stressMax;

		FreeSmallElementData();

		lastIt++;

		if (!adapted)
			break;
	}

	bool slopekinkWarnNeeded = false;
	if (maxStressElid != -1)
	{
		SRelement* elemAtMax = GetElement(maxStressElid);
		slopekinkWarnNeeded = elemAtMax->checkSlopeKink();
	}


	//clean up memory and disk space that is only needed during adaptive solution:
	CleanUp(PARTIALCLEANUP);

	allocateSmallElementData();
	OUTPRINT("\nAdaptive loop complete");
	LOGPRINT("\nAdaptive loop complete");
	LOGPRINT("Post Processing");
	post.PostProcess();

	bool flattenedWarningNeeded = checkFlattenedElementHighStress();

	input.numNodeFaces.Free();
	input.nodeFaces.Free();

	SRstring line;
	SRfile& frep = repFile;
	frep.Open(SRappendMode);
	if (breakout)
	{
		fileNameTail.Left('_', line);
		frep.PrintLine("Results of StressRefine breakout Analysis from model %s", line.str);
	}
	else
	{
		frep.PrintLine("Results of StressRefine full model solution");
		if (localAdapt)
			frep.PrintLine("Adapted near max stress in model");
	}
	if (useUnits)
	{
		frep.PrintLine("Units of Stress: %s", stressUnitstr.str);
		frep.PrintLine("Units of Length: %s", lengthUnitstr.str);
	}
	else
	{
		frep.PrintLine("\nUnits of Stress are the same as units for Young's Modulus in Nastran MAT1");
	}

	if (detectSacrificialElements)
	{
		if (SingStressCheck())
			frep.PrintLine("Singular");
	}
	if (errorChecker.smallMaxStressDetected)
		frep.PrintLine("SmallStressDetected");
	if (flattenedWarningNeeded)
		frep.PrintLine("Flattened High StressElement");
	if ((errForOutput) > 10.5)//10.5 instead of 10 is for roundoff
		frep.PrintLine("HighError: %lg", errForOutput);
	if (slopekinkWarnNeeded)
		frep.PrintLine("SlopeKinkAtMax");

	frep.Print("MaxVMIts");
	for (int i = 0; i < lastIt; i++)
		frep.Print(",%lg", passData.Get(i).maxsvm);
	frep.PrintReturn();
	frep.Print("MaxSP1Its");
	for (int i = 0; i < lastIt; i++)
		frep.Print(",%lg", passData.Get(i).maxsp1);
	frep.PrintReturn();
	frep.Print("MaxCustomIts");
	for (int i = 0; i < lastIt; i++)
		frep.Print(",%lg", passData.Get(i).maxCustom);
	frep.PrintReturn();

	if (breakout)
		frep.PrintLine("Maximum von Mises Stress in local region of Model: %12.3lg  ", stressUnitConversion*stressMax);
	else
		frep.PrintLine("Maximum von Mises Stress in Model: %12.3lg  ", stressUnitConversion*stressMax);
	frep.PrintLine("Maximum Principal Stress in Model: %12.3lg  ", stressUnitConversion*maxsp1);
	frep.PrintLine("Minimum Principal Stress in Model: %12.3lg  ", stressUnitConversion*minsp2);

	frep.PrintLine("Estimated Error in Stress Calculation: %6.2lg percent", errForOutput);
	double maxPercentYielded = 0.0;
	int nmat = materials.GetNum();
	int numactivemat = 0;
	for (int i = 0; i < nmat; i++)
	{
		SRmaterial* mat = GetMaterial(i);
		if (mat->isActive())
			numactivemat++;
	}
	if (numactivemat > 1 && !allMatsEquivElast)
	{
		for (int i = 0; i < nmat; i++)
		{
			SRmaterial* mat = GetMaterial(i);
			if (mat->isActive())
			{
				frep.PrintLine("For Elements with Material: %s", mat->GetName());
				double vm = mat->GetMaxSvm();
				frep.PrintLine("    Max von Mises Stress: %12.3lg", stressUnitConversion*vm);
				if (mat->isAllowableAssigned())
				{
					double percentYielded;
					percentYielded = mat->GetVolPercentYielded();
					double fs = 1000.0;
					double allow = mat->GetAllowableStress();
					if (vm > SMALL)
						fs = allow / vm;
					if (fs < 1.0)
					{
						frep.PrintLine("    Allowable Stress: %lg Factor of Safety: %5.2lg -- Yielded", stressUnitConversion*allow, fs);
						if (percentYielded > 1.0)
							percentYielded = 0.9999;
						frep.PrintLine("    Percentage of Material volume yielded: %5.2lf", (percentYielded*100.0));
						if (percentYielded > maxPercentYielded)
							maxPercentYielded = percentYielded;
					}
					else
						frep.PrintLine("    Allowable Stress: %lg Factor of Safety: %lg", stressUnitConversion*allow, fs);
				}
			}
		}
	}
	if (breakout)
		frep.PrintLine("Maximum displacement in local region of Model: %6.2lg  ", lengthUnitConversion*dmax.Magnitude());
	else
		frep.PrintLine("Maximum displacement in Model: %lg  ", lengthUnitConversion*dmax.Magnitude());
	frep.PrintLine("\nMaximum Polynomial Order in Model: %d\n", maxPinModel);
	frep.PrintLine("Number of in Equations in Model: %d\n", numEquations);
	frep.PrintReturn();
	if (nits != 1)
	{
		if (breakout)
		{
			frep.PrintLine("Adaptive Analysis     Max Polynomial Order    Estimated Stress Error    Max Von Mises Stress in");
			frep.PrintLine("Iteration                                                               local region of Model");
		}
		else
		{
			frep.PrintLine("Adaptive Analysis     Max Polynomial Order    Estimated Stress Error    Max Von Mises Stress");
			frep.PrintLine("Iteration                                                               in Model");

		}
		for (int its = 0; its < lastIt; its++)
		{
			SRPassData& pd = passData.Get(its);
			frep.PrintLine("         %d                      %d                  %-6.2lg                    %12.3lg", its + 1, pd.maxp, pd.err, stressUnitConversion*pd.maxsvm);
		}
	}

	frep.Close();
	
	OUTPRINT("\n\nFinal result for Model: %s\n", fileNameTail.str);
	PStats();
	OUTPRINT("Total Volume of Model : %lg\n", volume);
	sumFile.PrintLine("Final Result");
	sumFile.PrintLine("Volume: %lg\n", volume);

	OUTPRINT("Estimated Error in Stress Calculation: %6.2lg percent", errForOutput);
	OUTPRINT("Maximum Polynomial Order in Model: %d\n", maxPinModel);
	OUTPRINT("Maximum nodal displacement (at node %d)\n", nodeUidAtMaxDisp);
	OUTPRINT("ux: %lg\nuy: %lg\nuz: %lg\n", dmax.d[0], dmax.d[1], dmax.d[2]);
	OUTPRINT("Maximum von Mises Stress in Model: %lg  at position %lg %lg %lg\n", stressMax, maxStressPos.d[0], maxStressPos.d[1], maxStressPos.d[2]);
	for (int i = 0; i < 6; i++)
		OUTPRINT("Maximum Stress Component %s in Model: %lg\n", stressComp[i], stressMaxComp[i]);
	sumFile.PrintLine("MaxError: %lg\n", errorMax);
	sumFile.PrintLine("MaxVm %lg", stressMax);
	sumFile.PrintLine("max_nodal_disp: %lg %lg %lg\n", dmax.d[0], dmax.d[1], dmax.d[2]);
	sumFile.PrintLine("Maximum Stress Component");
	for (int i = 0; i < 6; i++)
		sumFile.PrintLine("%s %lg\n", stressComp[i], stressMaxComp[i]);

	if (maxFlattened > TINY)
	{
		OUTPRINT("max element flattening for bad mapping: %lg", maxFlattened);
		OUTPRINT("element: %d", maxFlattenedElUId);
	}

	CleanUp();

	sumFile.Close();
}


bool SRmodel::Adapt(int pIteration, bool checkErrorOnly)
{
	//adapt polynomial orders of all elements after error checking;
	//input:
		//pIteration = p-loop iteration number
		//checkErrorOnly = true to check error only but not increase p
	//notes:
		//adaptivity steps:
		//calculate raw strains and perform smoothing (in errorChecker.SetUp)
		//calculate errors (in errorChecker.SetUp):
			//1. smoothed vs raw strains
			//2. traction jumps across shared faces
			//3. traction jumps vs applied loads
		//use errors to determine next p level (in errorChecker.FindRequiredPOrder)

	int i, e;
	SRedge *edge;
	SRelement *elem;
	int p;
	bool pup, anypup = false;
	errorMax = 0.0;

	int maxp0 = maxPinModel;

	bool finalAdapt = (pIteration == adaptLoopMax - 2);

	if (finalAdapt)
		maxPorder = maxPorderFinalAdapt;
	else if (maxPorder > maxPorderFinalAdapt)
		maxPorder = maxPorderFinalAdapt;

	errorChecker.SetUp(finalAdapt);

	if (pIteration == 0 && detectSacrificialElements)
	{
		bool AnySacrElems = errorChecker.AutoSacrificialElements();
		if (AnySacrElems)
		{
			//this has changed sacrificial element status, so redo max calculation:
			post.CalculateMaxStress();
		}
	}

	//check for spurious nodal hot spot at point of max stress:
	//ttd: rethink before turning this on. it spuriously trips for multiclevisassy, and 
	//gives spurious downstream singular stress warning:
	bool doCheckForNodalHotSpot = false;
	if (breakout && doCheckForNodalHotSpot)
	{
		//possible hot spot:
		SRelement* elemAtMax = GetElement(maxStressElid);
		if (!elemAtMax->sacrificial != 0)
		{
			if (elemAtMax->checkForNodalHotSpot())
			{
				//this has changed sacrificial element status, so redo max calculation:
				post.CalculateMaxStress();
			}
			post.nodalStress.Free();
		}
	}

	if (!checkErrorOnly)
	{
		OUTPRINTRET;
		OUTPRINT("Maximum von Mises Stress in Model: %lg\n", stressMax);
		OUTPRINT("Maximum Stress Components in Model\n");
		PrintStresses(stressMaxComp);
	}


	if(pOrderUniform)
	{
		//AdaptUniform will increase p of all edges in model by 1 if any error exceeds tolerance
		anypup = AdaptUniform(pIteration, checkErrorOnly);
		if (checkErrorOnly)
		{
			errorChecker.CleanUp();
			return false;
		}
		for (e = 0; e < elements.GetNum(); e++)
		{
			elem = GetElement(e);
			elem->SetPChanged(true);
		}

		errorChecker.CleanUp();
		return anypup;
	}
	int maxpSacr = 0;
	for (e = 0; e < elements.GetNum(); e++)
	{
		elem = GetElement(e);
		elem->SetPChanged(false);
		pup = errorChecker.FindRequiredPOrder(elem, p);
		if (pup)
		{
			if (!elem->isSacrificial())
			{
				anypup = true;
				if (p > maxPinModel)
					maxPinModel = p;
			}
			else if (p > maxpSacr)
				maxpSacr = p;
			elem->PutNewPorder(p);
			elem->SetPChanged(true);
		}
		else
			elem->PutNewPorder(0);
	}

	if (checkErrorOnly)
	{
		errorChecker.CleanUp();
		maxPinModel = maxp0;
		return false;
	}
	else
	{
		if (anypup && maxpSacr > maxPinModel)
			maxPinModel = maxpSacr;
	}

	if (!anypup)
	{
		errorChecker.CleanUp();
		return anypup;
	}

	for (e = 0; e <elements.GetNum(); e++)
	{
		elem = GetElement(e);
		p = elem->GetNewPorder();
		//update P-order of all the element's edges:
		for (i = 0; i < elem->GetNumLocalEdges(); i++)
		{
			edge = elem->GetEdge(i);
			if (p > edge->GetPorder())
				edge->putPorder(p);
		}
	}

	for (e = 0; e < elements.GetNum(); e++)
	{
		elem = GetElement(e);
		if (elem->GetPChanged())
			continue;
		for (i = 0; i < elem->GetNumLocalEdges(); i++)
		{
			edge = elem->GetEdge(i);
			if (edge->PChanged())
			{
				elem->SetPChanged(true);
				break;
			}
		}
	}

	//for next to last p-pass, after calculating stresses,
	//set high stress elements near breakout boundary
	//to sacrificial
	if (finalAdapt)
		checkForHotSpotElemsNearBreakoutBoundary();

	errorChecker.CleanUp();
	return anypup;
}

bool SRmodel::AdaptUniform(int pIteration, bool checkErrorOnly)
{
	//uniform p adaptivity
	//input:
		//pIteration = p-loop iteration number
		//checkErrorOnly = true to check error only but not increase p
	//notes:
		//finds max error in model. if less than tolerance, return "any-p-increased" = false
		//else:
			//increase p of all edges in model by 1 unless they are already at max,
	//return:
		//"any-p-increased" = true

	int e;
	SRelement* elem;
	SRedge* edge;

	for (e = 0; e < elements.GetNum(); e++)
	{
		elem = GetElement(e);
		if (elem->isSacrificial())
			continue;
		errorChecker.FindError(elem);
	}

	errorMax = errorChecker.CalculateMaxErrorForOutput();
	errorMax *= 0.01;//converting from percentage to actual

	if (checkErrorOnly)
	{
		errorChecker.CleanUp();
		return false;
	}

	if (errorMax < ErrorTolerance)
		return false;

	SRintVector edgeupped;
	int id, i, p;
	edgeupped.Allocate(edges.GetNum());
	for (e = 0; e < elements.GetNum(); e++)
	{
		elem = GetElement(e);
		if (elem->isSacrificial())
			continue;
		for (i = 0; i <elem->GetNumLocalEdges(); i++)
		{
			id = elem->GetLocalEdgeGlobalId(i);
			if (edgeupped.Get(id))
				continue;
			edgeupped.Put(id, 1);
			edge = GetEdge(id);
			p = edge->GetPorder();
			p++;
			if (p <= maxPorder)
			{
				edge->putPorder(p);
				if (p > maxPinModel)
					maxPinModel = p;
			}
		}
	}
	edgeupped.Free();
	return true;
}

void SRmodel::NumberGlobalFunctions()
{
	//assign global function numbers to all functions in model

	SRnode* node;
	SRedge* edge;
	SRface* face;
	SRelement* elem;
	int i, n, nint;
	int fun;
	int j, k, pej;
	int lfun, gfun;

	//allocate space for function numbers for edges, faces, and elements:
	for (i = 0; i < edges.GetNum(); i++)
	{
		edge = GetEdge(i);
		n = (edge->pOrder) + 1;
		edge->AllocateGlobalFunctionNumbers(n);
	}
	for (i = 0; i < faces.GetNum(); i++)
	{
		face = GetFace(i);
		n = basis.CountFaceTotalFunctions(face, j);
		face->AllocateGlobalFunctionNumbers(n);
	}
	for (i = 0; i < elements.GetNum(); i++)
	{
		elem = GetElement(i);
		n = basis.CountElementFunctions(elem, nint);
		elem->AllocateGlobalFunctionNumbers(n);
	}

	//for nodes, assign a function number unless it is a midside node
	//or orphan:
	fun = 0;
	for (i = 0; i < nodes.GetNum(); i++)
	{
		node = GetNode(i);
		if (!node->isMidSide() && !node->isOrphan())
		{
			node->globalFunctionNumber = fun;
			fun++;
		}
	}

	//edges:
	//assign p2 funs 1st for easier mapping of p2 soln to adapted solns (e.g. itsolv):
	for (i = 0; i < edges.GetNum(); i++)
	{
		edge = GetEdge(i);
		edge->AssignGlobalFunctionNumbers(fun, 2, 2);
		//assign the same function number to the corresponding global node.
		//but fun was incremented in edge->AssignGlobalFunctionNumbers so subtract 1:
		SRnode* node = GetNode(edge->GetMidNodeId());
		node->PutGlobalFunctionNumber(fun - 1);
	}
	for (i = 0; i < edges.GetNum(); i++)
	{
		edge = GetEdge(i);
		edge->AssignGlobalFunctionNumbers(fun, 3, 8);
	}

	//faces:
	int nf = faces.GetNum();
	int nn, ne;
	for (i = 0; i < nf; i++)
	{
		lfun = 0;
		face = GetFace(i);
		//corner functions:
		nn = face->GetNumNodes();
		for (j = 0; j < nn; j++)
		{
			gfun = face->GetNode(j)->globalFunctionNumber;
			face->PutGlobalFunctionNumber(lfun, gfun);
			lfun++;
		}

		//edge p2 functions:
		ne = face->GetNumLocalEdges();
		for (j = 0; j < ne; j++)
		{
			edge = face->GetEdge(j);
			gfun = edge->GetGlobalFunctionNumber(2);
			face->PutGlobalFunctionNumber(lfun, gfun);
			lfun++;
		}

		//edge higher functions:
		ne = face->GetNumLocalEdges();
		for (j = 0; j < ne; j++)
		{
			edge = face->GetEdge(j);
			pej = edge->pOrder;
			for (k = 3; k <= pej; k++)
			{
				gfun = edge->GetGlobalFunctionNumber(k);
				face->PutGlobalFunctionNumber(lfun, gfun);
				lfun++;
			}
		}

		n = face->GetNumGlobalFunctions();
		for (j = lfun; j < n; j++)
		{
			face->PutGlobalFunctionNumber(lfun, fun);
			lfun++;
			fun++;
		}
	}
	//elements:
	maxNumElementFunctions = 0;
	for (i = 0; i < elements.GetNum(); i++)
	{
		lfun = 0;
		elem = GetElement(i);
		//nodes:
		n = elem->GetNumNodes();
		for (j = 0; j < n; j++)
		{
			gfun = elem->GetNode(j)->globalFunctionNumber;
			elem->PutGlobalFunctionNumbers(lfun, gfun);
			lfun++;
		}

		//edge p2 functions:
		for (j = 0; j < elem->GetNumLocalEdges(); j++)
		{
			edge = elem->GetEdge(j);
			gfun = edge->GetGlobalFunctionNumber(2);
			elem->PutGlobalFunctionNumbers(lfun, gfun);
			lfun++;
		}

		//edge higher functions:
		for (j = 0; j < elem->GetNumLocalEdges(); j++)
		{
			edge = elem->GetEdge(j);
			pej = edge->pOrder;
			for (k = 3; k <= pej; k++)
			{
				gfun = edge->GetGlobalFunctionNumber(k);
				elem->PutGlobalFunctionNumbers(lfun, gfun);
				lfun++;
			}
		}

		//face functions:
		for (j = 0; j < elem->GetNumLocalFaces(); j++)
		{
			face = elem->GetFace(j);
			n = basis.CountFaceTotalFunctions(face, nint);
			//skip number of edge and node functions on face:
			for (k = (n - nint); k < n; k++)
			{
				gfun = face->GetGlobalFunctionNumber(k);
				elem->PutGlobalFunctionNumbers(lfun, gfun);
				lfun++;
			}
		}

		//internal functions:
		n = basis.CountElementFunctions(elem, nint);
		if (n > maxNumElementFunctions)
			maxNumElementFunctions = n;
		for (j = 0; j < nint; j++)
		{
			elem->PutGlobalFunctionNumbers(lfun, fun);
			lfun++;
			fun++;
		}
	}

	numFunctions = fun;

}

void SRmodel::CalculateElementStiffnesses(bool anyLcsEnfd)
{
	//calculate elemental stiffness matrix for all elements in model and store on disk
	//unless elements fit in memory

	bool okStif = false;
	if (maxNumThreads > 1 && elementsInMemory)
	{
		OUTPRINT("parallel element calculation\n");

		okStif = CalculateElementStiffnessesOmp();
		if (!okStif)
		{
			//parallel element calculations failed. recover by setting num processors to 1
			try
			{

				FreeElementData();
				maxNumThreads = 1;
				allocateElementData();
			}
			catch (...)
			{
				ERROREXIT;
			}
		}
	}
	if (!okStif)
	{

		OUTPRINT("NOT parallel elements\n");

		//scratch space needed by elements:
		allocateElementData();

		int n = elements.GetNum();
		LOGPRINT("Calculating %d elements", n);

		openElStiff();
		int dprogTarget = (int)(0.1* (double)n);
		int progTarget = dprogTarget;
		int progout = 10;
		for (int i = 0; i < n; i++)
		{
			if (i > progTarget)
			{
				progTarget += dprogTarget;
				progout += 10;
			}

			SRelement* elem = GetElement(i);
			int len;
			double *stiff = elem->CalculateStiffnessMatrix(0, len);

			WriteElementStiffness(elem, len, stiff);
		}

		LOGPRINT("\n");

		for (int i = 0; i < n; i++)
		{
			SRelement* elem = GetElement(i);
			volume += elem->elVol;
		}

		closeElStiff();

		FreeElementData();
	}

	if (anyLcsEnfd)
	{
		for (int i = 0; i < maxNumThreads; i++)
		{
			double *enfdForceVec = elData.GetPointer(i)->enfdForceVec.d;
			for (int j = 0; j < numEquations; j++)
				solution.PlusAssign(j, enfdForceVec[j]);
		}
	}
}

bool SRmodel::CalculateElementStiffnessesOmp()
{
	//calculate elemental stiffness matrix for all elements in model using multiple processors
	//note:
	//only call this if there is more than one thread and elements fit in memory

	//scratch space needed by elements:
	try
	{
		allocateElementData();

		int n = elements.GetNum();
		LOGPRINT("Calculating %d elements\n", n);
		int processorNum = 0;

		int elnotProc0 = -1;

#pragma omp parallel
		{

#pragma omp for nowait
			for (int i = 0; i < n; i++)
			{
				processorNum = omp_get_thread_num();
				if (processorNum != 0 && elnotProc0 == -1)
					elnotProc0 = i;

				SRelement* elem = GetElement(i);
				int len;
				elem->CalculateStiffnessMatrix(processorNum, len);
			}
		}

		for (int i = 0; i < n; i++)
		{
			SRelement* elem = GetElement(i);
			volume += elem->elVol;
		}

		FreeElementData();
		return true;
	}
	catch (...)//... means any exception in c++
	{
		SCREENPRINT(" exception during parallel element calculation");
		return false;
	}
}

void SRmodel::Create()
{
	//create model by processing input file
	input.ReadModel();

	for (int i = 0; i < materials.GetNum(); i++)
	{
		SRmaterial* mat = GetMaterial(i);
		if (mat->isActive())
		{
			if(!mat->isAllowableAssigned())
				allMatsHaveAllowable = false;
			else
			{
				if (mat->GetAllowableStress() > maxAllowableAnyActiveMat)
					maxAllowableAnyActiveMat = mat->GetAllowableStress();
			}
		}
	}
}

void SRmodel::AllocateDofVectors()
{
	//allocate space for degree-of-freedom vectors for
	//function equation numbers and enforced displacements

	functionEquations.Allocate(numFunctions, 3);
	enforcedDisplacements.Allocate(numFunctions, 3);
}

void SRmodel::NumberEquations()
{
	//global equation numbering
	//note:
		//fills class variable functionEquations
		//ProcessConstraints has to be called 1st. 
		//it assigns a negative number to global dofs that are constrained
		//this routine will only assign an equation number to unconstrained
		//global dofs


	int eq = 0;
	for (int gfun = 0; gfun < numFunctions; gfun++)
	{
		for (int dof = 0; dof < 3; dof++)
		{
			//unconstrained dofs have temporarily been assigned to 0.
			//constrained dofs are assigned a negative number
			int eqt = GetFunctionEquation(gfun, dof);
			if (eqt == 0)
			{
				PutFunctionEquation(gfun, dof, eq);
				eq++;
			}
		}
	}
	numEquations = eq;
}

void SRmodel::SetStressMax(SRelement* elem, SRvec3& pos, double svm)
{
	//set the value of max stress in model, at corresponding element id and position number
	//input:
		//pos = position of stress that is greater than current max
		//svm = von Mises stress value
	maxStressPos.Copy(pos);
	stressMax = svm;
	maxStressElid = elem->id;
}

void SRmodel::SetStressMaxComp(double *stressComp)
{
	//check for any components of a stress tensor exceed max stress in model
	//input:
		//stressComp = stress tensor stored as vector
	//note:
		//updates class variable stressMaxComp

	for (int i = 0; i < 6; i++)
	{
		if (fabs(stressComp[i]) > fabs(stressMaxComp[i]))
			stressMaxComp[i] = stressComp[i];
	}
	double sp1, sp2;
	math.GetPrinStress(stressComp, sp1, sp2);
	if (sp1 > maxsp1)
		maxsp1 = sp1;
	if (sp2 < minsp2)
		minsp2 = sp2;
}




double SRmodel::GetDisplacementCoeff(int gfun, int dof)
{
	//look up displacement coefficient for a global function and dof
	//input:
		//gfun = global function number
		//dof = degree of freedom
	//return:
		//displacement coefficient

	int eq = GetFunctionEquation(gfun, dof);
	if (eq >= 0)
		return solution.Get(eq);
	else
		return GetEnforcedDisplacement(gfun, dof);
}

double SRmodel::GetEnforcedDisplacementCoeff(int gfun, int dof)
{
	//look up displacement coefficient for a global function and dof
	//input:
		//gfun = global function number
		//dof = degree of freedom
	//return:
		//displacement coefficient

	int eq = GetFunctionEquation(gfun, dof);
	if (eq >= 0)
		return 0.0;
	else
		return GetEnforcedDisplacement(gfun, dof);
}

void SRmodel::NumberEquationsSmooth()
{
	//number global equations for smoothing.
	//notes:
		//there is only one dof.
		//skip functions not belonging to elements in current material group
		//(determined by post.elSmooth flag)
	int i, e, gfun;
	SRelement* elem;
	for(i = 0; i < numFunctions; i++)
		skipFun.Put(i, 1);
	for (e = 0; e <elements.GetNum(); e++)
	{
		elem = GetElement(e);
		if (post.getElSmooth(e))
		{
			for(i = 0; i < elem->GetNumFunctions(); i++)
			{
				gfun = elem->GetFunctionNumber(i);
				skipFun.Put(gfun, 0);
			}
		}
	}


	int eq = 0;
	for(i = 0; i < numFunctions; i++)
	{
		if(skipFun.Get(i) == 0)
		{
			smoothFunctionEquations.Put(i, eq);
			eq++;
		}
		else
			smoothFunctionEquations.Put(i, -1);
	}
	numSmoothEquations = eq;
}


void SRmodel::ProcessConstraints()
{
	//process constraints
	//notes:
		//number functionEquations(fun,dof) as negative for a gcs constrained global dof.
		//the number is -1 if there is no enforced displacement, else
		//it is the number in the enforcedDisplacements matrix where the enforced displacement is stored
		//class variables functionEquations and enforcedDisplacements are updated
		//lcs constraints are skipped because they are handled via penalty:


	int i, dof, gfun;
	SRedge* edge;
	SRconstraint* con;
	SRnode* node;

	AllocateDofVectors();

	for (i = 0; i < constraints.GetNum(); i++)
	{
		con = GetConstraint(i);

		//skip lcs constraint, they are handled via penalty:
		if (!con->isGcs())
			continue;
		int eid = con->GetEntityId();
		if (con->GetType() == nodalCon)
		{
			node = GetNode(eid);
			gfun = node->GetGlobalFunctionNumber();
			for (dof = 0; dof < 3; dof++)
			{
				if (con->IsConstrainedDof(dof))
				{
					PutFunctionEquation(gfun, dof, -1);
					if (con->hasEnforcedDisp())
					{
						double enfd = con->getDisp(0, dof);
						PutEnforcedDisplacement(gfun, dof, enfd);
					}
				}
			}
		}
		else if (con->GetType() == faceCon)
		{
			if (!con->isBreakout())
				con->ProcessFaceConstraint();
		}
	}

	//process breakout constraints last so they override others at conflicting entities
	//(e.g. constraint and breakout constraint share same edge):
	for (i = 0; i < constraints.GetNum(); i++)
	{
		con = GetConstraint(i);
		if (!con->isGcs())
			continue;
		int eid = con->GetEntityId();
		if (con->isBreakout() && (con->GetType() != nodalCon))
		{
			//breakout constraints are always face constraints
			//(except nodal breakout constraints which are handled separately)
			con->ProcessFaceConstraint();
		}
	}

	for (i = 0; i < GetNumElements(); i++)
	{
		SRelement* elem = GetElement(i);
		for (int f = 0; f < elem->GetNumFunctions(); f++)
		{
			int gfun = elem->GetFunctionNumber(f);
			for (dof = 0; dof < 3; dof++)
			{
				if (functionEquations.Get(gfun, dof) < 0)
				{
					elem->SetConstrained();
					break;
				}
			}
			if (elem->isConstrained())
				break;
		}
	}

	int nfun = GetNumFunctions();
	funUncon.Allocate(nfun);
	funUncon.Set(1);
	//funUncon is used for solver bookkeeping. it is true only
	//if the function is not owned by any elements that have constraints, even
	//if they are not on this function
	for (i = 0; i < GetNumElements(); i++)
	{
		SRelement* elem = GetElement(i);
		if (elem->isConstrained())
		{
			for (int f = 0; f < elem->GetNumFunctions(); f++)
			{
				int gfun = elem->GetFunctionNumber(f);
				funUncon.Put(gfun, 0);
			}
		}
	}
}

bool SRmodel::PreProcessPenaltyConstraints()
{
	//see if any penalty constraints are needed; if so, calibrate the penalty constant

	bool anyLcsEnfd = false;

	SRconstraint* con;
	//fill faceLcsCoonstraints for elements with non-gcs constraints on boundary faces:
	for (int i = 0; i < constraints.GetNum(); i++)
	{
		con = GetConstraint(i);
		if (con->GetType() != faceCon)
			continue;
		if (!con->isGcs())
		{
			int f = con->GetEntityId();
			SRface* face = faces.GetPointer(f);
			if (!face->IsBoundaryFace() )
				ERROREXIT; //not boundary face, non-gcs not supported
			int eId = face->GetElementOwner(0);
			SRelement *elem = GetElement(eId);
			int lface = face->GetElementLocalFace(0);
			elem->CalibratePenalty();
			elem->AddfaceLCSConstraint(lface);
			if (con->hasEnforcedDisp())
				anyLcsEnfd = true;

		}
	}

	//fill nodeLcsCoonstraints for elements with non-gcs constraints on nodes:
	for (int i = 0; i < constraints.GetNum(); i++)
	{
		con = GetConstraint(i);
		if (con->GetType() != nodalCon)
			continue;
		if (!con->isGcs() && !con->breakout)
		{
			int n = con->GetEntityId();
			SRnode* node = nodes.GetPointer(n);
			bool found = false;
			for (int f = 0; f < faces.GetNum(); f++)
			{
				SRface* face = GetFace(f);
				if (!face->IsBoundaryFace())
					continue;
				int eid = face->GetElementOwner(0);
				SRelement* elem = GetElement(eid);
				//find local node corresponding to node:
				for (int l = 0; l < elem->numNodesTotal; l++)
				{
					if (elem->getNodeOrMidNodeId(l) == n)
					{
						elem->AddNodeLcSConstraints(l);
						elem->CalibratePenalty();
						found = true;
						break;
					}
				}
				if (found)
					break;
			}
			if (!found)
				ERROREXIT;//lcs constraint not on boundary- not supported
			if (con->hasEnforcedDisp())
				anyLcsEnfd = true;
		}
	}

	//breakout constraints on nodes:
	for (int i = 0; i < constraints.GetNum(); i++)
	{
		con = GetConstraint(i);
		if (!con->breakout || (con->GetType() != nodalCon))
			continue;
		int n = con->GetEntityId();
		SRnode* node = nodes.GetPointer(n);
		int eid = node->firstElementOwner;
		SRelement* elem = GetElement(eid);
		for (int l = 0; l < elem->numNodesTotal; l++)
		{
			if (elem->getNodeOrMidNodeId(l) == n)
			{
				elem->AddNodeLcSConstraints(l);
				elem->CalibratePenalty();
				break;
			}
		}
		anyLcsEnfd = true;
	}

	return anyLcsEnfd;
}

void SRmodel::FillGlobalFaces(bool needGlobalNodeOrderForBreakout)
{
	//fill global faces in the model
	//note:
		//fills class variable faces
		//faces are created with node numbers, elementOwners, elementLocalFace numbers,
		//and localEdges

	int el, lface, n1, n2, n3, n4, gn1, gn2, gn3, gn4, gface;
	SRface* face;
	SRlocalFace* ellocFace;
	int gej, direction;
	SRelement* elem;
	for (el = 0; el < elements.GetNum(); el++)
	{
		elem = GetElement(el);
		if (saveBreakout && !elem->saveForBreakout)
			continue;
		int nlocface = elem->GetNumLocalFaces();
		for (lface = 0; lface < nlocface; lface++)
		{
			elem->GetFaceNodes(lface, n1, n2, n3, n4);
			gface = SRfaceUtil::GlobalFaceMatch(gn1, gn2, gn3, gn4, n1, n2, n3, n4);
			ellocFace = GetElement(el)->GetLocalFace(lface);
			if (gface == -1)
			{
				//face not found, add new global face:
				gface = faces.GetNum();
				face = faces.Add();
				face->Create(gface,n1, n2, n3, n4);
				face->PutElementOwner(0, el);
				face->PutElementLocalFace(0, lface);
				//assign local edges to the new face:
				gej = SRedgeUtil::GlobalEdgeMatch(n1, n2, direction);
				face->AssignLocalEdge(0, gej, direction);
				gej = SRedgeUtil::GlobalEdgeMatch(n2, n3, direction);
				face->AssignLocalEdge(1, gej, direction);
				if (n4 == -1)
				{
					gej = SRedgeUtil::GlobalEdgeMatch(n1, n3, direction);
					face->AssignLocalEdge(2, gej, direction);
				}
				else
				{
					gej = SRedgeUtil::GlobalEdgeMatch(n4, n3, direction);
					face->AssignLocalEdge(2, gej, direction);

					gej = SRedgeUtil::GlobalEdgeMatch(n1, n4, direction);
					face->AssignLocalEdge(3, gej, direction);
				}
				gn1 = 0;
				gn2 = 1;
				gn3 = 2;
				if (n4 == -1)
					gn4 = -1;
				else
					gn4 = 3;
			}
			else
			{
				face = GetFace(gface);
				face->PutElementOwner(1, el);
				face->PutElementLocalFace(1, lface);
			}
			if (!saveBreakout || needGlobalNodeOrderForBreakout)
			{
				ellocFace->PutGlobalNodeOrder(gn1, 0);
				ellocFace->PutGlobalNodeOrder(gn2, 1);
				ellocFace->PutGlobalNodeOrder(gn3, 2);
				if (gn4 != -1)
					ellocFace->PutGlobalNodeOrder(gn4, 3);
			}
			ellocFace->PutGlobalFaceId(gface);
		}
	}

	if (saveBreakout)
		return;
	//assign boundaryFaceId to all edges of boundary faces;
	//overwrite the numnodefaces and nodefaces arrays with boundary faces only
	int nbdryface = 0;
	input.numNodeFaces.Zero();
	for (int f = 0; f < faces.GetNum(); f++)
	{
		SRface* face = GetFace(f);
		if (!face->IsBoundaryFace())
			continue;
		nbdryface++;
		for (int e = 0; e < face->GetNumLocalEdges(); e++)
		{
			SRedge* edge = face->GetEdge(e);
			edge->PutBoundaryFaceId(f,e);
		}
		for (int i = 0; i < face->GetNumNodes(); i++)
		{
			int nid = face->GetNodeId(i);
			int n = input.numNodeFaces.Get(nid);
			input.PutNodeFace(nid, n, f);
			input.numNodeFaces.PlusAssign(nid, 1);
		}
	}
}

void SRmodel::FindElemsAdjacentToBreakout()
{
	int el, lface, n1, n2, n3, n4, gn1, gn2, gn3, gn4, gface;
	SRface* face;
	int gej, direction;
	SRelement* elem;
	bool anyAdjElFound = false;
	for (int i = 0; i < elements.GetNum(); i++)
	{
		elem = GetElement(i);
		//elements that are not part of the breakout model are candidates for adjacent
		if (elem->saveForBreakout)
			continue;
		for (lface = 0; lface < elem->GetNumLocalFaces(); lface++)
		{
			elem->GetFaceNodes(lface, n1, n2, n3, n4);
			if (checkOrphanNode(n1))
				continue;
			else if (checkOrphanNode(n2))
				continue;
			else if (checkOrphanNode(n3))
				continue;
			else if (n4 != -1 && checkOrphanNode(n4))
				continue;
			gface = SRfaceUtil::GlobalFaceMatch(gn1, gn2, gn3, gn4, n1, n2, n3, n4);
			if (gface != -1)
			{
				face = faces.GetPointer(gface);
				if (face->elementOwners[1] == -1)
					face->elementOwners[1] = i;
				else if (face->elementOwners[0] == -1)
					face->elementOwners[0] = i;
			}
		}
	}
}


void SRmodel::CreateElemEdges(SRelement *elem, int nnodes, int inputNodes[])
{
	//Fill contribution of this element to global edges
	//input:
		//elem = pointer to eleme
		//nnodes = number of input nodes for this element
		//inputNodes = list of input nodes for this element
	//note: 
		//this routine also assigns the global edge ids to the element's local edges
		//fills class variable edges
		//new edges created with node numbers (corner and midside), and localEdges


	int ncorner = 4;
	int nej = 6;
	if (elem->GetType() == wedge)
	{
		ncorner = 6;
		nej = 9;
	}
	else if (elem->GetType() == brick)
	{
		ncorner = 8;
		nej = 12;
	}

	for (int lej = 0; lej < nej; lej++)
	{
		int n1, n2;
		elem->GetEdgeNodeIds(lej, n1, n2);
		int mid = inputNodes[lej + ncorner];
		int direction;
		SRedge* edge;
		int gej = SRedgeUtil::GlobalEdgeMatch(n1, n2, direction);
		if (gej == -1)
		{
			//edge not found:
			gej = GetNumEdges();
			edge = edges.Add();
			edge->Create(n1, n2, mid, gej);
			direction = 1;
			SRnode* node = GetNode(mid);
			node->SetAsMidside(gej);
			node->SetFirstElementOwner(elem->GetId());
		}
		elem->AssignLocalEdge(lej, gej, direction);
	}
}

void SRmodel::EnforcedDisplacementAssemble(bool doingPrevSolution)
{
	//assemble contribution of Enforced Displacements to global force vector
	//notes:
		//global force vector is stored in class variable model.solution

	if (!IsAnyEnforcedDisplacement())
		return;


	double* elstiff = NULL;
	SRfile* elfile = GetElementStiffnessFile();
	if (!elementsInMemory)
		elfile->Open(SRinbinaryMode);
	for (int el = 0; el < GetNumElements(); el++)
	{
		bool stiffhasbeenread = false;
		SRelement* elem = GetElement(el);
		for (int collfun = 0; collfun < elem->GetNumFunctions(); collfun++)
		{
			int colgfun = elem->GetFunctionNumber(collfun);
			for (int coldof = 0; coldof < 3; coldof++)
			{
				int eq = GetFunctionEquation(colgfun, coldof);
				//equation numbers less than zero are a flag that this
				//dof is constrained:
				if (eq < 0)
				{
					double enfdisp = GetEnforcedDisplacement(colgfun, coldof);
					if (!stiffhasbeenread)
					{
						if (doingPrevSolution)
						{
							int len;
							elstiff = elem->CalculateStiffnessMatrix(0, len);
						}
						else
							elstiff = ReadElementStiffness(elem);
						stiffhasbeenread = true;
					}
					//this column of element stiffness contributes to the global force vector:
					int colleq = collfun * 3 + coldof;
					for (int rowlfun = 0; rowlfun < elem->GetNumFunctions(); rowlfun++)
					{
						int rowgfun = elem->GetFunctionNumber(rowlfun);
						for (int rowdof = 0; rowdof < 3; rowdof++)
						{
							int roweq = GetFunctionEquation(rowgfun, rowdof);
							if (roweq >= 0)
							{
								int rowleq = rowlfun * 3 + rowdof;
								int stiffloc = elem->GetStiffnessLocation(rowleq, colleq);
								double enfdForceVal = -enfdisp*elstiff[stiffloc];
								AddToSolutionVector(roweq, enfdForceVal);
							}
						}
					}
				}
			}
		}
	}
	if (!elementsInMemory)
		elfile->Close();
}

void SRmodel::ProcessForces()
{
	//Process all forces in model: node, edge, face, volume, and thermal
	//and assemble into global force vector
	//notes:
		//global force vector is stored in class variable model.solution

	int nf, eid, dof, fun, eq;
	SRforce* force;
	double* globalForce = GetSolutionVector();
	SRvec3 resF;
	nf = forces.GetNum();
	for (int f = 0; f < nf; f++)
	{
		force = forces.GetPointer(f);
		eid = force->GetEntityId();
		if (force->GetType() == nodalForce)
		{
			SRnode* node = GetNode(eid);
			if (node->isOrphan() || node->checkUnsup())
				continue;
			if (!node->isMidSide())
				fun = node->globalFunctionNumber;
			else
			{
				int eid = node->midSideEdgeOwner;
				SRedge* edge = GetEdge(eid);
				fun = edge->globalFunctionNumbers.d[2];
			}
			double fv[3];
			for (dof = 0; dof < 3; dof++)
				fv[dof] = force->GetForceVal(0, dof);
			bool needRotate = false;
			SRmat33 R;
			double rf, sf;
			if (force->GetCoordId() != -1)
			{
				needRotate = true;
				SRcoord* coord = GetCoord(force->GetCoordId());
				coord->GetRotationMatrix(false, node->Position(), R);
			}

			if (needRotate)
			{
				SRvec3 t = fv;
				for (int i = 0; i < 3; i++)
				{
					fv[i] = 0.0;
					for (int j = 0; j < 3; j++)
						fv[i] += (R.rows[i].d[j] * t.d[j]);
				}
			}

			for (dof = 0; dof < 3; dof++)
			{
				eq = GetFunctionEquation(fun, dof);
				if (eq >= 0)
				{
					globalForce[eq] += (fv[dof]);
					resF.d[dof] += fv[dof];
				}
			}
		}
		else if (force->GetType() == edgeForce)
		{
			SRedge* edge = GetEdge(eid);
			edge->ProcessForce(force, resF);
		}
		else if (force->GetType() == faceForce)
		{
			SRface* face = GetFace(eid);
			face->ProcessForce(force, resF);
		}
	}

	if (thermalForce != NULL)
		thermalForce->Process();

	ProcessVolumeForces(resF);
	OUTPRINT("\nResultant load on model:\n");
	OUTPRINT("   Fx: %lg\n", resF.d[0]);
	OUTPRINT("   Fy: %lg\n", resF.d[1]);
	OUTPRINT("   Fz: %lg\n", resF.d[2]);
}


void SRmodel::ProcessVolumeForces(SRvec3& ResF)
{
	//fill up contribution of volume forces to global force vector
	//input:
		//ResF = resultant force vector on model
	//output:
		//ResF = updated
	//note:
		//global force vector is stored in class variable model.solution
		//volume force applies to entire model is it's element list elList is empty,
		//otherwise it only applies to the elements on the list.

	SRelement* elem;
	SRdoubleVector basvec;
	basvec.Allocate(maxNumElementFunctions);
	double* globalForce;
	globalForce = GetSolutionVector();
	SRvolumeForce* vf;
	int nvf = volumeForces.GetNum();
	for (int v = 0; v < nvf; v++)
	{
		vf = volumeForces.GetPointer(v);
		int ellistLen = vf->elList.GetNum();
		if (ellistLen == 0)
		{
			for (int i = 0; i < GetNumElements(); i++)
			{
				elem = GetElement(i);
				ProcessElemVolumeForce(vf, elem, ResF, basvec, globalForce);
			}
		}
		else
		{
			for (int i = 0; i < ellistLen; i++)
			{
				int eid = vf->elList.Get(i);
				elem = GetElement(eid);
				ProcessElemVolumeForce(vf, elem, ResF, basvec, globalForce);
			}
		}
	}
}

void SRmodel::ProcessElemVolumeForce(SRvolumeForce* vf, SRelement*elem, SRvec3& ResF, SRdoubleVector& basvec, double* globalForce)
{
	int nint, gp, fun, dof, eq, nfun, gfun;
	double r, s, t, w, detj;
	double force[3], bw;
	SRvec3 p;

	elem->PutThread(0);
	if (vf->GetType() == gravity)
		vf->GetForceValue(elem, p, force);
	nint = math.FillGaussPoints(elem);
	nfun = elem->GetNumFunctions();
	for (gp = 0; gp < nint; gp++)
	{
		math.GetGP3d(gp, r, s, t, w);
		detj = elem->FillMapping(r, s, t, true);
		w *= detj;
		basis.ElementBasisFuncs(r, s, t, elem, basvec.GetVector());
		if (vf->GetType() == centrifugal)
		{
			elem->Position(r, s, t, p);
			vf->GetForceValue(elem, p, force);
		}

		//ResF += force*w:
		p.Assign(force);
		p.Scale(w);
		ResF.PlusAssign(p);
		SRvec3 pos;
		elem->Position(r, s, t, pos);

		for (fun = 0; fun < nfun; fun++)
		{
			gfun = elem->GetFunctionNumber(fun);
			bw = basvec.Get(fun) * w;
			for (dof = 0; dof < 3; dof++)
			{
				eq = GetFunctionEquation(gfun, dof);
				if (eq >= 0)
					globalForce[eq] += (bw*force[dof]);
			}
		}
	}
}


void SRmodel::CleanUp(bool partial)
{
	//miscelaneous memory and disk cleanup at end of run

	FreeElementData();

	elementFileEven.Delete();
	elementFileOdd.Delete();

	if (partial)
		return;

	nodes.Free();
	edges.Free();
	faces.Free();
	constraints.Free();
	Coords.Free();
	materials.Free();
	forces.Free();
	for (int e = 0; e < GetNumElements(); e++)
	{
		SRelement *elem = GetElement(e);
		elem->Cleanup();
	}
	elements.Free();
	volumeForces.Free();
	if (thermalForce != NULL)
	{
		DELETEMEMORY thermalForce;
		thermalForce = NULL;
	}

	solution.Free();
	functionEquations.Free();
	smoothFunctionEquations.Free();
	enforcedDisplacements.Free();
	skipFun.Free();
	solver.Cleanup();
	post.Cleanup();
	input.Cleanup();
}



void SRmodel::allocateElementData()
{
	//allocate memory intensive data related to basis functions
	//for use by element stiffness routines
	//for multi-thread process there is one copy per processor

	//Note:
		//allocate and these data before element stiffness loop


	//worst case size scratch space for elements
	int nintMax = 0;
	for (int i = 0; i < elements.GetNum(); i++)
	{
		SRelement* elem = GetElement(i);
		int nint = math.CountGaussPoints(elem);
		if (nint > nintMax)
			nintMax = nint;
	}
	int nfunMax = maxNumElementFunctions;
	int neq = 3 * nfunMax;
	int stiffnessLength = (neq)*(neq + 1) / 2;


	for (int i = 0; i < maxNumThreads; i++)
	{
		SRElementData* eld = elData.GetPointer(i);
		eld->dbdx.Allocate(nintMax, nfunMax);
		eld->dbdy.Allocate(nintMax, nfunMax);
		eld->dbdz.Allocate(nintMax, nfunMax);
		eld->elementStiffness.Allocate(stiffnessLength);
	}
}

void SRmodel::FreeElementData()
{
	//free memory intensive data related to basis functions
	//for use by element stiffness routines
	//for multi-thread process there is one copy per processor

	//Note:
		//free these data after element stiffness loop
	if (elData.isEmpty())
		return;
	for (int i = 0; i < maxNumThreads; i++)
	{
		SRElementData* eld = elData.GetPointer(i);
		eld->dbdx.Free();
		eld->dbdy.Free();
		eld->dbdz.Free();
	}
}

void SRmodel::allocateSmallElementData(bool anyLcsEnfd)
{
	//allocate less memory intensive data 
	//for use by element stiffness and stress routines
	//for multi-thread process there is one copy per processor

	//note:
		//these data must persist through element stiffness, error checking, and postprocessing,
		//so allocate top of adaptivity loop.
		//but must be allocated after call to NumberGlobalFunctions so nfunmax is known
		//data must also be allocated before final postprocessing outside adapt loop

	//worst case size scratch space for elements
	int nintMax = 0;
	for (int i = 0; i < elements.GetNum(); i++)
	{
		SRelement* elem = GetElement(i);
		int nint = math.CountGaussPoints(elem);
		if (nint > nintMax)
			nintMax = nint;
	}
	int nfunMax = maxNumElementFunctions;

	elData.Allocate(maxNumThreads);
	for (int i = 0; i < maxNumThreads; i++)
	{
		SRElementData* eld = elData.GetPointer(i);
		eld->basisVec.Allocate(nfunMax);
		eld->dbasisdr.Allocate(nfunMax);
		eld->dbasisds.Allocate(nfunMax);
		eld->dbasisdt.Allocate(nfunMax);
		eld->intWt.Allocate(nintMax);
		if (anyLcsEnfd)
			eld->enfdForceVec.Allocate(numEquations);
	}
}

void SRmodel::FreeSmallElementData()
{
	//note:
		//these data must persist through element stiffness, error checking, and postprocessing,
		//so free at bottom of adaptivity loop.
		//also free after final postprocessing outside adapt loop
	if (elData.isEmpty())
		return;
	for (int i = 0; i < maxNumThreads; i++)
	{
		SRElementData* eld = elData.GetPointer(i);
		eld->basisVec.Free();
		eld->dbasisdr.Free();
		eld->dbasisds.Free();
		eld->dbasisdt.Free();
		eld->intWt.Free();
		eld->enfdForceVec.Free();
	}
	elData.Free();
}


void SRmodel::NodalMaxDisp(int& nodeUidAtMaxDisp, SRvec3& umax)
{
	//find max nodal displacement in model
	//output:
		//nodeUidAtMaxDisp = node Uid at which max disp occurred
	//output:
		//umax =3 dof displacement vector with highest magnitude in model
	double umagmax = 0.0;
	umax.Zero();
	SRvec3 u;
	for (int i = 0; i < nodes.GetNum(); i++)
	{
		SRnode* node = nodes.GetPointer(i);
		node->GetDisp(u);
		double umag = u.Length();
		if (umag > umagmax)
		{
			umagmax = umag;
			umax.Copy(u);
			nodeUidAtMaxDisp = node->GetUserid();
		}
	}
}

void SRmodel::PrintLoadVector()
{
	//print the load vector to the model output file
	OUTOPEN();
	OUTPRINTNORET("resultant applied load vector\n  function    load (dof 1,2,3)\n");
	int nfun = GetNumFunctions();
	for (int fun = 0; fun < nfun; fun++)
	{
		OUTPRINTNORET("%d ", fun);
		for (int dof = 0; dof < 3; dof++)
		{
			int eq = GetFunctionEquation(fun, dof);
			double force = 0.0;
			if (eq >= 0)
				force = solution.Get(eq);
			OUTPRINTNORET("     %lg ", force);
		}
		OUTPRINTNORET("\n");
	}
	OUTPRINTNORET("\n");
	OUTCLOSE();
}

void SRmodel::setMinPorder(int minp)
{
	//set model minumum p-order to minp
	//input:
		//minp = minumum p-order
	if (minp < 2)
		OUTPRINT(" Warning- min polynomial order in model should be at least 2 for accurate solution\n");
	minPorder = maxPinModel = minp;
}


void SRmodel::PrintStresses(double *stress)
{
	//print stress tensor to model the model output file
	//input:
		//stress = stress tensor stored as vector
	for (int i = 0; i < 6; i++)
		OUTPRINT("%s: %lg\n", stressComp[i], stress[i]);
}

void SRmodel::CheckElementsFitInMemory()
{
	//check if all element stiffness matrices will fit in memory
	//note:
		//sets class variable elementsInMemory
	double elStiflen = 0.0;
	int numel = GetNumElements();
	for (int e = 0; e < numel; e++)
	{
		SRelement *elem = GetElement(e);
		int nfun = elem->GetNumFunctions();
		int neq = 3 * nfun;
		int len = neq*(neq + 1) / 2;
		elStiflen += (double) len;
	}
	elStiflen *= sizeof(double);
	if (elStiflen > MaxElementMem)
	{
		OUTPRINT(" writing elements to disk\n");
		elementsInMemory = false;
	}
	else
		elementsInMemory = true;

	double gigFactor = 1024.0*1024.0*1024.0;
	double elGigs = elStiflen / gigFactor;
	OUTPRINT("Storage for all element stiffnesses: %lg GBytes\n", elGigs);
}

void SRmodel::openElStiff()
{
	//open element stiffness file if elements don't fit in memory
	if (!elementsInMemory)
	{
		//open file for writing:
		SRfile *f = GetElementStiffnessFile();
		f->Open(SRoutbinaryMode);
		//open previous-P file for reading unchanged elements
		if (adaptIt != 0)
		{
			f = GetElementStiffnessFile(READPREVP);
			f->Open(SRinbinaryMode);
		}
	}
}


void SRmodel::closeElStiff()
{
	//close element stiffness file if elements don't fit in memory
	if (!elementsInMemory)
	{
		SRfile *f = GetElementStiffnessFile();
		f->Close();
		f = GetElementStiffnessFile(READPREVP);
		f->Close();
	}
}

double* SRmodel::ReadElementStiffness(SRelement* elem, bool readPrevP, int processorNum)
{
	//read element stiffness vector from disk if elements don't fit in memory,
	//else look it up
	//input:
		//elem = pointer to element
		//readPrevP if unchanged element from previous p-pass is to be read, else false
		//processorNum = processor number, 0 if multi-threaded
	//return:
		//start of element stiffness matrix for this element
	double *stiff = NULL;
	if (elementsInMemory)
	{
		int id = elem->GetId();
		stiff = elem->GetStiffnessMatrix();
	}
	else
	{
		stiff = GetElementStiffnessVector(processorNum);
		SRfile* f = GetElementStiffnessFile(readPrevP);
		//if readPrevP and this is 1st adaptive it, f will be NULL
		if (f == NULL)
			return NULL;
		int fpos = elem->GetFilePos();
		f->SeekBinary(fpos);
		f->ReadBinary(elem->GetStiffLength(), stiff);
	}
	return stiff;
}

void SRmodel::WriteElementStiffness(SRelement* elem, int len, double *stiff)
{
	//write element stiffness matrix for an element to a file (model.elementFile)
	//if elements don't fit in memory,
	//input:
		//elem = pointer to element
		//len = length of element stiffness matrix
		//stiff = element stiffness matrix stored as a vector

	//no action required if in memory 
	if (elementsInMemory)
		return;
	SRfile* f = GetElementStiffnessFile();
	elem->SetFilePos(f->GetFilePos());
	f->WriteBinary(len, stiff);
}

void SRmodel::PStats()
{
	//print p-order statistics to model output file
	OUTPRINT("P-order    percentage of model");
	int nedge = edges.GetNum();
	for (int p = minPorder; p <= maxPinModel; p++)
	{
		int edgesAtP = 0;
		for (int i = 0; i < nedge; i++)
		{
			if (GetEdge(i)->GetPorder() == p)
				edgesAtP++;
		}
		double percentage = 100.0* ((double)edgesAtP) / ((double)nedge);
		OUTPRINT("      %d    %4.2lf", p, percentage);
	}
}

void SRmodel::SetErrorMax(double error, int eluid, double errorSmoothRaw, double errorFaceJumps)
{
	//set the maximum error in model
	//input:
		//error = error value in an element
		// eluid = element userid
		//errorSmoothRaw = smooth vs raw error value in the element
		//errorFaceJumps = face traction jump error value in the element
	//note:
	//if error is greater than current max, sets class variables errorMax, errorSmoothRawAtMax, errorFaceJumpAtMax, and elUidAtMaxError
	if (error > errorMax)
	{
		errorMax = error;
		errorSmoothRawAtMax = errorSmoothRaw;
		errorFaceJumpAtMax = errorFaceJumps;
		elUidAtMaxError = eluid;
	};
}

void SRmodel::checkElementMapping()
{

	//check mapping of all elements at all gauss points needed this p-pass.
	//straighten edges and mark associated elements as sacrificial if mapping is bad.

	int nel = GetNumElements();
	bool anyFail = false;
	int e;
	for (e = 0; e < nel; e++)
	{
		SRelement* elem = GetElement(e);
		elem->FillMappingNodes();
		elem->PutThread(0);
		if (!elem->testMapping())
			anyFail = true;
	}

	if (!anyFail)
		return;

	for (e = 0; e < nel; e++)
	{
		SRelement* elem = GetElement(e);
		elem->checkForStraightenedEdges();
	}

	for (e = 0; e < nel; e++)
	{
		SRelement* elem = GetElement(e);
		if (!elem->flattened)
			continue;
		elem->FillMappingNodes();
		elem->thread = 0;
		if (!elem->testMapping())
		{
			//OUTPRINT("checkElementMapping elem failed: %d", elem->userId);
			anyFail = true;
		}
	}
	//flattening of an element may have invalidated an adjacent element
	if (!anyFail)
		return;

	//flattening of an element may have invalidated an adjacent element
	for (e = 0; e < nel; e++)
	{
		SRelement* elem = GetElement(e);
		if (!elem->flattened)
			continue;
		elem->FillMappingNodes();
		elem->thread = 0;
		if (!elem->testMapping())
			anyFail = true;
	}
	if (anyFail)
	{
		partialFlatten = false;
		for (e = 0; e < nel; e++)
		{
			SRelement* elem = GetElement(e);
			if (!elem->flattened)
				continue;
			elem->FillMappingNodes();
			elem->thread = 0;
			if (!elem->testMapping())
				anyFail = true;
		}
		for (e = 0; e < nel; e++)
		{
			SRelement* elem = GetElement(e);
			elem->checkForStraightenedEdges();
		}
	}
}

bool SRmodel::checkFlattenedElementHighStress()
{
	//check if an element whose mapping was flattened is adjacent to max- possible "error pollution"
	bool adjacentElemFlattened = false;

	double scaledFlattenTol = 0.1; // this is like flattening an edge that spans 45 degrees of a circle

	//rescale edge straightenfraction using edge chord length, then update elements flattenfraction
	for (int e = 0; e < GetNumEdges(); e++)
	{
		SRedge* edge = GetEdge(e);
		edge->ScaleStraightenVsEdgeLength();
	}
	for (int e = 0; e < GetNumElements(); e++)
	{
		SRelement* elem = GetElement(e);
		elem->rescaleFlattenFraction();
	}

	if (maxStressElid == -1)
		return false;
	SRelement* elemAtMax = GetElement(maxStressElid);
	double maxFlattenHighStress = 0.0;
	SRintVector badAdjElems;
	for (int n = 0; n < elemAtMax->GetNumNodesTotal(); n++)
	{
		int nid = elemAtMax->getNodeOrMidNodeId(n);
		int nf = input.numNodeFaces.Get(nid);
		for (int f = 0; f < nf; f++)
		{
			int fid = input.GetNodeFace(nid, f);
			SRface* face = GetFace(fid);
			for (int e = 0; e < 2; e++)
			{
				int elid = face->elementOwners[e];
				if (elid == -1 || elid == maxStressElid)
					continue;
				SRelement* adjElem = GetElement(elid);
				bool bf = false;
				for (int f = 0; f < adjElem->GetNumLocalFaces(); f++)
				{
					SRface* face = adjElem->GetFace(f);
					if (face->IsBoundaryFace())
					{
						bf = true;
						break;
					}
				}
				if (!bf)
					continue;

				if (adjElem->flattenFraction > scaledFlattenTol)
				{
					badAdjElems.PushBack(elid);
					adjacentElemFlattened = true;
					if (adjElem->flattenFraction > maxFlattenHighStress)
					{
						maxFlattenHighStress = adjElem->flattenFraction;
					}
				}
			}
		}
	}
	if (adjacentElemFlattened)
	{
		OUTPRINT(" max flattening adjacent to max stress element: %lg", maxFlattenHighStress);
		//ttd: comment out after dbg
		//plot faces of elements that were flattend more than scaledFlattenTol, but only if they have bdry faces:
		//post.PlotElems(badAdjElems.GetNum(), badAdjElems.d);
		return true;
	}
	return false;
}


SRElementData* SRmodel::GetElementData(int i)
{
	//look up the element data for a processor number
	//input:
		//i = processor number, ignored if not multi-threaded
	if (maxNumThreads == 1)
		return elData.GetPointer(0);
	else
		return elData.GetPointer(i);
};

SRfile* SRmodel::GetElementStiffnessFile(bool readPrevP)
{
	//return reference to element stiffness file for reading or writing
	//input:
		// readPrevP = true for read stiffness from previous p-pass, else false
	bool odd = math.Odd(adaptIt);
	if (readPrevP)
		odd = !odd;
	if (odd)
		return &elementFileOdd;
	else
		return &elementFileEven;
}

double* SRmodel::GetElementStiffnessVector(int processorNum)
{
	//look up the element stiffness Vector for a processor number
	//input:
		//processorNum = processor number, 0 if not multi-threaded

	SRElementData* eld = GetElementData(processorNum);
	return eld->elementStiffness.GetVector();
}

void SRmodel::WriteAllElements()
{
	//write all elements to disk if they are currently in memory and the memory is needed for the solver

	OUTPRINT(" solver needes more memory. writing elements to disk\n");
	//solver needs more memory and elements are in memory, so write elements to disk:
	elementsInMemory = false;
	SRfile* f = GetElementStiffnessFile();
	f->Open(SRoutbinaryMode);
	for (int e = 0; e < GetNumElements(); e++)
	{
		SRelement* elem = GetElement(e);
		double* stiff = elem->GetStiffnessMatrix();
		int len = elem->GetStiffLength();
		WriteElementStiffness(elem, len, stiff);
		elem->deleteStiffnessMaxtrix();
	}
	f->Close();
	MaxElementMem = 0.0; //this will force OOC elements henceforth
}

void SRmodel::dumpFileOpen(bool append)
{
	SRstring name;
	inputFile.filename.Left('.', name);
	name.Cat("_dumpFile.txt");
	if (append)
		dumpFile.Open(name, SRappendMode);
	else
		dumpFile.Open(name, SRoutputMode);
}


void SRmodel::SetEdgesToPorder(int p)
{
	SRedge* edge;
	for (int i = 0; i < edges.GetNum(); i++)
	{
		edge = GetEdge(i);
		edge->putPorder(p);
	}
}

void SRmodel::CreateElemEdgesSBO(SRelement *elem)
{
	//Fill contribution of this element to global edges
	//version for save breakout only- elem was previously created
	//input:
		//elem = pointer to element
	//this routine also assigns the global edge ids to the element's local edges
	//fills class variable edges
	//new edges created with node numbers (corner and midside), and localEdges

	int ncorner = 4;
	int nej = 6;
	if (elem->GetType() == wedge)
	{
		ncorner = 6;
		nej = 9;
	}
	else if (elem->GetType() == brick)
	{
		ncorner = 8;
		nej = 12;
	}

	for (int lej = 0; lej < nej; lej++)
	{
		int n1, n2;
		elem->GetEdgeNodeIds(lej, n1, n2);
		int mid = elem->nodeIds.Get(lej + ncorner);
		int direction;
		SRedge* edge;
		int gej = SRedgeUtil::GlobalEdgeMatch(n1, n2, direction);
		if (gej == -1)
		{
			//edge not found:
			gej = GetNumEdges();
			edge = edges.Add();
			edge->Create(n1, n2, mid, gej);
			direction = 1;
			SRnode* node = GetNode(mid);
			node->SetAsMidside(gej);
			node->SetFirstElementOwner(elem->GetId());
		}
		elem->AssignLocalEdge(lej, gej, direction);
	}
}

void SRmodel::GetFileNameFromExtension(char* ext, SRstring& name)
{
	outputFile.filename.Left('.', name);
	name += ".";
	name += ext;
}

SRnode* SRmodel::GetNodeFromUid(int i)
{
	int id = input.NodeFind(i);
	if (id == -1)
		return NULL;
	else
		return nodes.GetPointer(id);
}

bool SRmodel::SingStressCheck()
{
	//return true if distance from location of max stress to any sacrficial element is small compared to
	//size of element at max
	double distTol = BIG;
	if (maxStressElid == -1)
		return false;
	SRelement* elem = GetElement(maxStressElid);
	distTol = elem->size;
	if (breakout)
		distTol = MATHMIN(distTol, post.maxRadNearOrigin);
	double sacrDistMin = BIG;
	int elAtMinDist;
	for (int e = 0; e < GetNumElements(); e++)
	{
		SRelement* elem = GetElement(e);
		if (!elem->isSacrificial() || elem->flattened)
			continue;
		SRvec3 elcent;
		map.ElementCentroid(elem, elcent);
		double dist = maxStressPos.Distance(elcent);
		if (dist < sacrDistMin)
		{
			sacrDistMin = dist;
			elAtMinDist = e;
		}
	}
	if (sacrDistMin < distTol)
	{
		int maxStressEluid = GetElement(maxStressElid)->userId;
		int elUidAtmindist = GetElement(elAtMinDist)->userId;
		return true;
	}
	else
		return false;
}

void SRmodel::AddSoftSprings()
{
	//add soft springs to each node for model stabilization.
	//approx model average stiffness: K = AE/L. A ~ L*L
	//L = model size. so K = LE. apply ~1e-08 times this to each element
	//face, so 1/3 that to each node
	double kmult = 1.e-08;
	double L = size;
	double E = GetMaterial(0)->E;
	double nodestiff = ONETHIRD*kmult*L*E;
	SRintVector nodeDone;
	int nnodes = GetNumNodes();
	nodeDone.Allocate(nnodes);
	int nel = GetNumElements();
	for (int elnum = 0; elnum < nel; elnum++)
	{
		SRelement* elem = GetElement(elnum);
		if (!elem->pChanged)
			continue;
		double *elStiff = ReadElementStiffness(elem);
		int nn = elem->GetNumCorners();
		for (int n = 0; n < nn; n++)
		{
			int nid = elem->GetNodeId(n);
			if (nodeDone.d[nid] == 1)
				continue;
			nodeDone.d[nid] = 1;
			int lfun = n;
			for (int dof = 0; dof < 3; dof++)
			{
				int eq = lfun * 3 + dof;
				int rowcolloc = elem->stiffDiag.Get(eq);
				elStiff[rowcolloc] += nodestiff;
			}
		}
	}
}

void SRmodel::readResultsSrr()
{
	SRfile f;
	SRstring line, tok;
	inputFile.filename.Left('.', line);
	line += ".srr";
	if (!f.Open(line, SRinputMode))
	{
		REPPRINT("Error: Nastran displacement results not available or incomplete");
		ERROREXIT;
	}
	f.GetLine(line);
	if (!line.CompareUseLength("displacements"))
	{
		REPPRINT("Error: Nastran displacement results not available or incomplete");
		ERROREXIT;
	}
	int nuid;
	SRvec3 disp;
	int nnode = nodes.GetNum();
	post.nodeDisps.Allocate(nnode);


	int numDispsRead = 0;
	SRstring linesav;
	SRnode* node;

	int iStart = 0;
	int nid;

	bool eidmaxLineRead = false;

	while (1)
	{
		if (!f.GetLine(line))
			break;
		if (line.CompareUseLength("eidAtMax"))
		{
			eidmaxLineRead = true;
			break;
		}
		line.setTokSep(",");
		linesav = line;
		if (!line.TokRead(nuid))
			break;
		if (!line.TokRead(disp.d[0]))
			break;
		if (!line.TokRead(disp.d[1]))
			break;
		if (!line.TokRead(disp.d[2]))
			break;
		int nid = input.NodeFind(nuid);
		if (nid != -1)
		{
			post.nodeDisps.Put(nid, disp);
			node = GetNode(nid);
			node->hasDisp = true;
			numDispsRead++;
		}
	}

	//sanity check: if not fromPartialDispFile, all non-orphan nodes have to have disps:
	if (!saveBreakout && !saveBreakoutData.fromPartialDispFile)
	{
		if (anyUnsupNode || saveBreakout)
		{
			for (int n = 0; n < nnode; n++)
			{
				SRnode* node = GetNode(n);
				if (!node->isOrphan() && !node->hasDisp)
				{
					REPPRINT("Error: Nastran Displacement results file incompatible with model file");
					ERROREXIT;
				}
			}
		}
	}


	feasvmmax = 0.0;

	if (saveBreakout && saveBreakoutData.atMax && !saveBreakoutData.fromPartialDispFile)
	{
		if (!eidmaxLineRead)
			f.GetLine(line);
		if (!line.CompareUseLength("eidAtMax"))
		{
			//fatal error if stresses are needed:
			REPPRINT("Error: Nastran stress results not available. Cannot solve local region at maximum stress in model");
			ERROREXIT;
		}

		int eluidatmax = -1;
		double svmmax;
		line.Right(':', tok);
		tok.TokRead(eluidatmax);
		tok.TokRead(svmmax);
		int eid = input.elemFind(eluidatmax);
		if (eid == -1)
		{
			REPPRINT("Error: Nastran stress results not available. Cannot solve local region at maximum stress in model");
			ERROREXIT;
		}
		SRelement* elem = GetElement(eid);
		elem->approxCentroid(feavmmaxpos);
		saveBreakoutData.origin.Copy(feavmmaxpos);
	}
	f.Close();
}

double SRmodel::findMaxStressNearBreakoutOrigin()
{
	double svmmaxnearorgin = 0.0;
	for (int i = 0; i < GetNumElements(); i++)
	{
		SRelement* elem = GetElement(i);
		if (elem->isSacrificial())
			continue;
		double rad = elem->centroidDistance(saveBreakoutData.origin);
		if (rad < post.maxRadNearOrigin)
		{
			if (elem->svmMax > svmmaxnearorgin)
				svmmaxnearorgin = elem->svmMax;
		}
	}
	return svmmaxnearorgin;
}

void SRmodel::setAdaptLoopMax(int nits)
{
	adaptLoopMax = nits;
	if (nits != 3 && maxPJump == 10)
		maxPJump = 2;
}

void SRmodel::setElemsAllNodesHaveDisp()
{
	int nelall = 0;
	if (saveBreakoutData.fromPartialDispFile)
	{
		for (int e = 0; e < GetNumElements(); e++)
		{
			SRelement* elem = GetElement(e);
			elem->checkAllNodesHaveDisp();
			if (elem->allNodesHaveDisp)
				nelall++;
		}
		if (nelall < 1000)
		{
			REPPRINT("Warning: Very few elements were found in the local region of your model");
			REPPRINT("Please review your node / displacement data table");
			ERROREXIT;
		}

	}
	else
	{
		int nel = GetNumElements();
		for (int e = 0; e < nel; e++)
		{
			SRelement* elem = GetElement(e);
			elem->allNodesHaveDisp = true;
		}
		nelall = nel;
	}
	post.nelAllNodesDisp = nelall;
}

void SRmodel::checkForHotSpotElemsNearBreakoutBoundary()
{
	for (int e = 0; e < GetNumElements(); e++)
	{
		SRelement* elem = GetElement(e);
		double elsvmmax = elem->GetSvmMax();
		if (elsvmmax > 0.8*stressMax)
		{
			if (elem->checkNearBktBdry())
				elem->sacrificial = 1;
		}
	}
}

void SRmodel::SetNumThreads()
{
#pragma omp parallel
#pragma omp master
	{
		maxNumThreads = omp_get_num_threads();
		SCREENPRINT(" Using %d processors for parallel calculations\n", maxNumThreads);
	}
}

void SRmodel::UpdateCustomCriterion(double stress[])
{
	double custom;

	//Replace this with whatever stress criterion is desired
	{
		custom = math.GetSvm(stress);
	}


	if (custom > maxCustom)
		maxCustom = custom;
}

void SRmodel::zeroStressMax()
{
	stressMax = 0.0;
	for (int i = 0; i < 6; i++)
		stressMaxComp[i] = 0.0;
	maxsp1 = 0.0;
	minsp2 = BIG;
	maxCustom = 0.0;
};











