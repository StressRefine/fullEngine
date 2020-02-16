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
// SRpostProcess.cpp: implementation of the SRpostProcess class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <stdlib.h>
#include <search.h>
#include <string>
#include "SRmodel.h"
#include "SRelement.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

extern SRmodel model;

SRpostProcess::SRpostProcess()
{
	noTopoFilter = false;
	minElForBreakout = 1000;
	elIdAtRad0 = -1;
}

void SRpostProcess::Cleanup()
{
	for (int i = 0; i < 6; i++)
		smoothedStrains[i].Free();
	elSmooth.Free();
}

double SRpostProcess::GlobalStrainSmooth()
{
	//global strain smoothing
	//since strain continuity is not assured across material interfaces,
	//smooth strains in groups of elements with same material
	//double sInit = dsecnd();

	int n = model.GetNumFunctions();
	model.AllocateSmoothFunctionEquations(n);
	model.AllocateSkipFun(n);

	int numel = model.GetNumElements();
	elSmooth.Allocate(numel);
	double strainMax = 0.0;
	SRelement* elem;
	if (model.allMatsEquivElast)
	{
		for (int j = 0; j < numel; j++)
		{
			elem = model.GetElement(j);
			elSmooth.Put(j, 1);
		}
		strainMax = ElementStrainSmoothByMaterial();
		for (int j = 0; j < numel; j++)
		{
			elem = model.GetElement(j);
			elem->DownloadSmoothedStrains();
			elem->checkMaterialStrainMax(strainMax);
		}
		return strainMax;
	}

	for (int i = 0; i < model.GetNumMaterials(); i++)
	{
		double strainMaxMat = 0.0;
		elSmooth.Zero();
		bool anyElToSmooth = false;
		for (int j = 0; j < numel; j++)
		{
			elem = model.GetElement(j);
			if (elem->GetMaterialId() == i)
			{
				elSmooth.Put(j, 1);
				anyElToSmooth = true;
			}
		}

		if (anyElToSmooth)
		{
			strainMaxMat = ElementStrainSmoothByMaterial();
			if (strainMaxMat > strainMax)
				strainMax = strainMaxMat;

			for (int j = 0; j < numel; j++)
			{
				elem = model.GetElement(j);
				if (elSmooth.Get(j) == 1)
				{
					elem->DownloadSmoothedStrains();
					elem->checkMaterialStrainMax(strainMaxMat);
				}
			}
		}

	}

	//double sElapsed = (dsecnd() - sInit);
	//OUTPRINT("global smooth elapsed secs: %lg", sElapsed);

	return strainMax;
}

double SRpostProcess::ElementStrainSmoothByMaterial()
{
	//interpolate Strain with same basis functions as were used in last solution
	//for displacement. Determine coefficients by least-squares fitting to "raw
	//strains"
	//return:
		//maximum "raw" strain in model
	//notes:
		//this is called once for each material. The elements that have the currently active material
		//are stored in model.elSmooth.

		//sacrificial elements are smoothed together with other elements of same material to avoid discontinuity at 
		//boundary. sacrificial elements are frozen at p2 to minimize effects of singularities

	int e, i, j, nint, gp, gi, gj, stiffLoc, comp, symloc, elmax;
	double r, s, t, w, bi, bj, detJ, strain[6], ae, strainMax = 0.0, stress[6], svm, svmmax = 0.0;
	double eT, eTmax = 0.0;
	SRelement* elem;

	model.NumberEquationsSmooth();
	
	int n = model.GetNumSmoothEquations();
	if(n == 0)
		return 0.0;

	SRsolver* sol = &(model.solver);
	SRpardiso* par = NULL;
	par = sol->parDisoPtr();
	par->smoothBookkeep();
	double* strainRhsvec[6];

	for(i = 0; i < 6; i++)
	{
		smoothedStrains[i].Allocate(n);
		strainRhsvec[i] = smoothedStrains[i].GetVector();
	}
	int nfunel;
	int len, eqi, eqj;

	for(e = 0; e < model.GetNumElements(); e++)
	{
		double elSvmMax = 0.0;
		double* stiff = model.GetElementStiffnessVector();
		elem = model.GetElement(e);
		elem->PutThread(0);
		elem->SetBasisData();
		elem->SetSvmMax(0.0);
		nfunel = elem->GetNumFunctions();
		elem->FillStiffDiag(nfunel);
		nint = model.math.FillGaussPoints(elem);
		if(elSmooth.Get(e) == 0)
			continue;
		elem->AllocateRawStrains(nint);
		len = nfunel*(nfunel + 1)/2;
		for(i = 0; i < len; i++)
			stiff[i] = 0.0;
		double etx, ety, etz;
		for (gp = 0; gp < nint; gp++)
		{
			model.math.GetGP3d(gp, r, s, t, w);
			detJ = elem->FillMapping(r, s, t);
			double* basisvec = elem->FillBasisFuncs(r, s, t, both);

			eT = elem->CalculateRawStrain(r, s, t, strain, etx, ety, etz);

			if(eT > eTmax)
				eTmax = eT;
			elem->StraintoStress(r, s, t, strain, stress);
			svm = model.math.GetSvm(stress);
			if (svm > elSvmMax)
				elSvmMax = svm;
			if(svm > elem->GetSvmMax())
				elem->SetSvmMax(svm);
			if (!elem->isSacrificial())
			{
				for (comp = 0; comp < 6; comp++)
				{
					ae = fabs(strain[comp]);
					if (ae > strainMax)
						strainMax = ae;
				}
			}

			elem->PutRawStrain(gp, strain);
			w *= detJ;
			for(i = 0; i < nfunel; i++)
			{
				gi = elem->GetFunctionNumber(i);
				eqi = model.GetSmoothFunctionEquation(gi);
				SRASSERT(eqi >= 0);
				bi = basisvec[i] * w;
				for (comp = 0; comp < 6; comp++)
					strainRhsvec[comp][eqi] += (strain[comp]*bi);
				for (j = i; j < nfunel; j++)
				{
					bj = basisvec[j];
					symloc = elem->GetStiffnessLocationAboveDiag(i, j);
					stiff[symloc] += bi*bj;
				}
			}
		}

		if (elSvmMax > svmmax)
		{
			svmmax = elSvmMax;
			elmax = e;
		}

		par->smoothAssemble(e, stiff);
	}

	bool smoothing = true;
	par->solve(smoothing);
	par->clear();

	return strainMax;
}

void SRpostProcess::PostProcess()
{
	//postprocessing:
	//output mesh, displacement results, and stress results

	if (model.breakout)
	{
		PostProcessForBreakout();
		return;
	}

	model.zeroStressMax();
	for (int i = 0; i < model.GetNumMaterials(); i++)
	{
		SRmaterial* mat = model.GetMaterial(i);
		mat->PutMaxSvm(0.0);
	}

	int i, j, dof;
	SRnode* node;

	int numnode = model.GetNumNodes();

	//default output: Cgx (Calculix) .frd format

	if (!model.cgxFrdFile.OutOpenNoFail())
	{
		OUTPRINT("unable to open cgx frd file for postprocessing");
		OUTPRINT("this model may be open already in Cgx postprocessor");
		REPPRINT("unable to open cgx frd file for postprocessing");
		REPPRINT("this model may be open already in Cgx postprocessor");
		ERROREXIT;
	}


	//mesh definition:
	int numNodeOut = MeshToFrd();

	//stresses:
	PostProcessElementStresses();

	//stress at each node is now stored in nodalStress (node num, component num)
	FRDPRINT("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	FRDPRINT(" -4  STRESS      6    1");
	FRDPRINT(" -5  SXX         1    4    1    1");
	FRDPRINT(" -5  SYY         1    4    2    2");
	FRDPRINT(" -5  SZZ         1    4    3    3");
	FRDPRINT(" -5  SXY         1    4    1    2");
	FRDPRINT(" -5  SYZ         1    4    2    3");
	FRDPRINT(" -5  SZX         1    4    3    1");
	double stressComp;
	double stressConv = model.stressUnitConversion;
	for (i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		if (!node->checkUnsupOrBktCon())
		{
			FRDPRINTNORET(" -1%10d", node->GetUserid());
			FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, xxComponent));
			FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, yyComponent));
			FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, zzComponent));
			FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, xyComponent));
			FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, yzComponent));
			FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, xzComponent));
			FRDPRINTRET;
		}
		else
		{
			//node touches boundary of breakout model or is  unsupported:
			FRDPRINT(" -1%10d 0.0 0.0 0.0 0.0 0.0 0.0", node->GetUserid());
		}
	}
	FRDPRINT(" -3");

	//displacements
	FRDPRINT("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	FRDPRINT(" -4  DISP        4    1");
	FRDPRINT(" -5  D1          1    2    1    0");
	FRDPRINT(" -5  D2          1    2    2    0");
	FRDPRINT(" -5  D3          1    2    2    0");
	FRDPRINT(" -5  ALL         1    2    0    0    1ALL");
	SRvec3 disp;
	for (i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		node->GetDisp(disp);
		disp.Scale(model.lengthUnitConversion);
		FRDPRINT(" -1%10d%12.5le%12.5le%12.5le", node->GetUserid(), disp.d[0], disp.d[1], disp.d[2]);
	}
	FRDPRINT(" -3");

	//PFringeOut();

	FRDPRINT("9999"); //end of data

	model.cgxFrdFile.Close();

	if (model.outputf06)
		OutputF06();

	nodalStress.Free();
	nodeDisps.Free();
}

void SRpostProcess::PostProcessForBreakoutPartial()
{
	//postprocessing:
	//output mesh, displacement results, and stress results to "partial frd" file
	//(breakout model only)

	model.zeroStressMax();
	for (int i = 0; i < model.GetNumMaterials(); i++)
	{
		SRmaterial* mat = model.GetMaterial(i);
		mat->PutMaxSvm(0.0);
	}

	int i, j, dof;
	SRstring line, tok;

	//default output: Cgx (Calculix) .frd format

	if (!model.cgxFrdFile.OutOpenNoFail())
	{
		OUTPRINT("unable to open cgx frd file for postprocessing");
		OUTPRINT("this model may be open already in Cgx postprocessor");
		REPPRINT("unable to open cgx frd file for postprocessing");
		REPPRINT("this model may be open already in Cgx postprocessor");
		ERROREXIT;
	}

	int nuid;
	SRnode* node;
	int numNodes = model.GetNumNodes();

	//mesh definition:
	int numNodeOut = MeshToFrd();

	//stresses:
	PostProcessElementStresses();

	clipBreakoutSacrStreses();

	//stress at each node for breakout model is now stored in nodalStress (node num, component num)
	FRDPRINT("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	FRDPRINT(" -4  STRESS      6    1");
	FRDPRINT(" -5  SXX         1    4    1    1");
	FRDPRINT(" -5  SYY         1    4    2    2");
	FRDPRINT(" -5  SZZ         1    4    3    3");
	FRDPRINT(" -5  SXY         1    4    1    2");
	FRDPRINT(" -5  SYZ         1    4    2    3");
	FRDPRINT(" -5  SZX         1    4    3    1");
	double stressComp;
	double conv = model.stressUnitConversion;
	for (i = 0; i < numNodes; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		nuid = node->userId;
		if (!node->checkUnsupOrBktCon())
		{
			FRDPRINTNORET(" -1%10d", nuid);
			FRDPRINTNORET("%12.5le", conv*nodalStress.Get(i, xxComponent));
			FRDPRINTNORET("%12.5le", conv*nodalStress.Get(i, yyComponent));
			FRDPRINTNORET("%12.5le", conv*nodalStress.Get(i, zzComponent));
			FRDPRINTNORET("%12.5le", conv*nodalStress.Get(i, xyComponent));
			FRDPRINTNORET("%12.5le", conv*nodalStress.Get(i, yzComponent));
			FRDPRINTNORET("%12.5le", conv*nodalStress.Get(i, xzComponent));
			FRDPRINTRET;
		}
		else
		{
			//node touches boundary of breakout model or is  unsupported:
			FRDPRINT(" -1%10d 0.0 0.0 0.0 0.0 0.0 0.0", nuid);
		}
	}
	FRDPRINT(" -3");

	//displacements
	FRDPRINT("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	FRDPRINT(" -4  DISP        4    1");
	FRDPRINT(" -5  D1          1    2    1    0");
	FRDPRINT(" -5  D2          1    2    2    0");
	FRDPRINT(" -5  D3          1    2    2    0");
	FRDPRINT(" -5  ALL         1    2    0    0    1ALL");
	conv = model.lengthUnitConversion;
	SRvec3 disp;
	for (i = 0; i < numNodes; i++)
	{
		SRnode* node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		nuid = node->userId;
		node->GetDisp(disp);
		disp.Scale(conv);
		FRDPRINT(" -1%10d%12.5le%12.5le%12.5le", nuid, disp.d[0], disp.d[1], disp.d[2]);
	}
	FRDPRINT(" -3");

	FRDPRINT("9999"); //end of data

	model.cgxFrdFile.Close();

	if (model.outputf06)
		OutputF06();

	nodalStress.Free();
	nodeDisps.Free();
}

void SRpostProcess::PostProcessForBreakout()
{
	//postprocessing:
	//output mesh, displacement results, and stress results to "full frd" file

	//mark nodes of elements on bktelsacrlist as "ownedByBktEl"
	for (int e = 0; e < bktElSacrList.GetNum(); e++)
	{
		int euid = bktElSacrList.Get(e);
		int eid = model.input.elemFind(euid);
		if (eid == -1)
			continue;
		SRelement* elem = model.GetElement(eid);
		for (int n = 0; n < elem->GetNumNodesTotal(); n++)
		{
			SRnode* node = elem->getNodeOrMidNode(n);
			node->sacrificialType = ownedByBktEl;
		}
	}
	model.zeroStressMax();
	PostProcessForBreakoutPartial();
}

void SRpostProcess::PostProcessElementStresses()
{
	//output stress results for each element
	for (int i = 0; i < model.GetNumMaterials(); i++)
	{
		SRmaterial* mat = model.GetMaterial(i);
		mat->SetVolPercentYielded(0.0);
	}

	int e;
	SRelement* elem;

	int nn = model.GetNumNodes();
	nodalStress.Allocate(nn, 6);
	nodalStressCount.Allocate(nn);

	for (e = 0; e < model.GetNumElements(); e++)
	{
		elem = model.GetElement(e);
		elem->PutThread(0);
		CalculateElementStresses(elem);
	}

	int nel = model.GetNumElements();
	bool doMaxClipping = (nel - model.errorChecker.GetnumSacrificialElements() > 0);
	for (e = 0; e < nel; e++)
	{
		elem = model.GetElement(e);
		elem->PutThread(0);
		if (elem->isSacrificial())
			fillSacricialElementNodalStress(elem, doMaxClipping);
	}

	nodalStressCount.Free();
}

void SRpostProcess::PFringeOut()
{
	SRnode *node;
	int numnode = model.GetNumNodes();
	nodalPOrders.Allocate(numnode);
	nodalPOrders.Set(2);
	for (int e = 0; e < model.GetNumEdges(); e++)
	{
		SRedge* edge = model.GetEdge(e);
		int p = edge->pOrder;
		int nid = edge->midnodeId;
		if (p > nodalPOrders.d[nid])
			nodalPOrders.d[nid] = p;
		for (int n = 0; n < 2; n++)
		{
			nid = edge->nodeIds[n];
			if (p > nodalPOrders.d[nid])
				nodalPOrders.d[nid] = p;
		}
	}

	int numNodeOut = 0;
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		numNodeOut++;
	}

	FRDPRINT("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	FRDPRINT(" -4  PORDER  1    1");
	FRDPRINT(" -5  P1      1    1    1    0");
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		FRDPRINTNORET(" -1%10d", node->userId);
		double xp = (double)nodalPOrders.d[i];
		FRDPRINTNORET("%12.5le", xp);
		FRDPRINTRET;
	}
	FRDPRINT(" -3");
}

void SRpostProcess::CalculateElementStresses(SRelement* elem, bool checkMaxOnly)
{
	//calculate stress results for an element

	//input:
		//elem = pointer to element
		//checkMaxOnly = true for check max only, false to store stress in nodalStress 
	//note:
		//nodalStress stores stress at each node in model
		//if global stress smoothing has not been performed, stress contribution
		//from each element can be different so averaging is performed

	static int numNodesTotalBrick = 27; //nodes, centroid, face centroids
	static int numNodesTotalTet = 15; //nodes, centroid, face centroids
	static int numNodesTotalWedge = 21; //nodes, centroid, face centroids
	int i, nn, ne, nnt;
	double r, s, t;
	double stress[6];

nn = elem->GetNumNodes();
	ne = elem->GetNumLocalEdges();
	nnt = nn + ne;
	SRnode* node;
	SRedge* edge;
	int id;
	double elSvmMax = 0.0;
	elem->SetSvmMax(0.0);
	double svm;
	int threadnum = elem->GetThread();
	elem->SetBasisData(threadnum);
	bool hasAllowable = !checkMaxOnly && model.allMatsHaveAllowable && elem->GetMaterial()->isAllowableAssigned();
	double Allowable = 0.0;
	int matid = elem->GetMaterial()->getId();
	if(hasAllowable)
		Allowable = elem->GetMaterial()->GetAllowableStress();
	int numNodesYielded = 0;
	bool yielded = false;
	for (i = 0; i < nn; i++)
	{
		id = elem->GetNodeId(i);
		node = model.GetNode(id);
		elem->NodeNaturalCoords(i, r, s, t);
		elem->GetStress(r, s, t, stress, !checkMaxOnly);
		if (!elem->isSacrificial())
		{
			svm = StressMaxCheck(elem, node->Position(), stress);
			if (svm > elSvmMax)
				elSvmMax = svm;
			if (hasAllowable && !checkMaxOnly)
			{
				if (svm > Allowable)
					numNodesYielded++;
			}
			if (!checkMaxOnly)
			{
				if (nodalStressCount.d[id] == 0)
				{
					for (int c = 0; c < 6; c++)
						nodalStress.Put(id, c, stress[c]);
					nodalStressCount.d[id]++;
				}
			}
		}
	}
	for(i = 0; i < ne; i++)
	{
		edge = elem->GetEdge(i);
		id = edge->GetMidNodeId();
		node = model.GetNode(id);
		elem->NodeNaturalCoords(i + nn, r, s, t);
		elem->GetStress(r, s, t, stress, !checkMaxOnly);
		if (!elem->isSacrificial())
		{
			svm = StressMaxCheck(elem, node->Position(), stress);
			if (svm > elSvmMax)
				elSvmMax = svm;
			if (hasAllowable && !checkMaxOnly)
			{
				if (svm > Allowable)
					numNodesYielded++;
			}
			if (!checkMaxOnly)
			{
				if (nodalStressCount.d[id] == 0)
				{
					for (int c = 0; c < 6; c++)
						nodalStress.Put(id, c, stress[c]);
					nodalStressCount.d[id]++;
				}
			}
		}
	}

	//also sample face centroids and element centroids:
	SRvec3 pos;
	double rf, sf;
	for (i = 0; i < elem->GetNumLocalFaces(); i++)
	{
		SRface* face = elem->GetFace(i);
		model.map.FaceCentroid(face, rf, sf);
		model.map.ElementNaturalCoordsFromFace(elem, i, rf, sf, r, s, t);
		elem->GetStress(r, s, t, stress, !checkMaxOnly);
		if (!elem->isSacrificial())
		{
			face->Position(rf, sf, pos);
			svm = StressMaxCheck(elem, pos, stress);
			if (hasAllowable && !checkMaxOnly)
			{
				if (svm > Allowable)
					numNodesYielded++;
			}

			if (svm > elSvmMax)
				elSvmMax = svm;
		}
	}
	model.map.ElementCentroid(elem, r, s, t);
	elem->GetStress(r, s, t, stress, !checkMaxOnly);
	if (!elem->isSacrificial())
	{
		elem->Position(r, s, t, pos);
		svm = StressMaxCheck(elem, pos, stress);
		if (hasAllowable && !checkMaxOnly)
		{
			if (svm > Allowable)
				numNodesYielded++;
		}
		if (svm > elSvmMax)
			elSvmMax = svm;
	}

	elem->SetSvmMax(elSvmMax);
	StressMaxCheckvsAllowable(elem, elSvmMax);

	if (!checkMaxOnly)
	{
		int nnTotal = numNodesTotalTet;
		if (elem->GetType() == brick)
			nnTotal = numNodesTotalBrick;
		else if (elem->GetType() == wedge)
			nnTotal = numNodesTotalWedge;
		double percentYielded = ((double)numNodesYielded) / ((double)nnTotal);
		elem->GetMaterial()->AddToVolPercentYielded(percentYielded);
	}
}

void SRpostProcess::fillSacricialElementNodalStress(SRelement* elem, bool doMaxClipping)
{
	//for sacrifical elements only, fill up stress in nodes not owned by nonsacrifical elements.
	//but clip back to max in material the element belongs to

	//input:
	//elem = pointer to element
	//note:
	//fills contribution of this element to nodalStress
	int i, nn, ne, nnt;
	double r, s, t;
	double stress[6];
	nn = elem->GetNumNodes();
	ne = elem->GetNumLocalEdges();
	nnt = nn + ne;
	SRedge* edge;
	int id;
	int threadnum = elem->GetThread();
	elem->SetBasisData(threadnum);
	SRmaterial* elmat = elem->GetMaterial();
	double maxSvm = elmat->GetMaxSvm();
	double targetRatio = 0.9;
	//sanity check maxSvm should be <= max in model
	if (maxSvm > model.stressMax)
		maxSvm = model.stressMax;
	//keep it down to targetRatio so it doesn't show up as hot spot in fringes:
	maxSvm *= targetRatio;
	double ratc[6];
	for (i = 0; i < nn; i++)
	{
		id = elem->GetNodeId(i);
		SRnode* node = model.GetNode(id);
		if (doMaxClipping && nodalStressCount.Get(id) != 0)
			continue;
		elem->NodeNaturalCoords(i, r, s, t);
		elem->GetStress(r, s, t, stress);
		double vm = model.math.GetSvm(stress);
		double ratvm = 1.0;
		for (int j = 0; j < 6; j++)
			ratc[j] = 1.0;
		if (doMaxClipping)
		{
			if (vm > maxSvm && vm > TINY)
				ratvm = maxSvm / vm;
			for (int c = 0; c < 6; c++)
			{
				double maxstress = targetRatio*fabs(model.GetStressMaxComp(c));
				double fstress = fabs(stress[c]);
				if (fstress > TINY)
				{
					ratc[c] = maxstress / fstress;
					if (ratvm < ratc[c])
						ratc[c] = ratvm;
				}
			}
		}

		nodalStressCount.Put(id, 1);
	}
	for (i = 0; i < ne; i++)
	{
		edge = elem->GetEdge(i);
		if (edge->GetPorder() < 2)
			continue;
		id = edge->GetMidNodeId();
		SRnode* node = model.GetNode(id);
		if (doMaxClipping && nodalStressCount.Get(id) != 0)
			continue;

		elem->NodeNaturalCoords(i + nn, r, s, t);
		elem->GetStress(r, s, t, stress);
		double vm = model.math.GetSvm(stress);
		double ratvm = 1.0;
		for (int j = 0; j < 6; j++)
			ratc[j] = 1.0;
		if (doMaxClipping)
		{
			if (vm > maxSvm && vm > TINY)
				ratvm = maxSvm / vm;
			for (int c = 0; c < 6; c++)
			{
				double maxstress = targetRatio*fabs(model.GetStressMaxComp(c));
				double fstress = fabs(stress[c]);
				if (fstress > TINY)
				{
					ratc[c] = maxstress / fstress;
					if (ratvm < ratc[c])
						ratc[c] = ratvm;
				}
			}
			for (int c = 0; c < 6; c++)
			{
				stress[c] *= ratc[c];
				nodalStress.Put(id, c, stress[c]);
			}

			for (int c = 0; c < 6; c++)
			{
				stress[c] *= ratc[c];
				nodalStress.Put(id, c, stress[c]);
			}

			nodalStressCount.Put(id, 1);
		}
		nodalStressCount.Put(id, 1);
	}
}

double SRpostProcess::StressMaxCheck(SRelement* elem, SRvec3& pos, double stress[])
{
	//check if stress exceeds max in model
	//input:
		//pos = xyz position
		//stress = stress tensor stored as vector

	double svm = model.math.GetSvm(stress);
	if (svm> model.GetStressMax())
		model.SetStressMax(elem, pos, svm);
	model.SetStressMaxComp(stress);
	model.UpdateCustomCriterion(stress);
	return svm;
}

void SRpostProcess::StressMaxCheckvsAllowable(SRelement* elem, double svmMax)
{
	//check the max stress in an element vs the allowable stress for the material assigned to the element
	//print warning if it is exceeded
	SRmaterial* mat = elem->GetMaterial();
	if (!elem->isSacrificial())
	{
		if (svmMax > mat->GetMaxSvm())
			mat->PutMaxSvm(svmMax);
	}
	if (mat->GetHighStressWarned())
		return;
	double allowstress = mat->GetAllowableStress();
	if (allowstress < TINY)
		return;
	if (svmMax > allowstress)
	{
		OUTPRINT(" warning. allowable stress exceeded for material: %s", mat->GetName());
		mat->SetHighStressWarned(true);
		model.anyMatYielded = true;
	}
}

void SRpostProcess::CalculateMaxStress()
{
	//check for max stress during p loop
	//note:
	//this uses checkMaxOnly with output set to true
	//so stresses are sampled but not output.
	//CalculateElementStresses samples each element node (corner and midside)
	//check for max stress during p loop
	int i;
	SRelement* elem;
	model.zeroStressMax();

	bool checkmaxonly = true;
	if (model.breakout && model.checkForHotSpotatMax)
	{
		checkmaxonly = false;
		int nn = model.GetNumNodes();
		nodalStress.Allocate(nn, 6);
		nodalStressCount.Allocate(nn);
	}

	for (i = 0; i < model.GetNumElements(); i++)
	{
		elem = model.GetElement(i);
		elem->PutThread(0);
		if (!elem->isSacrificial())
			CalculateElementStresses(elem, checkmaxonly);
	}

	if (model.breakout && model.checkForHotSpotatMax)
		nodalStressCount.Free();
}


void SRpostProcess::printNodalForce(int nodeUid, SRforce* force, SRfile&f)
{
	f.Print(" %d", nodeUid);
	if (force->isGcs())
		f.Print(" gcs");
	else if (force->isPressure())
		f.Print(" pressure %lg", force->GetForceVal(0, 0));
	else
	{
		int cid = force->GetCoordId();
		f.Print(" coord %s", model.GetCoord(cid)->GetName());
	}
	if (!force->isPressure())
	{
		for (int dof = 0; dof < 3; dof++)
			f.Print(" %lg", force->GetForceVal(0, dof));
	}
	f.PrintReturn();
}

double SRpostProcess::GetMaxFaceRotation()
{
	double maxRot = 0.0;
	for (int i = 0; i < model.GetNumFaces(); i++)
	{
		double rot = model.GetFace(i)->normalRotation();
		if (rot > maxRot)
			maxRot = rot;
	}
	return maxRot;
}

int SRpostProcess::MeshToFrd()
{
	//mesh definition, nodes and elements:
	//local edge number mapping between SR and cgx:
	int brickSRtoCgx[12] = { 0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7 };
	int wedgeSRtoCgx[9] = { 0, 1, 2, 6, 7, 8, 3, 4, 5 };
	int tetSRtoCgx[6] = { 0, 1, 2, 3, 4, 5 };

	FRDPRINT("    1C"); // ''1C'' defines a new calc
	FRDPRINT("    1UDATE   26.march.2000"); // ''1U'' stores user job - information, can be any string, ttd put real date, add more lines, e.g. PGM stressrefine, model file path 
	//nodes:
	int i, j, dof;
	SRnode* node;
	int numnode = model.GetNumNodes();
	int numNodeOut = 0;
	for (i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		numNodeOut++;
	}
	FRDPRINT("    2C                  %12d                                     1", numNodeOut);
	for (i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		//http://bconverged.com/calculix/doc/cgx/html/node167.html
		//long Format:(1X,'-1',I10,3E12.5)
		FRDPRINT(" -1%10d%12.5le%12.5le%12.5le", node->GetUserid(), node->GetXyz(0), node->GetXyz(1), node->GetXyz(2));
	}
	FRDPRINT(" -3");
	//elements. eltype: 10 node tet 6, wedge 5, brick 4
	int numel = model.GetNumElements();
	FRDPRINT("    3C                  %12d                                     1", numel);
	int *srToCgx;
	for (i = 0; i < numel; i++)
	{
		SRelement* elem = model.GetElement(i);
		int type;
		bool needContinue = false;
		if (elem->GetType() == brick)
		{
			needContinue = true;
			type = 4;
			srToCgx = brickSRtoCgx;
		}
		else if (elem->GetType() == wedge)
		{
			needContinue = true;
			type = 5;

			srToCgx = wedgeSRtoCgx;
		}
		else
		{
			type = 6;
			srToCgx = tetSRtoCgx;
		}
		FRDPRINT(" -1%10d%5d%5d%5d", elem->GetUserid(), type, 0, 1);
		FRDPRINTNORET(" -2");
		int nn = elem->GetNumNodes();
		for (j = 0; j < nn; j++)
			FRDPRINTNORET("%10d", elem->GetNode(j)->GetUserid());
		if (!model.saveBreakout)
		{
			for (j = 0; j < elem->GetNumLocalEdges(); j++)
			{
				int edgeNum = srToCgx[j];
				FRDPRINTNORET("%10d", elem->GetEdge(edgeNum)->GetMidNodeUserId());
				nn++;
				if (nn == 10 && needContinue)
				{
					FRDPRINTRET;
					FRDPRINTNORET(" -2");
				}
			}
		}
		FRDPRINTRET;
	}
	FRDPRINT(" -3");
	return numNodeOut;
}

void SRpostProcess::writeFrd()
{
	if (!model.cgxFrdFile.OutOpenNoFail())
	{
		return;
	}

	//mesh definition:
	int numNodeOut = MeshToFrd();

	//displacements
	FRDPRINT("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	FRDPRINT(" -4  DISP        4    1");
	FRDPRINT(" -5  D1          1    2    1    0");
	FRDPRINT(" -5  D2          1    2    2    0");
	FRDPRINT(" -5  D3          1    2    2    0");
	FRDPRINT(" -5  ALL         1    2    0    0    1ALL");
	SRvec3 disp;
	for (int i = 0; i < model.GetNumNodes(); i++)
	{
		SRnode* node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		disp = nodeDisps.Get(i);
		disp.Scale(model.lengthUnitConversion);
		FRDPRINT(" -1%10d%12.5le%12.5le%12.5le", node->GetUserid(), disp.d[0], disp.d[1], disp.d[2]);
	}
	FRDPRINT(" -3");

	FRDPRINT("9999"); //end of data

	model.cgxFrdFile.Close();

}

void SRpostProcess::readPreviousResultsF06()
{
	SCREENPRINT("reading displacment and stress (f06) results.\nsource:\n");
	SCREENPRINT(" %s\n", model.f06File.filename.str);
	OUTPRINT(" displacment and stress (f06) results source:");
	OUTPRINT(" %s", model.f06File.filename.str);
	SRfile& f = model.f06File;
	bool dispok = true;
	if (!f.Open(SRinputMode))
		dispok = false;
	else
		dispok = readPreviousDisplacementF06();
	if (!dispok)
	{
		OUTPRINT(" displacment and stress (f06) results are needed for breakout models and");
		OUTPRINT(" unsupported entities.");
		REPPRINT("Nastran displacement results not available or incomplete");
		ERROREXIT;
	}

	if (model.saveBreakout && model.saveBreakoutData.atMax)
		readPreviousStressF06();
	f.Close();
}

bool SRpostProcess::getF06DispOrStressLineCheckForSkip(SRfile& f, SRstring& line)
{
	//get next line to parse for disp or stress.
	//catch occurrence of new page header line- these
	//are only lines in file that have "page" in them
	if (!f.GetLine(line))
		return false;
	while (1)
	{
		if (!line.isAllBlank())
			break;
		if (!f.GetLine(line))
			return false;
	}
	SRstring linesav;
	linesav = line;
	string sline;
	sline.assign(line.str);
	int pageLoc = sline.find("  PAGE");
	if (pageLoc >= 0)
	{
		//f06 file has new page so skip header:
		for (int i = 0; i < 6; i++)
		{
			if (!f.GetLine(line))
				return false;
		}
		if (!f.GetLine(line))
			return false;
	}
	else
		line = linesav;
	return true;
}

bool SRpostProcess::readPreviousDisplacementF06()
{
	//read displacement results for the model from a previous run from a Nastran f06 file
	//notes:
	//if minPorder =1, reads displacements at corners, converts to p1 solution
	//else reads displacements at corners and midnodes, converts to p2 solution
	bool debugF06 = false;

	SRfile& f = model.f06File;
	int dof, fun;
	SRnode* node;
	SRvec3 disp;
	SRstring line, tok;

	int nnode = model.GetNumNodes();

	bool header1Found = false;
	while (1)
	{
		//read header lines until find "D I S P L A C E M E N T   V E C T O R" header:
		if (!f.GetLine(line))
		{
			return false;
		}
		if (!header1Found)
		{
			char* str = line.FirstChar('D');
			if (str != NULL)
			{
				tok.Copy(str);
				if (tok.CompareUseLength("D I S P L A C E M E N T   V E C T O R"))
					header1Found = true;
			}
		}
		else
		{
			tok = line.Token();
			if (tok.Compare("POINT"))
				break;
		}

	}

	if (!header1Found)
		return false;

	//displacements are in lines following the header until a blank is found
	double umagmax = 0.0;

	int numnonorphan = 0;
	for (int n = 0; n < nnode; n++)
	{
		SRnode* node = model.GetNode(n);
		if (!node->isOrphan())
			numnonorphan++;
	}

	int numdispAssigned = 0;
	while (1)
	{
		if (!getF06DispOrStressLineCheckForSkip(f, line))
			break;
		SRstring lineIn;
		lineIn = line;
		int nuid;
		line.TokRead(nuid);
		tok = line.Token(); //type field. make sure it is a grid:
		if (!tok.Compare("G"))
			continue;
		SRnode* node = model.GetNodeFromUid(nuid);
		if (node == NULL)
			continue;

		if (node->isOrphan())
			continue;

		node->hasDisp = true;
		line.TokRead(disp.d[0]);
		line.TokRead(disp.d[1]);
		line.TokRead(disp.d[2]);
		nodeDisps.Put(node->id, disp);
		double umag = disp.Magnitude();
		if (umag > umagmax)
		{
			umagmax = umag;
			model.umaxInModel.Copy(disp);
			model.nodeuidAtumax = nuid;
		}
		numdispAssigned++;
		if (numdispAssigned == numnonorphan)
			break;
	}

	//sanity check: all non-orphan nodes have to have disps:
	for (int n = 0; n < nnode; n++)
	{
		SRnode* node = model.GetNode(n);
		if (!node->isOrphan() && !node->hasDisp)
		{
			REPPRINT("Error: Nastran Displacement results file (.f06) incompatible with model file");
			ERROREXIT;
		}
	}
	return true;
}


void SRpostProcess::readPreviousStressF06()
{
	//make sure stress records exist in F06 file
	SRfile& f = model.f06File;
	int nnode = model.nodes.GetNum();
	SRstring line, tok;
	line.setTokSep(" ");
	int numel = model.GetNumElements();
	//read header lines until find "S T R E S S" header.
	// Note have to also look for "C H E X A" or "C P E N T A" or "C T E T R A" or  for mixed meshes
	//so shells, etc are skipped
	bool stressHeaderFound = false;
	while (1)
	{
		if (!f.GetLine(line))
			break;
		char* str = line.FirstChar('S');
		if (str != NULL)
		{
			tok.Copy(str);
			if (tok.CompareUseLength("S T R E S S E S   I N"))
			{
				str = tok.FirstChar('(');
				line.Copy(str);
				if (line.CompareUseLength("( C H E X A") || line.CompareUseLength("( C P E N T A") || line.CompareUseLength("( C T E T R A"))
				{
					stressHeaderFound = true;
					break;
				}
			}
		}
	}
	int skip;
	bool stressesOK = true;
	if (!stressHeaderFound)
		stressesOK = false;
	else
	{
		//two more header lines
		for (skip = 0; skip < 2; skip++)
		{
			if (!f.GetLine(line))
			{
				stressesOK = false;
				break;
			}
		}
	}

	if (!stressesOK)
	{
		//fatal error if stresses are needed:
		if (model.saveBreakout && model.saveBreakoutData.atMax)
		{
			REPPRINT("Error: Nastran stress results not available. Cannot solve local region at maximum stress in model");
			REPPRINT("or maximum stress in material");
			ERROREXIT;
		}
		else
			return;
	}


	double svm;

	int numelRead = 0;
	bool stressReadError = false;
	while (1)
	{
		//next line contains eluid or blank line if there is a new page:
		if (!getF06DispOrStressLineCheckForSkip(f, line))
			break;
		SRstring lineSav;
		lineSav = line;
		tok = line.Token();
		int elid, eluid;
		line.TokRead(eluid);
		elid = model.input.elemFind(eluid);
		SRelement* elem = model.GetElement(elid);
		int ncorner = elem->GetNumCorners();
		if (!getF06DispOrStressLineCheckForSkip(f, line))
			break;
		//see if this line is centroid:
		line.Token();
		tok = line.Token();
		bool readCenter = false;
		double stress[6];
		double x, y, z;
		if (tok.CompareUseLength("CENTER"))
		{
			for (skip = 0; skip < 2; skip++)
			{
				if (!getF06DispOrStressLineCheckForSkip(f, line))
					break;
			}
			readCenter = true;
		}
		for (int corner = 0; corner < ncorner; corner++)
		{
			if (readCenter || corner != 0)
			{
				if (!getF06DispOrStressLineCheckForSkip(f, line))
					break;
			} //else don't read line because line just read was 1st corner
			lineSav = line;
			for (skip = 0; skip < 3; skip++)
				tok = line.Token();
			if (!tok.CompareUseLength("X"))
			{
				stressReadError = true;
				break;
			}
			line.TokRead(stress[xxComponent]);
			tok = line.Token();
			if (!tok.CompareUseLength("XY"))
			{
				stressReadError = true;
				break;
			}
			line.TokRead(stress[xyComponent]);
			if (!getF06DispOrStressLineCheckForSkip(f, line))
				break;
			lineSav = line;
			tok = line.Token();
			if (!tok.CompareUseLength("Y"))
			{
				stressReadError = true;
				break;
			}
			line.TokRead(stress[yyComponent]);
			tok = line.Token();
			if (!tok.CompareUseLength("YZ"))
			{
				stressReadError = true;
				break;
			}
			line.TokRead(stress[yzComponent]);
			if (!getF06DispOrStressLineCheckForSkip(f, line))
				break;
			lineSav = line;
			tok = line.Token();
			if (!tok.CompareUseLength("Z"))
			{
				stressReadError = true;
				break;
			}
			line.TokRead(stress[zzComponent]);
			tok = line.Token();
			if (!tok.CompareUseLength("ZX"))
			{
				stressReadError = true;
				break;
			}
			line.TokRead(stress[xzComponent]);
			svm = model.math.GetSvm(stress);
			SRvec3 elcen;
			elem->approxCentroid(elcen);
			if (svm > model.feasvmmax)
			{
				model.feasvmmax = svm;
				model.feavmmaxpos.Copy(elcen);
			}
		}
		if (stressReadError)
		{
			OUTPRINT("Stress Read Error. element %d", eluid);
			OUTPRINT("line from f06 file:");
			OUTPRINT("%s", lineSav.str);
			break;
		}
		numelRead++;
		if (numelRead == numel)
			break;
	}
	if (stressReadError)
	{
		REPPRINT("Error: Nastran stress results not readable. Cannot solve local region at maximum stress in model");
		REPPRINT("or maximum stress in material");
		ERROREXIT;
	}
	if (model.saveBreakout && model.saveBreakoutData.atMax)
		model.saveBreakoutData.origin.Copy(model.feavmmaxpos);
}

void SRpostProcess::PlotEdges(int numEdge, int edgeId[])
{
	SRstring name;
	name.Copy(model.outdir);
	name.Cat("\\PlotEdges");
	name.Cat(".frd");
	SRfile frdf;
	frdf.Open(name, SRoutputMode);

	//mesh definition:
	int numNodeOut = 0;
	int brickSRtoCgx[12] = { 0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7 };
	int wedgeSRtoCgx[9] = { 0, 1, 2, 6, 7, 8, 3, 4, 5 };
	int tetSRtoCgx[6] = { 0, 1, 2, 3, 4, 5 };
	frdf.PrintLine("    1C"); // ''1C'' defines a new calc
	frdf.PrintLine("    1UDATE   26.march.2000"); // ''1U'' stores user job - information, can be any string, ttd put real date, add more lines, e.g. PGM stressrefine, model file path 
	//nodes:
	SRnode* node = NULL;
	int numnode = model.GetNumNodes();
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		numNodeOut++;
	}
	frdf.PrintLine("    2C                  %12d                                     1", numNodeOut);
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		//http://bconverged.com/calculix/doc/cgx/html/node167.html
		//long Format:(1X,'-1',I10,3E12.5)
		frdf.PrintLine(" -1%10d%12.5le%12.5le%12.5le", node->userId, node->pos.d[0], node->pos.d[1], node->pos.d[2]);
	}
	frdf.PrintLine(" -3");

	int numel = model.GetNumElements();
	frdf.PrintLine("    3C                  %12d                                     1", numel);
	int *srToCgx;
	for (int e = 0; e < numel; e++)
	{
		SRelement* elem = model.GetElement(e);
		int type;
		bool needContinue = false;
		if (elem->type == brick)
		{
			needContinue = true;
			type = 4;
			srToCgx = brickSRtoCgx;
		}
		else if (elem->type == wedge)
		{
			needContinue = true;
			type = 5;

			srToCgx = wedgeSRtoCgx;
		}
		else
		{
			type = 6;
			srToCgx = tetSRtoCgx;
		}
		frdf.PrintLine(" -1%10d%5d%5d%5d", elem->userId, type, 0, 1);
		frdf.Print(" -2");
		int nn = elem->GetNumNodes();
		for (int j = 0; j < nn; j++)
			frdf.Print("%10d", elem->GetNode(j)->userId);
		if (!model.saveBreakout)
		{
			for (int j = 0; j < elem->GetNumLocalEdges(); j++)
			{
				int edgeNum = srToCgx[j];
				frdf.Print("%10d", elem->GetEdge(edgeNum)->GetMidNodeUserId());
				nn++;
				if (nn == 10 && needContinue)
				{
					frdf.PrintLine("\n");
					frdf.Print(" -2");
				}
			}
		}
		frdf.PrintLine("\n");
	}
	frdf.PrintLine(" -3");

	SRintVector tmpNodeDisps;
	tmpNodeDisps.Allocate(model.GetNumNodes());
	for (int f = 0; f < numEdge; f++)
	{
		int eid = edgeId[f];
		SRedge* edge = model.GetEdge(eid);
		int nn = 3;
		if (model.saveBreakout)
			nn = 2;
		for (int n = 0; n < nn; n++)
		{
			int nid = edge->GetNodeOrMidNodeId(n);
			tmpNodeDisps.Put(nid, 1);
		}
	}
	//displacements

	frdf.PrintLine("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	frdf.PrintLine(" -4  DISP        4    1");
	frdf.PrintLine(" -5  D1          1    2    1    0");
	frdf.PrintLine(" -5  D2          1    2    2    0");
	frdf.PrintLine(" -5  D3          1    2    2    0");
	frdf.PrintLine(" -5  ALL         1    2    0    0    1ALL");

	SRvec3 disp;
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		disp.Zero();
		if (tmpNodeDisps.Get(i) == 1)
			disp.Assign(1.0, 1.0, 1.0);
		frdf.PrintLine(" -1%10d%12.5le%12.5le%12.5le", node->userId, disp.d[0], disp.d[1], disp.d[2]);
	}
	frdf.PrintLine(" -3");

	frdf.PrintLine("9999"); //end of data

	frdf.Close();

}


void SRpostProcess::PlotFaces(int numFace, SRface* facev[])
{
	SRintVector fidv;
	fidv.Allocate(numFace);
	for (int i = 0; i < numFace; i++)
		fidv.d[i] = facev[i]->id;
	PlotFaces(numFace, fidv.d);
}
void SRpostProcess::PlotFaces(int numFace, int fidv[])
{
	SRstring name;
	name.Copy(model.outdir);
	name.Cat("\\PlotFaces");
	name.Cat(".frd");
	SRfile frdf;
	frdf.Open(name, SRoutputMode);

	//mesh definition:
	int numNodeOut = 0;
	int brickSRtoCgx[12] = { 0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7 };
	int wedgeSRtoCgx[9] = { 0, 1, 2, 6, 7, 8, 3, 4, 5 };
	int tetSRtoCgx[6] = { 0, 1, 2, 3, 4, 5 };
	frdf.PrintLine("    1C"); // ''1C'' defines a new calc
	frdf.PrintLine("    1UDATE   26.march.2000"); // ''1U'' stores user job - information, can be any string, ttd put real date, add more lines, e.g. PGM stressrefine, model file path 
	//nodes:
	SRnode* node = NULL;
	int numnode = model.GetNumNodes();
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		numNodeOut++;
	}
	frdf.PrintLine("    2C                  %12d                                     1", numNodeOut);
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		//http://bconverged.com/calculix/doc/cgx/html/node167.html
		//long Format:(1X,'-1',I10,3E12.5)
		frdf.PrintLine(" -1%10d%12.5le%12.5le%12.5le", node->userId, node->pos.d[0], node->pos.d[1], node->pos.d[2]);
	}
	frdf.PrintLine(" -3");

	int numel = model.GetNumElements();
	frdf.PrintLine("    3C                  %12d                                     1", numel);
	int *srToCgx;
	for (int e = 0; e < numel; e++)
	{
		SRelement* elem = model.GetElement(e);
		int type;
		bool needContinue = false;
		if (elem->type == brick)
		{
			needContinue = true;
			type = 4;
			srToCgx = brickSRtoCgx;
		}
		else if (elem->type == wedge)
		{
			needContinue = true;
			type = 5;

			srToCgx = wedgeSRtoCgx;
		}
		else
		{
			type = 6;
			srToCgx = tetSRtoCgx;
		}
		frdf.PrintLine(" -1%10d%5d%5d%5d", elem->userId, type, 0, 1);
		frdf.Print(" -2");
		int nn = elem->GetNumNodes();
		for (int j = 0; j < nn; j++)
			frdf.Print("%10d", elem->GetNode(j)->userId);
		if (!model.saveBreakout)
		{
			for (int j = 0; j < elem->GetNumLocalEdges(); j++)
			{
				int edgeNum = srToCgx[j];
				frdf.Print("%10d", elem->GetEdge(edgeNum)->GetMidNodeUserId());
				nn++;
				if (nn == 10 && needContinue)
				{
					frdf.PrintLine("\n");
					frdf.Print(" -2");
				}
			}
		}
		frdf.PrintLine("\n");
	}
	frdf.PrintLine(" -3");

	SRintVector tmpNodeDisps;
	tmpNodeDisps.Allocate(model.GetNumNodes());
	for (int f = 0; f < numFace; f++)
	{
		int fid = fidv[f];
		SRface* face = model.GetFace(fid);
		int nn = face->GetNumNodesTotal();
		if (model.saveBreakout)
			nn = face->GetNumNodes();
		for (int n = 0; n < nn; n++)
		{
			int nid = face->GetNodeOrMidNodeId(n);
			tmpNodeDisps.Put(nid, 1);
		}
	}
	//displacements

	frdf.PrintLine("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	frdf.PrintLine(" -4  DISP        4    1");
	frdf.PrintLine(" -5  D1          1    2    1    0");
	frdf.PrintLine(" -5  D2          1    2    2    0");
	frdf.PrintLine(" -5  D3          1    2    2    0");
	frdf.PrintLine(" -5  ALL         1    2    0    0    1ALL");

	SRvec3 disp;
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		disp.Zero();
		if (tmpNodeDisps.Get(i) == 1)
			disp.Assign(1.0, 1.0, 1.0);
		frdf.PrintLine(" -1%10d%12.5le%12.5le%12.5le", node->userId, disp.d[0], disp.d[1], disp.d[2]);
	}
	frdf.PrintLine(" -3");

	frdf.PrintLine("9999"); //end of data

	frdf.Close();

}

void SRpostProcess::PlotElems(int numelToPlot, int* elems, char* inName, int filenum)
{
	SRstring name;
	name.Copy(model.outdir);
	if (inName == NULL)
		name.Cat("\\plotelems");
	else
	{
		name.Cat("\\");
		name.Cat(inName);
	}
	if (filenum != -1)
	{
		char buf[20];
		sprintf_s(buf, "%d", filenum);
		name.Cat(buf);
	}
	name.Cat(".frd");
	SRfile frdf;
	frdf.Open(name, SRoutputMode);

	//mesh definition:
	int numNodeOut = 0;
	int brickSRtoCgx[12] = { 0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7 };
	int wedgeSRtoCgx[9] = { 0, 1, 2, 6, 7, 8, 3, 4, 5 };
	int tetSRtoCgx[6] = { 0, 1, 2, 3, 4, 5 };
	frdf.PrintLine("    1C"); // ''1C'' defines a new calc
	frdf.PrintLine("    1UDATE   26.march.2000"); // ''1U'' stores user job - information, can be any string, ttd put real date, add more lines, e.g. PGM stressrefine, model file path 
	//nodes:
	SRnode* node = NULL;
	int numnode = model.GetNumNodes();
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		numNodeOut++;
	}
	frdf.PrintLine("    2C                  %12d                                     1", numNodeOut);
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		//http://bconverged.com/calculix/doc/cgx/html/node167.html
		//long Format:(1X,'-1',I10,3E12.5)
		frdf.PrintLine(" -1%10d%12.5le%12.5le%12.5le", node->userId, node->pos.d[0], node->pos.d[1], node->pos.d[2]);
	}
	frdf.PrintLine(" -3");

	int numel = model.GetNumElements();
	frdf.PrintLine("    3C                  %12d                                     1", numel);
	int *srToCgx;
	for (int e = 0; e < numel; e++)
	{
		SRelement* elem = model.GetElement(e);
		int type;
		bool needContinue = false;
		if (elem->type == brick)
		{
			needContinue = true;
			type = 4;
			srToCgx = brickSRtoCgx;
		}
		else if (elem->type == wedge)
		{
			needContinue = true;
			type = 5;

			srToCgx = wedgeSRtoCgx;
		}
		else
		{
			type = 6;
			srToCgx = tetSRtoCgx;
		}
		frdf.PrintLine(" -1%10d%5d%5d%5d", elem->userId, type, 0, 1);
		frdf.Print(" -2");
		int nn = elem->GetNumNodes();
		for (int j = 0; j < nn; j++)
			frdf.Print("%10d", elem->GetNode(j)->userId);
		if (!model.saveBreakout)
		{
			for (int j = 0; j < elem->GetNumLocalEdges(); j++)
			{
				int edgeNum = srToCgx[j];
				frdf.Print("%10d", elem->GetEdge(edgeNum)->GetMidNodeUserId());
				nn++;
				if (nn == 10 && needContinue)
				{
					frdf.PrintLine("\n");
					frdf.Print(" -2");
				}
			}
		}
		frdf.PrintLine("\n");
	}
	frdf.PrintLine(" -3");

	SRintVector tmpNodeDisps;
	tmpNodeDisps.Allocate(model.GetNumNodes());
	for (int i = 0; i < numelToPlot; i++)
	{
		int eid = elems[i];
		SRelement* elem = model.GetElement(eid);
		for (int n = 0; n < elem->nodeIds.GetNum(); n++)
		{
			int nid = elem->nodeIds.d[n];
			tmpNodeDisps.Put(nid, 1);
		}
		if (!model.saveBreakout)
		{
			for (int e = 0; e < elem->localEdges.GetNum(); e++)
			{
				int mid = elem->GetEdge(e)->midnodeId;
				tmpNodeDisps.Put(mid, 1);
			}
		}
		int jjj = 0;
	}
	//displacements

	frdf.PrintLine("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	frdf.PrintLine(" -4  DISP        4    1");
	frdf.PrintLine(" -5  D1          1    2    1    0");
	frdf.PrintLine(" -5  D2          1    2    2    0");
	frdf.PrintLine(" -5  D3          1    2    2    0");
	frdf.PrintLine(" -5  ALL         1    2    0    0    1ALL");

	SRvec3 disp;
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		disp.Zero();
		if (tmpNodeDisps.Get(i) == 1)
			disp.Assign(1.0, 1.0, 1.0);
		frdf.PrintLine(" -1%10d%12.5le%12.5le%12.5le", node->userId, disp.d[0], disp.d[1], disp.d[2]);
	}
	frdf.PrintLine(" -3");

	frdf.PrintLine("9999"); //end of data

	frdf.Close();
}

void SRpostProcess::PlotSBElems(char* inName, int filenum)
{
	SRintVector ellist;
	ellist.Allocate(model.GetNumElements());
	int nelsbo = 0;
	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		if (elem->saveForBreakout)
		{
			ellist.Put(nelsbo, elem->id);
			nelsbo++;
		}
	}
	PlotElems(nelsbo, ellist.d, inName, filenum);
}

void SRpostProcess::PlotElemsSBOnly(int numelToPlot, int* elems, char* inName, int filenum)
{
	SRstring name;
	name.Copy(model.outdir);
	if (inName == NULL)
		name.Cat("\\plotelems");
	else
	{
		name.Cat("\\");
		name.Cat(inName);
	}
	if (filenum != -1)
	{
		char buf[20];
		sprintf_s(buf, "%d", filenum);
		name.Cat(buf);
	}
	name.Cat(".frd");
	SRfile frdf;
	frdf.Open(name, SRoutputMode);

	//mesh definition:
	int numNodeOut = 0;
	int brickSRtoCgx[12] = { 0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7 };
	int wedgeSRtoCgx[9] = { 0, 1, 2, 6, 7, 8, 3, 4, 5 };
	int tetSRtoCgx[6] = { 0, 1, 2, 3, 4, 5 };
	frdf.PrintLine("    1C"); // ''1C'' defines a new calc
	frdf.PrintLine("    1UDATE   26.march.2000"); // ''1U'' stores user job - information, can be any string, ttd put real date, add more lines, e.g. PGM stressrefine, model file path 
	//nodes:
	SRnode* node = NULL;
	int numnode = model.GetNumNodes();
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		numNodeOut++;
	}
	frdf.PrintLine("    2C                  %12d                                     1", numNodeOut);
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		//http://bconverged.com/calculix/doc/cgx/html/node167.html
		//long Format:(1X,'-1',I10,3E12.5)
		frdf.PrintLine(" -1%10d%12.5le%12.5le%12.5le", node->userId, node->pos.d[0], node->pos.d[1], node->pos.d[2]);
	}
	frdf.PrintLine(" -3");

	int numel = model.GetNumElements();
	int nelout = 0;
	for (int e = 0; e < numel; e++)
	{
		SRelement* elem = model.GetElement(e);
		if (elem->saveForBreakout)
			nelout++;
	}

	frdf.PrintLine("    3C                  %12d                                     1", nelout);
	int *srToCgx;
	for (int e = 0; e < numel; e++)
	{
		SRelement* elem = model.GetElement(e);
		if (!elem->saveForBreakout)
			continue;
		int type;
		bool needContinue = false;
		if (elem->type == brick)
		{
			needContinue = true;
			type = 4;
			srToCgx = brickSRtoCgx;
		}
		else if (elem->type == wedge)
		{
			needContinue = true;
			type = 5;

			srToCgx = wedgeSRtoCgx;
		}
		else
		{
			type = 6;
			srToCgx = tetSRtoCgx;
		}
		frdf.PrintLine(" -1%10d%5d%5d%5d", elem->userId, type, 0, 1);
		frdf.Print(" -2");
		int nn = elem->GetNumNodes();
		for (int j = 0; j < nn; j++)
			frdf.Print("%10d", elem->GetNode(j)->userId);
		if (!model.saveBreakout)
		{
			for (int j = 0; j < elem->GetNumLocalEdges(); j++)
			{
				int edgeNum = srToCgx[j];
				frdf.Print("%10d", elem->GetEdge(edgeNum)->GetMidNodeUserId());
				nn++;
				if (nn == 10 && needContinue)
				{
					frdf.PrintLine("\n");
					frdf.Print(" -2");
				}
			}
		}
		frdf.PrintLine("\n");
	}
	frdf.PrintLine(" -3");

	SRintVector tmpNodeDisps;
	tmpNodeDisps.Allocate(model.GetNumNodes());
	for (int i = 0; i < numelToPlot; i++)
	{
		int eid = elems[i];
		SRelement* elem = model.GetElement(eid);
		for (int n = 0; n < elem->nodeIds.GetNum(); n++)
		{
			int nid = elem->nodeIds.d[n];
			tmpNodeDisps.Put(nid, 1);
		}
		if (!model.saveBreakout)
		{
			for (int e = 0; e < elem->localEdges.GetNum(); e++)
			{
				int mid = elem->GetEdge(e)->midnodeId;
				tmpNodeDisps.Put(mid, 1);
			}
		}
		int jjj = 0;
	}
	//displacements

	frdf.PrintLine("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	frdf.PrintLine(" -4  DISP        4    1");
	frdf.PrintLine(" -5  D1          1    2    1    0");
	frdf.PrintLine(" -5  D2          1    2    2    0");
	frdf.PrintLine(" -5  D3          1    2    2    0");
	frdf.PrintLine(" -5  ALL         1    2    0    0    1ALL");

	SRvec3 disp;
	for (int i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		disp.Zero();
		if (tmpNodeDisps.Get(i) == 1)
			disp.Assign(1.0, 1.0, 1.0);
		frdf.PrintLine(" -1%10d%12.5le%12.5le%12.5le", node->userId, disp.d[0], disp.d[1], disp.d[2]);
	}
	frdf.PrintLine(" -3");

	frdf.PrintLine("9999"); //end of data

	frdf.Close();
}


void SRpostProcess::PlotSacr()
{
	SRintVector elList;
	int ns = 0;
	for (int i = 0; i < model.GetNumElements(); i++)
	{
		SRelement* elem = model.GetElement(i);
		if (elem->isSacrificial())
		{
			ns++;
			elList.PushBack(i);
		}
	}
	if (ns > 0)
		PlotElems(ns, elList.d, "sacrEls");
}

void SRpostProcess::clipBreakoutSacrStreses()
{
	if (model.breakout)
	{
		//don't let nodalStress of any sacrificial node be high compared to max in model
		double target = 0.8*model.GetStressMax();
		double stress[6], svmClipped;

		for (int i = 0; i < model.GetNumNodes(); i++)
		{
			SRnode* node = model.GetNode(i);
			if (!node->checkSacr())
				continue;
			for (int c = 0; c < 6; c++)
				stress[c] = nodalStress.Get(i, c);
			double svm = model.math.GetSvm(stress);
			if (svm > target)
			{
				double ratvm = target / svm;
				for (int c = 0; c < 6; c++)
				{
					//ratio down each component:
					nodalStress.Put(i, c, ratvm*stress[c]);
				}
			}
		}
	}
}

