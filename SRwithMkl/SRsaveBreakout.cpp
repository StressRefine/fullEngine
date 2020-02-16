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
// SRsaveBreakout.cpp: implementation of savebreakout routine of the SRpostProcess class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <stdlib.h>
#include <search.h>
#include "SRmodel.h"
#include "SRelement.h"

extern SRmodel model;

static int SRElemWithRadiusCompareFunc(const void *v1, const void *v2)
{
	//compare function for sorting and searching of list elements
	//input:
	//v1 = pointer to elem  1
	//v21 = pointer to elem  2
	//return:
	//-1 if node1 elem 1's centroid is less than elem 2 from pos
	//0 if node1 elem 1's centroid is eqidistant than elem 2 from pos
	//+1 if node1 elem 1's centroid is greater than elem 2 from pos
	SRElemWithRadius* elc1 = (SRElemWithRadius *)v1;
	SRElemWithRadius* elc2 = (SRElemWithRadius *)v2;

	if (elc1->Rad < elc2->Rad)
		return -1;
	else if (elc1->Rad == elc2->Rad)
		return 0;
	else
		return 1;
}


int SRpostProcess::findSaveBreakoutElementsTopo(int nel, int& nelAllNodesDisp, bool& topoIsRagged)
{
	//for model that needs topo filtering. do topo filtering first. don't do autobreakoutsphere before topo filtering to trim
	//model unless it is very large. The only reason for trimming is so createelemedges and fillglobalfaces is not slow.
	//I timed them for full model, multithinslot5, 296K elements and it took 0.018 secs. So use 300K as trim threshhold.
	SRBreakoutData *bdat = &model.saveBreakoutData;
	int nelMaxNeedTrim = 300000;
	int nelBreakoutMaxTopo = 10000; //this is larger than for single-part models to make it more likely that entire part
									//is in breakout model. Testing shows this will still solve quickly
	int nelnonsacr = MATHMAX(5000, (2 * nelBreakoutMaxTopo) / 3);
	bool trimNeeded = (nel > nelMaxNeedTrim);
	int nelClipped = nel;
	if (bdat->fromPartialDispFile)
	{
		nelAllNodesDisp = 0;
		//temporarily set the breakout origin to the centroid of the displacement data table for autobreakoutsphere;
		//it will be reset to the one estimated from fea stresses after topofiltering
		if (bdat->atMax)
			nelAllNodesDisp = fillBreakoutModelOriginFromDispCentroid(bdat);
	}
	if (trimNeeded)
	{
		nelClipped = autoBreakoutSphere(bdat, nelMaxNeedTrim);
		//PlotElemsSBOnly(1, &elIdAtRad0, "elidAtRad0");
	}
	else
		findElidAtRad0andSetSB(bdat);

	int nelBreakout = 0;

	if (bdat->fromPartialDispFile && nelAllNodesDisp < minElForBreakout)
	{
		REPPRINT("Error: Too small of a region was chosen in your Data Table to achieve accurate stress results");
		REPPRINT("in a breakout analyis. Please choose a larger region");
		ERROREXIT;
	}

	model.edges.Free();
	model.edges.Allocate(12 * nelClipped);

	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		if (elem->saveForBreakout)
			model.CreateElemEdgesSBO(elem);
	}

	int needGlobalNodeOrderForBreakout = bdat->fromPartialDispFile;
	model.FillGlobalFaces(needGlobalNodeOrderForBreakout);

	bool elidAtRad0FromStress = false;
	if (bdat->fromPartialDispFile)
	{
		//fillBreakoutModelOriginFromDispCentroid can be inadequate- see model thinslotbox, it broke out a part
		//that was not at max stress.
		//but fillBreakoutModelOriginFromHCentroidStresses is not safe, 
		//see multithinslot5 twopart version, found a wierd singular place at
		//edge of data table. two remedies:
		//1. since we've filled faces, there is enough info to catch reentrant corner.
		//   make those sacrificial and don't search them for stresses:
		findSacrificialElementsAtKinks();

		//2. only look for max stress in the vicinitiy of the centroid. already sorted elements
		//   by radius from centroid, search just 1st nel/2.
		int nsb = nel / 2;
		fillBreakoutModelOriginFromHCentroidStressesElemRadList(bdat, nsb);
		PlotElems(1, &elIdAtRad0, "elidatRad0fromHstresses");
		elidAtRad0FromStress = true;
		//sacrificial status shouln't matter during saveBreakout run but 
		//clean it up for neatness:
		for (int i = 0; i < model.GetNumElements(); i++)
		{
			SRelement* elem = model.GetElement(i);
			if (!elem->saveForBreakout)
				continue;
			elem->sacrificial = 0;
		}
	}

	topoIsRagged = false;
	//do topo filtering with a higher allowable number of elements to make it more likely to find 
	//entire part
	nelBreakoutMax = nelClipped;
	nelBreakout = topoFilterSaveBreakoutElems(topoIsRagged);
	nelBreakoutMax = nelBreakoutMaxTopo;

	bool did2ndAutoSphere = false;
	if (nelBreakout < minElForBreakout)
	{
		REPPRINT("Warning: Breakout Model is small for accurate stress calculation.");
		REPPRINT("         Please verify that was your intent, for example it is ok");
		REPPRINT("         if the local region is a small part in an assembly");
	}
	else if (nelBreakout > nelBreakoutMax)
	{
		if (bdat->atMax && bdat->fromPartialDispFile)
		{
			if (!elidAtRad0FromStress)
			{
				//reset origin to elem with max stress determined from element centroidal stresses:
				fillBreakoutModelOriginFromHCentroidStresses(bdat, true);//SBonlyNoBsurfNoSacr is true
				//PlotElems(1, &elIdAtRad0, "BdatOriginFromHCentroidStresses");
				elidAtRad0FromStress = true;
			}
			//another pass with autosphere to trim model:
			for (int e = 0; e < model.GetNumElements(); e++)
			{
				SRelement* elem = model.GetElement(e);
				elem->saveForBreakout = false;
			}
			nelBreakout = autoBreakoutSphere(bdat, nelBreakoutMax);
			//PlotSBElems("AfterAutoBreakoutSphere");
			//even if entire part had been identified, this trimming pass will make breakout
			//ragged again:
			topoIsRagged = true;
			did2ndAutoSphere = true;
		}
	}

	if (topoIsRagged)
	{
		//only let the first nelBreakoutMaxNoTopo savebreakout elements be nonsacrificial.
		//This will speed things up and make it less likely to have hot spots near edge of breakout region
		bktElSacrList.Allocate(nelBreakout);
		bktElSacrList.Set(-1);
		SRintVector bktelSacrElIds;
		bktelSacrElIds.Allocate(nelBreakout);
		if (elemRadiusList.GetNum() == 0)
		{
			elemRadiusList.Allocate(nel);
			SRvec3 pos;
			SRvec3* posOrigin = &bdat->origin;
			int elemListLength = nel;

			for (int i = 0; i < nel; i++)
			{
				SRelement* elem = model.GetElement(i);
				elem->saveForBreakout = false;
				SRElemWithRadius& elc = elemRadiusList.Get(i);
				elc.Rad = elem->minCornerDistance(*posOrigin);
				elc.elId = i;
			}

			SortElemsByRadius();

			elIdAtRad0 = elemRadiusList.Get(0).elId;
		}

		int nelsbo = 0;
		if (nelBreakout > nelnonsacr)
		{
			for (int e = 0; e < elemRadiusList.GetNum(); e++)
			{
				SRElemWithRadius elwr = elemRadiusList.Get(e);
				int eid = elwr.elId;
				SRelement *elem = model.GetElement(eid);
				if (!elem->saveForBreakout)
					continue;
				if (nelsbo > nelnonsacr)
				{
					bktElSacrList.Put(nelsbo, elem->userId);
					bktelSacrElIds.Put(nelsbo, eid);
				}
				nelsbo++;
				if (nelsbo >= nelBreakout)
					break;
			}
		}

		//PlotElemsSBOnly(bktelSacrElIds.GetNum(), bktelSacrElIds.d, "bktElSacrList");
	}
	else if (!elidAtRad0FromStress && bdat->atMax && bdat->fromPartialDispFile)
	{
		fillBreakoutModelOriginFromHCentroidStresses(bdat);
		findElidAtRad0AndMaxRadNearOrigin(bdat);
		//PlotElemsSBOnly(1, &elIdAtRad0, "elIdAtRad0AfterHCentroidStresses");
	}

	elemOnTopoList.Free();
	return nelBreakout;
}

void SRpostProcess::saveBreakoutModel()
{
	int nel = model.GetNumElements();
	needTopoFilter = model.anyBsurf && !noTopoFilter;
	maxRadNearOrigin = model.size;
	SRBreakoutData *bdat = &model.saveBreakoutData;
	int nelBreakout = 0;
	if (bdat->bktNode != -1)
	{
		int nid = model.input.NodeFind(bdat->bktNode);
		if (nid == -1)
		{
			if (bdat->fromPartialDispFile)
				REPPRINT("Error: node specified for local region is not within the region for which displacement table was output");
			else
				REPPRINT("Error: node specified for local region is not in model");
			ERROREXIT;
		}
		SRnode* node = model.GetNode(nid);
		if (node != NULL)
		{
			if (bdat->fromPartialDispFile && !node->hasDisp)
			{
				REPPRINT("Error: Node specified for breakout region does not exist in displacement results file");
				ERROREXIT;
			}
			bdat->origin.Copy(node->pos);
		}
	}
	//if breakout at max, origin is determined from position of max from fea solution (feavmmaxpos)
	//during readPreviousDisplacementForBreakout or from eluidatmax in readResultsSrr.
	//if fromPartialDispFile stress from fea solution is not known. it will be estimated
	//using fillBreakoutModelOriginFromHCentroidStresses below
	 
	int nnode = model.GetNumNodes();
	model.input.numNodeEdges.Allocate(nnode);
	model.input.nodeEdges.Allocate(nnode, 20);
	model.input.numNodeFaces.Allocate(nnode);
	model.input.nodeFaces.Allocate(nnode, MAXNODEFACEOWNERS);
	if(model.ReadDispStressSRR)
		SCREENPRINT("reading FEA displacements\n");
	readPreviousDisplacementForBreakout();

	SCREENPRINT("determining breakout region\n");

	minElForBreakout = 1000;
	nelNearOrigin = 500;
	bool topoIsRagged = false;

	if (needTopoFilter)
		nelBreakout = findSaveBreakoutElementsTopo(nel, nelAllNodesDisp, topoIsRagged);
	else
	{
		nelBreakoutMax = 5000;


		if (nelBreakoutMax > nel)
			nelBreakoutMax = nel;


		SRvec3 localPos;
		double localRad;
		SRvec3 bktorigin;

		nelAllNodesDisp = 0;
		if (bdat->fromPartialDispFile && bdat->atMax)
			nelAllNodesDisp = fillBreakoutModelOriginFromHCentroidStresses(bdat);
		else
		{
			for (int e = 0; e < nel; e++)
			{
				SRelement* elem = model.GetElement(e);
				elem->saveForBreakout = false;
				if (elem->allNodesHaveDisp)
					nelAllNodesDisp++;
			}
		}

		if (nelAllNodesDisp < minElForBreakout)
		{
			REPPRINT("Error: Too small of a region was chosen in your Data Table to achieve accurate stress results");
			REPPRINT("in a breakout analyis. Please choose a larger region");
			ERROREXIT;
		}

		nelBreakout = autoBreakoutSphere(bdat, nelBreakoutMax);
		//PlotSBElems("SBelementsAfterAutoBreakoutSphere");

		model.edges.Free();
		model.edges.Allocate(12 * nelBreakout);

		for (int e = 0; e < model.GetNumElements(); e++)
		{
			SRelement* elem = model.GetElement(e);
			if (elem->saveForBreakout)
				model.CreateElemEdgesSBO(elem);
		}

		model.FillGlobalFaces();
	}


	elemRadiusList.Free();

	if (model.numShellOrBeamNodes != 0)
		model.input.applyShellOrBeamBreakoutCons();

	SRnode* node;
	//set all nodes that are not owned by the saved elements as orphans. They will not be output:
	for (int n = 0; n < model.GetNumNodes(); n++)
	{
		node = model.GetNode(n);
		node->SetFirstElementOwner(-1);
	}

	int nelSbkt = 0;
	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		if (elem->saveForBreakout)
		{
			nelSbkt++;
			for (int n = 0; n < elem->GetNumNodes(); n++)
			{
				node = elem->GetNode(n);
				node->SetFirstElementOwner(e);
			}
		}
	}

	model.FindElemsAdjacentToBreakout();

	int numfaces = model.GetNumFaces();
	breakoutConstraints.Allocate(numfaces);
	breakoutElemUids.Allocate(numfaces);
	int numNewConstraints = 0;
	int numBreakoutCon = 0;
	SRvec3 enfDisp;
	for (int f = 0; f < numfaces; f++)
	{
		SRface* face = model.GetFace(f);
		bool elbreakout0 = false, elbreakout1 = false;
		SRelement *elem = model.GetElement(face->GetElementOwner(0));
		SRelement *elem1 = NULL;
		if (elem->saveForBreakout)
			elbreakout0 = true;
		if (!face->IsBoundaryFace())
		{
			elem1 = model.GetElement(face->GetElementOwner(1));
			if (elem1->saveForBreakout)
				elbreakout1 = true;
		}

		bool needBreakoutCon = false;
		if (elbreakout0)
		{
			face->SetSaveForBreakout(0);
			if (!elbreakout1 && !face->IsBoundaryFace())
				needBreakoutCon = true;
		}
		if (elbreakout1)
		{
			elem = elem1;
			face->SetSaveForBreakout(1);
			if (!elbreakout0)
				needBreakoutCon = true;
		}

		if (needBreakoutCon)
		{
			breakoutElemUids.Put(numBreakoutCon, elem->GetUserid());
			numBreakoutCon++;
			SRconstraint* con = breakoutConstraints.Add();
			con->SetBreakout();
			con->SetType(faceCon);
			con->SetEntityId(f);
			int ncorner = face->GetNumNodes();
			int nnodestotal = 2 * ncorner;
			con->allocateEnforcedDisplacementData(nnodestotal);
			for (int n = 0; n < nnodestotal; n++)
			{
				int nid;
				if (n >= ncorner)
				{
					nid = face->GetMidNodeId(n - ncorner);
					node = model.GetNode(nid);
				}
				else
				{
					node = face->GetNode(n);
					nid = face->GetNodeId(n);
				}
				if (node->constraintId != -1)
				{
					//prevent a possible conflict
					SRconstraint* nodalcon = model.GetConstraint(node->constraintId);
					if (nodalcon->breakout)
					{
						nodalcon->type = inactiveCon;
						node->constraintId = -1;
					}
				}
				enfDisp.Copy(nodeDisps.Get(nid));
				for (int dof = 0; dof < 3; dof++)
				{
					con->SetConstrainedDof(dof);
					con->PutEnforcedDisplacementData(n, dof, enfDisp.d[dof]);
				}
			}
		}
	}

	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		SRconstraint* con = model.GetConstraint(i);
		if (con->type == nodalCon && con->breakout)
		{
			int nid = con->entityId;
			con->allocateEnforcedDisplacementData(1);
			enfDisp.Copy(nodeDisps.Get(nid));
			for (int dof = 0; dof < 3; dof++)
			{
				con->SetConstrainedDof(dof);
				con->PutEnforcedDisplacementData(0, dof, enfDisp.d[dof]);
			}
		}
	}
	SCREENPRINT("writing breakout model\n");

	outputBreakout(nelBreakout, numBreakoutCon);
}

void SRpostProcess::outputBreakout(int nelBreakout, int numBreakoutCon)
{

	SRBreakoutData *bdat = &model.saveBreakoutData;
	SRnode* node;
	OUTPRINT("Number of elements in breakout model: %d\n", nelBreakout);

	//output breakout model's settings file:
	SRstring line;
	bdat->saveBreakoutMshFile.filename.Left('\\', line);
	line += "\\report.txt";
	if (SRfile::Existcheck(line))
		SRfile::Delete(line.str);
	bdat->saveBreakoutMshFile.filename.Left('.', line);
	line += ".srs";
	SRfile bktSrs;
	bktSrs.Open(line, SRoutputMode);
	bktSrs.PrintLine("breakout %lg %d //maxradnearorigin, atnode, from model %s", maxRadNearOrigin, bdat->bktNode, model.fileNameTail.str);
	bktSrs.PrintLine("//Number of elements in Breakout model: %d", nelBreakout);
	int nbktElSacr = 0;
	for (int i = 0; i < bktElSacrList.GetNum(); i++)
	{
		if (bktElSacrList.Get(i) >= 0)
			nbktElSacr++;
	}
	if (nbktElSacr > 0)
		bktSrs.PrintLine("nbktElSacr %d", nbktElSacr);

	if (model.useUnits)
	{
		bktSrs.PrintLine("stress_conversion %lg %s", model.stressUnitConversion, model.stressUnitstr.str);
		bktSrs.PrintLine("length_conversion  %lg %s", model.lengthUnitConversion, model.lengthUnitstr.str);
	}
	else
		bktSrs.PrintLine("NOUNITS");
	if (!model.detectSacrificialElements)
		bktSrs.PrintLine("NOSACRIFICIAL");
	if (model.minPorder != 2 || model.input.customPmax != 0)
	{
		int maxp = MAXP;
		if (model.input.customPmax != 0)
			maxp = model.input.customPmax;
		bktSrs.PrintLine("PLIMITS %d %d", model.minPorder, maxp);
	}
	if (model.adaptLoopMax != 3)
		bktSrs.PrintLine("MAXITERATIONS %d", model.adaptLoopMax);
	if (model.ErrorTolerance != 0.05)
		bktSrs.PrintLine("ERRORTOL %lg", model.ErrorTolerance);
	if (model.errorChecker.lowStressTol != LOWSTRESSTOL)
		bktSrs.PrintLine("LOWSTRESSTOL %lg", model.errorChecker.lowStressTol);
	if (model.outputf06)
		bktSrs.PrintLine("outputf06");

	bktSrs.Close();


	//output breakout model:
	SRfile& bf = bdat->saveBreakoutMshFile;
	bf.Open(SRoutputMode);
	bf.PrintLine("breakout //from model %s", model.fileNameTail.str);
	bf.PrintLine("//Number of elements in Breakout model: %d", nelBreakout);

	bf.PrintLine("nodes");
	for (int n = 0; n < model.GetNumNodes(); n++)
	{
		node = model.GetNode(n);
		if (!node->isOrphan())
			bf.PrintLine(" %d %lg %lg %lg", node->GetUserid(), node->Position().d[0], node->Position().d[1], node->Position().d[2]);
	}
	bf.PrintLine("end nodes");

	bf.PrintLine("elements");
	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		if (elem->saveForBreakout)
		{
			bf.Print("%d ", elem->GetUserid());
			bf.Print(" %s", elem->GetMaterial()->GetName());
			for (int n = 0; n < elem->GetNumNodes(); n++)
				bf.Print(" %d", elem->GetNode(n)->GetUserid());
			bf.PrintReturn();
		}
	}
	bf.PrintLine("end elements");

	bf.PrintLine("constraints");
	int ncon = model.GetNumConstraints();
	for (int c = 0; c < ncon; c++)
	{
		SRconstraint* con = model.GetConstraint(c);
		if (con->isBreakout())
			continue;// breakout constraints handled below
		int eid = con->GetEntityId();
		if (eid == -1)
			continue;
		if (con->GetType() == nodalCon)
		{
			node = model.GetNode(eid);
			if (node->isOrphan())
				continue;
			bf.Print(" %d", node->GetUserid());
			for (int dof = 0; dof < 3; dof++)
			{
				if (con->IsConstrainedDof(dof))
				{
					double enfd = con->GetEnforcedDisp(0, dof);
					bf.Print(" %lg", enfd);
				}
				else bf.Print(" -");
			}
			if (!con->isGcs())
				bf.Print(" coord %s", con->GetCoord()->GetName());
			bf.PrintReturn();
		}
		else if (con->GetType() == faceCon)
		{
			SRface *face = model.GetFace(eid);
			if (!face->isSaveForBreakout())
				continue;
			int ncorner = face->GetNumNodes();
			int nnodestotal = ncorner + face->GetNumLocalEdges();
			for (int n = 0; n < nnodestotal; n++)
			{
				int nid;
				if (n >= ncorner)
				{
					nid = face->GetMidNodeId(n - ncorner);
					node = model.GetNode(nid);
				}
				else
				{
					node = face->GetNode(n);
					nid = face->GetNodeId(n);
				}
				if (node->isOrphan())//this can't happen
					continue;
				bf.Print(" %d", node->GetUserid());
				for (int dof = 0; dof < 3; dof++)
				{
					if (con->IsConstrainedDof(dof))
					{
						if (con->hasEnforcedDisp())
						{
							double enfd = con->GetEnforcedDisp(c, dof);
							bf.Print(" %lg", enfd);
						}
						else
							bf.Print(" 0");
					}
					else bf.Print(" -");
				}
				bf.PrintReturn();
			}
		}
	}

	bf.PrintLine("end constraints");
	if (nbktElSacr > 0)
	{
		bf.PrintLine("bktElSacrList");
		for (int i = 0; i < bktElSacrList.GetNum(); i++)
		{
			int euid = bktElSacrList.Get(i);
			if (euid >= 0)
				bf.PrintLine("%d", euid);
		}
		bf.PrintLine("end bktElSacrList");
	}

	if (model.input.numUnSupportedNode != 0)
	{
		bf.PrintLine("nodalbreakoutconstraints");
		int ncon = model.GetNumConstraints();
		for (int c = 0; c < ncon; c++)
		{
			SRconstraint* con = model.GetConstraint(c);
			if (!con->isBreakout())
				continue;
			int eid = con->GetEntityId();
			if (eid == -1)
				continue;
			if (con->GetType() == nodalCon)
			{
				node = model.GetNode(eid);
				if (node->isOrphan())
					continue;
				bf.Print(" %d", node->GetUserid());
				for (int dof = 0; dof < 3; dof++)
				{
					double enfd = con->GetEnforcedDisp(0, dof);
					bf.Print(" %lg", enfd);
				}
				bf.PrintReturn();
			}
		}

		bf.PrintLine("end nodalbreakoutconstraints");
	}

	if (numBreakoutCon > 0)
	{
		bf.PrintLine("breakoutConstraints");
		for (int i = 0; i < numBreakoutCon; i++)
		{
			SRconstraint* con = breakoutConstraints.GetPointer(i);
			int fid = con->GetEntityId();
			SRface* face = model.GetFace(fid);
			bf.Print(" %d", breakoutElemUids.Get(i));
			int nnodes = face->GetNumNodes();
			for (int n = 0; n < nnodes; n++)
				bf.Print(" %d", face->GetNode(n)->GetUserid());
			for (int l = 0; l < face->GetNumLocalEdges(); l++)
			{
				SRnode* node = model.GetNode(face->GetMidNodeId(l));
				bf.Print(" %d", node->GetUserid());
			}
			bf.PrintReturn();
			int numenfd = con->GetNumEnforcedDisp();
			for (int n = 0; n < numenfd; n++)
			{
				bf.Print("#"); //continue character
				for (int dof = 0; dof < 3; dof++)
					bf.Print(" %lg", con->getDisp(n, dof));
				bf.PrintReturn();
			}
		}
		bf.PrintLine("end breakoutConstraints");
	}
	breakoutConstraints.Free();
	breakoutElemUids.Free();


	int nfacePressure = 0;
	int nfaceTraction = 0;
	bf.PrintLine("forces");
	for (int i = 0; i < model.GetNumForces(); i++)
	{
		SRforce* force = model.GetForce(i);
		int eid = force->GetEntityId();
		if (force->GetType() == nodalForce)
		{
			if (eid == -1)
				continue;
			node = model.GetNode(eid);
			if (node->isOrphan())
				continue;
			printNodalForce(node->GetUserid(), force, bf);
		}
		else if (force->GetType() == faceForce) //handled separately as facePressures or faceTractions
		{
			int fid = -1;
			SRelement* elem = model.GetElement(force->elemId);
			if (!elem->saveForBreakout)
				continue;
			if (force->nv[1] == -1)
				fid = elem->GetFace(force->nv[0])->GetId();
			else
			{
				int gno[4];
				fid = model.input.elemFaceFind(elem, force->nv, gno);
			}
			force->entityId = fid;

			if (!model.GetFace(fid)->isSaveForBreakout())
				continue;
			if (force->isPressure())
				nfacePressure++;
			else
				nfaceTraction++;
			continue;
		}
	}
	bf.PrintLine("end forces");

	if (nfacePressure != 0)
	{
		bf.PrintLine("facePressures");
		for (int i = 0; i < model.GetNumForces(); i++)
		{
			SRforce* force = model.GetForce(i);

			if (!force->isPressure())
				continue;
			if (force->GetType() == faceForce)
			{
				int fid = force->GetEntityId();
				if (fid == -1)
					continue;
				SRface *face = model.GetFace(fid);
				if (!face->isSaveForBreakout())
					continue;
				int eid = face->GetElementOwner(0);
				int euid = model.GetElement(eid)->GetUserid();
				bf.Print(" %d", euid);
				int nnodes = face->GetNumNodes();
				for (int n = 0; n < nnodes; n++)
					bf.Print(" %d", face->GetNode(n)->GetUserid());
				if (nnodes == 3)
					bf.Print(" -1");
				for (int n = 0; n < nnodes; n++)
					bf.Print(" %lg", force->GetForceVal(n, 0));
				bf.PrintReturn();
			}
		}
		bf.PrintLine("end facePressures");
	}

	if (nfaceTraction != 0)
	{
		bf.PrintLine("faceTractions");
		for (int i = 0; i < model.GetNumForces(); i++)
		{
			SRforce* force = model.GetForce(i);
			int fid = force->GetEntityId();
			if (fid == -1)
				continue;
			if (force->isPressure())
				continue;
			if (force->GetType() == faceForce)
			{
				int fid = force->GetEntityId();
				SRface *face = model.GetFace(fid);
				if (!face->isSaveForBreakout())
					continue;
				int eid = face->GetElementOwner(0);
				int euid = model.GetElement(eid)->GetUserid();
				bf.Print(" %d", euid);
				int nnodes = face->GetNumNodes();
				for (int n = 0; n < nnodes; n++)
					bf.Print(" %d", face->GetNode(n)->GetUserid());
				if (nnodes == 3)
					bf.Print(" -1");
				bf.PrintReturn();
				for (int dof = 0; dof < 3; dof++)
				{
					bf.Print("#"); //continuation character
					for (int n = 0; n < nnodes; n++)
						bf.Print(" %lg", force->GetForceVal(n, dof));
					bf.PrintReturn();
				}
			}
		}
		bf.PrintLine("end faceTractions");
	}

	bf.PrintLine("materials");
	for (int i = 0; i < model.GetNumMaterials(); i++)
		model.GetMaterial(i)->printToFile(bf);
	bf.PrintLine("end materials");

	bf.PrintLine("coordinates");

	//skip the dummy default gcs coord system because it is automatically created:
	for (int i = 1; i < model.GetNumCoords(); i++)
	{
		SRcoord* coord = model.GetCoord(i);
		coord->PrintToFile(bf);
	}
	bf.PrintLine("end coordinates");

	if (model.volumeForces.GetNum() > 0)
	{
		bf.PrintLine("volumeForces");
		for (int i = 0; i < model.volumeForces.GetNum(); i++)
		{
			SRvolumeForce* vf = model.volumeForces.GetPointer(i);
			if (vf->type == gravity)
			{
				bf.Print("gravity");
				bf.PrintLine(" %lg %lg %lg", vf->g1, vf->g2, vf->g3);
			}
			else if (vf->type == centrifugal)
			{
				bf.Print("centrifugal");
				bf.PrintLine(" %lg %lg %lg %lg %lg %lg %lg %lg", sqrt(vf->omega2),
					vf->axis.d[0], vf->axis.d[1], vf->axis.d[2],
					vf->origin.d[0], vf->origin.d[1], vf->origin.d[2],
					vf->alpha);
			}
		}
		bf.PrintLine("end volumeForces");
	}
	if (model.thermalForce != NULL)
	{
		bf.PrintLine("thermalForce");
		if (model.thermalForce->constantTemp)
			bf.PrintLine("constant %lg", model.thermalForce->temp);
		else
		{
			bf.PrintLine("variable");
			for (int j = 0; j < model.thermalForce->nodalTemp.GetNum(); j++)
				bf.Print(" %d %lg", model.GetNode(j)->userId, model.thermalForce->nodalTemp.Get(j));
		}
		bf.PrintLine("end thermalForce");
	}


	bf.Close();

	model.input.numNodeEdges.Free();
	model.input.nodeEdges.Free();
	model.input.numNodeFaces.Free();
	model.input.nodeFaces.Free();
	nodeDisps.Free();
	model.sumFile.Close();
}

void SRpostProcess::readPreviousDisplacementForBreakout()
{
	if (!nodeDisps.isEmpty())
		return;
	nodeDisps.Allocate(model.GetNumNodes());
	//read displacement results for the model from a previous run for use in creating breakout models
	if (model.ReadDispStressSRR)
		model.readResultsSrr();
	else
		readPreviousResultsF06();
	if (model.saveBreakout)
		model.setElemsAllNodesHaveDisp();

}

void SRpostProcess::SortElemsByRadius()
{
	int elemListLength = elemRadiusList.GetNum();
	qsort(elemRadiusList.GetVector(), elemListLength, sizeof(SRElemWithRadius), SRElemWithRadiusCompareFunc);
}

void SRpostProcess::findElidAtRad0AndMaxRadNearOrigin(SRBreakoutData *bdat)
{
	int nel = model.GetNumElements();
	elemRadiusList.Allocate(nel);
	SRvec3 pos;
	SRvec3* posOrigin = &bdat->origin;

	int nelsb = 0;
	for (int i = 0; i < nel; i++)
	{
		SRelement* elem = model.GetElement(i);
		if (!elem->saveForBreakout)
			continue;
		SRElemWithRadius& elc = elemRadiusList.Get(nelsb);
		elc.Rad = elem->minCornerDistance(*posOrigin);
		elc.elId = nelsb;
		nelsb++;
	}
	elemRadiusList.Allocate(nelsb);

	SortElemsByRadius();

	elIdAtRad0 = elemRadiusList.Get(0).elId;

	//for maxRadNearOrigin use centroid distance not mincorner in case of long slivery elements
	//(maxRadNearOrigin is used in SingStressCheck)

	maxRadNearOrigin = 0.0;
	for (int i = 0; i < nelNearOrigin; i++)
	{
		SRElemWithRadius& elc = elemRadiusList.Get(i);
		int elid = elc.elId;
		SRelement *elem = model.GetElement(elid);
		double rad = elem->centroidDistance(*posOrigin);
		if (rad > maxRadNearOrigin)
			maxRadNearOrigin = rad;
	}
}

void SRpostProcess::findElidAtRad0andSetSB(SRBreakoutData *bdat)
{
	//note: it's not necessary to sort by radius to find elIdAtRad0,
	//but sorted elemradlist is needed to safely search for stresses
	//using fillBreakoutModelOriginFromHCentroidStresses.
	int nel = model.GetNumElements();

	elemRadiusList.Allocate(nel);
	int elemListLength = nel;


	double minDist = BIG;
	for (int i = 0; i < nel; i++)
	{
		SRelement* elem = model.GetElement(i);
		if (bdat->fromPartialDispFile && !elem->allNodesHaveDisp)
			continue;
		elem->saveForBreakout = true;
		double dist = elem->minCornerDistance(bdat->origin);
		SRElemWithRadius& elc = elemRadiusList.Get(i);
		elc.Rad = dist;
		elc.elId = i;
	}

	SortElemsByRadius();
	elIdAtRad0 = elemRadiusList.Get(0).elId;
}

int SRpostProcess::autoBreakoutSphere(SRBreakoutData *bdat, int NumElBreakoutCandidate)
{
	//save elements within a sphere as breakout model.
	//sphere is automatically determined by working outwards in radius from point of max stress in model
	//until desired number of elements
	int nelBreakout = 0;


	int nel = model.GetNumElements();
	elemRadiusList.Allocate(nel);
	SRvec3 pos;
	SRvec3* posOrigin = &bdat->origin;
	int elemListLength = nel;

	for (int i = 0; i < nel; i++)
	{
		SRelement* elem = model.GetElement(i);
		elem->saveForBreakout = false;
		SRElemWithRadius& elc = elemRadiusList.Get(i);
		elc.Rad = elem->minCornerDistance(*posOrigin);
		elc.elId = i;
	}

	SortElemsByRadius();

	elIdAtRad0 = elemRadiusList.Get(0).elId;

	//for maxRadNearOrigin use centroid distance not mincorner in case of long slivery elements
	//(maxRadNearOrigin is used in SingStressCheck)
	maxRadNearOrigin = 0.0;

	maxRadSaveBreakout = 0.0;
	for (int i = 0; i < NumElBreakoutCandidate; i++)
	{
		SRElemWithRadius& elc = elemRadiusList.Get(i);
		int elid = elc.elId;
		SRelement *elem = model.GetElement(elid);
		double rad = elem->centroidDistance(*posOrigin);
		if (i < nelNearOrigin)
		{
			if (rad > maxRadNearOrigin)
				maxRadNearOrigin = rad;
		}
		if (elem->allNodesHaveDisp)
		{
			elem->saveForBreakout = true;
			if (rad > maxRadSaveBreakout)
				maxRadSaveBreakout = rad;
			nelBreakout++;
		}
	}

	return nelBreakout;
}

void SRpostProcess::autoLocalAdapt()
{
	//set all elements to sacrificial that are not near point of max stress in model
	SRvec3 origin = model.maxStressPos;
	int nel = model.GetNumElements();
	int nelLocal = 5000;
	if (nel < nelLocal)
		return;

	elemRadiusList.Allocate(nel);
	SRvec3 pos;
	int elemListLength = nel;

	for (int i = 0; i < nel; i++)
	{
		SRelement* elem = model.GetElement(i);
		elem->saveForBreakout = false;
		SRElemWithRadius& elc = elemRadiusList.Get(i);
		elc.Rad = elem->centroidDistance(origin);
		elc.elId = i;
	}
	SortElemsByRadius();

	for (int i = nelLocal; i < nel; i++)
	{
		SRElemWithRadius& elc = elemRadiusList.Get(i);
		int elid = elc.elId;
		SRelement *elem = model.GetElement(elid);
		elem->sacrificial = 1;
	}

	elemRadiusList.Free();
}

int SRpostProcess::topoFilterSaveBreakoutElems(bool& isRagged)
{
	//traverse model topologically starting at elem closest to bkt origin
	//stop when encounter element faces with onPartBdry or elements withou allnodeshavedisp
	//elements get set to saveforbreakout if inside the filter region.

	SRintVector incheckList;
	SRintVector outcheckList;

	isRagged = false;
	int numInEl = 0;
	int numElOnTopoList = 0;
	int nel = model.GetNumElements();
	elemOnTopoList.Allocate(nel);
	elemOnTopoList.Set(-1);
	incheckList.Allocate(nel);
	outcheckList.Allocate(nel);
	numInEl = 1;
	SRelement* elem = model.GetElement(elIdAtRad0);
	if (!elem->saveForBreakout)
	{

		REPPRINT("Error: unable to create breakout model using position of max stress");
		REPPRINT("Please specify a node for the center of the local region");
		ERROREXIT;
	}
	elemOnTopoList.Put(elIdAtRad0, 0);
	numElOnTopoList++;
	incheckList.Put(0, elIdAtRad0);
	int numLevels = 0;
	while (1)
	{
		int nextoutEl = 0;
		for (int e = 0; e < numInEl; e++)
		{
			SRelement* elem = model.GetElement(incheckList.Get(e));
			checkElemFacesForPartBdry(elem, nextoutEl, outcheckList, numLevels, numElOnTopoList);
			if (numElOnTopoList > nelBreakoutMax)
			{
				isRagged = true;
				break;
			}
		}
		numInEl = nextoutEl;
		if (numInEl == 0)
			break;
		for (int e = 0; e < numInEl; e++)
			incheckList.Put(e, outcheckList.Get(e));
		//if(numLevels%5 == 0)
			//PlotTopoLevelElems(numLevels);
		numLevels++;
	}
	//PlotTopoLevelElems(numLevels);

	int nelBreakout = 0;
	for (int e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		if (elemOnTopoList.Get(elem->id) >= 0)
		{
			elem->saveForBreakout = true;
			nelBreakout++;
		}
		else
			elem->saveForBreakout = false;
	}
	return nelBreakout;
}

void SRpostProcess::checkElemFacesForPartBdry(SRelement* elem, int& nextOutEl, SRintVector& outchecklist, int numLevels, int& numElOnTopoList)
{
	if (!elem->saveForBreakout)
		return;//this element did not pass autobreakoutsphere. It won't have corresponding faces 

	for (int i = 0; i < elem->localFaces.GetNum(); i++)
	{
		SRface* face = elem->GetFace(i);
		int elemid2 = face->elementOwners[1];
		if (elemid2 == -1)
			continue;
		int elemid1 = face->elementOwners[0];
		if (elemid1 != elem->id)
			elemid2 = elemid1;
		SRelement *elem2 = model.GetElement(elemid2);
		if (!elem2->saveForBreakout)
			continue;//this element did not pass autobreakoutsphere. It won't have corresponding faces
		if (!face->onPartBdry)
		{
			int eid2OnTopoList = elemOnTopoList.Get(elem2->id);
			if (eid2OnTopoList == -1)
			{
				elemOnTopoList.Put(elem2->id, numLevels);
				numElOnTopoList++;
				if (numElOnTopoList > nelBreakoutMax || nextOutEl > nelBreakoutMax)
				{
					nextOutEl = 0;
					return;
				}
				outchecklist.Put(nextOutEl, elem2->id);
				nextOutEl++;
			}
		}
	}
}

int SRpostProcess::fillBreakoutModelOriginFromHCentroidStressesElemRadList(SRBreakoutData *bdat, int nelToCheck)
{
	//find max element centroid svm in any element
	double svmmax = 0.0;
	SRelement* elatmax = NULL;

	for (int e = 0; e < nelToCheck; e++)
	{
		int elid = elemRadiusList.Get(e).elId;
		SRelement* elem = model.GetElement(elid);
		if (elem->sacrificial != 0)
			continue;
		if (elem->allNodesHaveDisp)
		{
			nelAllNodesDisp++;
			bool skipEl = false;
			//check if element has any node that touches a bsurf, if so skip:
			for (int n = 0; n < elem->GetNumCorners(); n++)
			{
				SRnode* node = elem->GetNode(n);
				if (node->sacrificialType == onBsurf)
				{
					skipEl = true;
					break;
				}
			}
			if (skipEl)
				continue;
			double svm = elem->CalculateSvmHCentroid();
			if (svm > svmmax)
			{
				elatmax = elem;
				svmmax = svm;
			}
		}
	}
	SRvec3 centroidPos;
	if (elatmax != NULL)
	{
		model.map.ElementCentroid(elatmax, centroidPos);
		elIdAtRad0 = elatmax->id;
	}
	bdat->origin.Copy(centroidPos);

	return nelAllNodesDisp;
}


int SRpostProcess::fillBreakoutModelOriginFromHCentroidStresses(SRBreakoutData *bdat, bool SBonlyNoBsurfNoSacr)
{
	//find max element centroid svm in any element
	double svmmax = 0.0;
	SRelement* elatmax = NULL;

	int nel = model.GetNumElements();

	for (int e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		if (SBonlyNoBsurfNoSacr && !elem->saveForBreakout || (elem->sacrificial != 0) )
			continue;
		if (elem->allNodesHaveDisp)
		{
			nelAllNodesDisp++;
			if (SBonlyNoBsurfNoSacr)
			{
				bool skipEl = false;
				//check if element has any node that touches a bsurf, if so skip:
				for (int n = 0; n < elem->GetNumCorners(); n++)
				{
					SRnode* node = elem->GetNode(n);
					if (node->sacrificialType == onBsurf)
					{
						skipEl = true;
						break;
					}
				}
				if (skipEl)
					continue;
			}
			double svm = elem->CalculateSvmHCentroid();
			if (svm > svmmax)
			{
				elatmax = elem;
				svmmax = svm;
			}
		}
	}
	SRvec3 centroidPos;
	if (elatmax != NULL)
	{
		model.map.ElementCentroid(elatmax, centroidPos);
		elIdAtRad0 = elatmax->id;
	}
	bdat->origin.Copy(centroidPos);

	return nelAllNodesDisp;
}

int SRpostProcess::fillBreakoutModelOriginFromDispCentroid(SRBreakoutData *bdat)
{
	double xcentroid = 0.0;
	double ycentroid = 0.0;
	double zcentroid = 0.0;

	int nelAllNodesDisp = 0;
	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		elem->saveForBreakout = false;
		elem->FillMappingNodesBreakoutOnly();
		if (elem->allNodesHaveDisp)
		{
			nelAllNodesDisp++;
			int nc = elem->GetNumCorners();
			double oneOverNc = 1.0 / nc;
			double elcentroidx = 0.0;
			double elcentroidy = 0.0;
			double elcentroidz = 0.0;
			for (int n = 0; n < nc; n++)
			{
				elcentroidx += elem->xnode[n];
				elcentroidy += elem->ynode[n];
				elcentroidz += elem->znode[n];
			}
			elcentroidx *= oneOverNc;
			elcentroidy *= oneOverNc;
			elcentroidz *= oneOverNc;
			xcentroid += elcentroidx;
			ycentroid += elcentroidy;
			zcentroid += elcentroidz;
		}
	}
	xcentroid /= nelAllNodesDisp;
	ycentroid /= nelAllNodesDisp;
	zcentroid /= nelAllNodesDisp;
	bdat->origin.Assign(xcentroid, ycentroid, zcentroid);
	return nelAllNodesDisp;
}

void SRpostProcess::PlotTopoLevelElems(int numLevels)
{
	bool dbgPlotTopoElems = false;//ttd!!
	if (!dbgPlotTopoElems)
		return;
	SRintVector ellist;
	for (int i = 0; i < model.elements.GetNum(); i++)
	{
		int level = elemOnTopoList.Get(i);
		if (level != -1 && level <= numLevels)
			ellist.PushBack(i);
	}
	int nelOnList = ellist.GetNum();
	PlotElems(nelOnList, ellist.d, "PlotTopoLevelElems", numLevels);
}

void SRpostProcess::PlotElemsByRadius(double Rad, int filenum)
{
	bool dbgPlotRadElems = false;//ttd!!
	if (!dbgPlotRadElems)
		return;
	SRintVector ellist;
	for (int i = 0; i < elemRadiusList.GetNum(); i++)
	{
		SRElemWithRadius& elwr = elemRadiusList.Get(i);
		double irad = elwr.Rad;
		if (irad <= Rad)
		{
			int eid = elwr.elId;
			if (model.GetElement(eid)->saveForBreakout)
				ellist.PushBack(eid);
		}
		else
			break;
	}
	PlotElems(ellist.GetNum(), ellist.d, "PlotElemsByRadius", filenum);
}

bool SRpostProcess::findSacrificialElementsAtKinks()
{
	//mark elements as sacrificial that have kinks at face intersections

	SRintVector singNode;
	int nnode = model.GetNumNodes();
	singNode.Allocate(nnode);
	SRintVector singEdge;
	singEdge.Allocate(model.GetNumEdges());
	bool anySingNode = false;
	bool anySingEdge = false;

	//check for reentrant corners faces on surface; singular unless
	//all 3 dofs of both faces contrained
	int nej = model.GetNumEdges();
	SRintVector edgeface1;
	SRintVector edgeface2;
	SRintVector edgefacelej1;
	SRintVector edgefacelej2;
	edgeface1.Allocate(nej);
	edgeface2.Allocate(nej);
	edgefacelej1.Allocate(nej);
	edgefacelej2.Allocate(nej);
	for (int i = 0; i < nej; i++)
	{
		edgeface1.Put(i, -1);
		edgeface2.Put(i, -1);
	}

	SRface* face;
	int n = model.GetNumFaces();
	for (int i = 0; i < n; i++)
	{
		face = model.GetFace(i);
		if (!face->IsBoundaryFace())
			continue;
		for (int lej = 0; lej < face->GetNumLocalEdges(); lej++)
		{
			int gej = face->GetLocalEdgeGlobalId(lej);
			if (edgeface1.Get(gej) == -1)
			{
				edgeface1.Put(gej, i);
				edgefacelej1.Put(gej, lej);
			}
			else
			{
				edgeface2.Put(gej, i);
				edgefacelej2.Put(gej, lej);
			}
		}
	}

	SRface* face2;
	double dotmin = 1.0;

	for (int i = 0; i < nej; i++)
	{
		SRedge* edge = model.GetEdge(i);

		int f1 = edgeface1.Get(i);
		int f2 = edgeface2.Get(i);
		if (f1 == -1 || f2 == -1)
			continue;
		face = model.GetFace(f1);
		face2 = model.GetFace(f2);
		int lej = edgefacelej1.Get(i);
		int lej2 = edgefacelej2.Get(i);

		SRvec3 norm, norm2;
		SRvec3 pos1, pos2, v21;
		double r, s;
		face->NaturalCoordinatesNearMidedge(lej, r, s);
		face->Position(r, s, pos1);
		face->OutwardNormal(r, s, norm);
		face2->NaturalCoordinatesNearMidedge(lej2, r, s);
		face2->Position(r, s, pos2);
		face2->OutwardNormal(r, s, norm2);
		pos1.Subtract(pos2, v21);
		v21.Normalize();
		double normdot = norm.Dot(norm2);
		if (normdot < 0.866)//150 degrees
		{
			//kink in the normals to the faces at the corners
			if ((norm2.Dot(v21) > 0.0) && (norm.Dot(v21) < 0.0))
			{
				//reentrant if vector from point on face2 to point on face 1
				//points in same direction as normal to face2 and in opposite direction to normal to face1
				//OUTPRINT(" sing edge face %d norm dot %lg\n", face->GetId(), dot);
				anySingEdge = true;
				if (normdot < dotmin)
				{
					dotmin = normdot;
					SRelement* elem;
					int eid = face->GetElementOwner(0);
					elem = model.GetElement(eid);
					eid = face2->GetElementOwner(0);
					elem = model.GetElement(eid);
				}
				if (normdot >= -0.1)
				{
					// < -0.1 is suspicious. could be crack-like feature, but could also be spurious, e.g. outward normal sign
					//fooled by slivery element.
					singEdge.Put(i, 1);
				}
			}
		}
	}
	if (!anySingEdge)
		return false;

	//mark nodes owned by singular edges as singular:
	for (int i = 0; i < model.GetNumEdges(); i++)
	{
		int pf = singEdge.Get(i);
		if (pf != 0)
		{
			SRedge* edge = model.GetEdge(i);
			for (int j = 0; j < 2; j++)
			{
				int id = edge->GetNodeId(j);
				singNode.Put(id, pf);
			}
		}
	}

	//mark all elements that own a singular node as sacrificial:

	int numsacr = 0;
	for (int i = 0; i < model.GetNumElements(); i++)
	{
		SRelement* elem = model.GetElement(i);
		if (!elem->saveForBreakout)
			continue;
		for (int j = 0; j < elem->GetNumNodes(); j++)
		{
			int id = elem->GetNodeId(j);
			int pf = singNode.Get(id);
			if (pf > 0)
			{
				elem->SetSacrificial(pf);
				numsacr++;
				break;
			}
		}
	}

	bool dbgPlotSacr = true;
	if (dbgPlotSacr)
		model.post.PlotSacr();

	return (numsacr != 0);
}


