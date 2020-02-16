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
// SRinputUtilities.cpp: implementation of utility routines for the SRinput class.
//
//////////////////////////////////////////////////////////////////////



#include "stdafx.h"
#include <stdlib.h>
#include <search.h>
#include "SRmodel.h"
#include "SRmachDep.h"
#include "SRinput.h"

extern SRmodel model;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

int SRuidDataCompareFunc(const void *v1, const void *v2)
{
	//compare function for sorting and searching of node userids
	//input:
		//v1 = pointer to uid data for node 1
		//v2 = pointer to uid data for node 2
	//return:
		//-1 if node1 uid is less than node2 uid
		//0  if node1 uid is equal to node2 uid
		//+1 if node1 uid is greater than node2 uid

	SRuidData* node1 = (SRuidData*)v1;
	SRuidData* node2 = (SRuidData*) v2;
	if (node1->uid < node2->uid)
		return -1;
	else if (node1->uid == node2->uid)
		return 0;
	else
		return 1;
}

void SRinput::Cleanup()
{
	nodeUids.Free();
	elemUids.Free();
}


int SRinput::NodeFind(int uid)
{
	//find node with user Id uid
	//input:
	//uid = user id to match
	//return:
	//number of the node that matches uid, -1 if not found

	//binary search:
	if (uid == -1)
		return -1;
	else if (nodeUidOffset != -1)
	{
		int nid = uid - nodeUidOffset;
		//make sure nid is in bounds:
		if (nid < 0 || nid >= model.GetNumNodes())
			return -1;
		return nid;
	}
	else
	{
		SRuidData* nuid;
		SRuidData uidt;
		//"search key" has to be same data type expected by compare function see SRuidDataCompareFunc:
		uidt.uid = uid;
		nuid = (SRuidData *)bsearch(&uidt, nodeUids.GetVector(), nodeUids.GetNum(), sizeof(SRuidData), SRuidDataCompareFunc);
		if (nuid == NULL)
			return -1;
		else
			return nuid->id;
	}
}


int SRinput::elemFind(int uid)
{
	//find elem with user Id uid
	//input:
		//uid = user id to match
	//return:
		//number of the node that matches uid, -1 if not found

	//binary search:
	if (elemUidOffset != -1)
		return uid - elemUidOffset;
	else
	{
		SRuidData* euid;
		SRuidData uidt;
		//"search key" has to be same data type expected by compare function see SRuidDataCompareFunc:
		uidt.uid = uid;
		euid = (SRuidData *)bsearch(&uidt, elemUids.GetVector(), elemUids.GetNum(), sizeof(SRuidData), SRuidDataCompareFunc);
		if (euid == NULL)
			return -1;
		else
			return euid->id;
	}
}

void SRinput::InputBreakoutSacrElems()
{
	SRstring line;
	int numelread = 0;
	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank(SKIPCONTINATION))
			continue;
		if (line == "end")
			break;
		int euid;
		line.TokRead(euid);
		model.post.bktElSacrList.Put(numelread, euid);
		numelread++;
		if (numelread >= model.post.bktElSacrList.GetNum())
			break;
	}
}

void SRinput::CountEntities(int &num)
{	
	//count the entities currently being read in input file by counting 
	//lines (except blank or comment) until end is encountered

	SRstring line;
	while (1)
	{
		if(!model.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank(SKIPCONTINATION))
            continue;
		if(line=="end")
			break;
		num++;
	}
}

void SRinput::CountElements(int &num)
{
	//count the elements currently being read in input file by counting 
	//lines (except blank or comment) until end is encountered
	//note:
		//this is the same as CountEntities but also checks for linear mesh

	SRstring line, tok;
	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (num == 0)
		{
			tok = line.Token();//skip uid field
			tok = line.Token();//skip material name field
			int nt = 0, nuid;
			while (1)
			{
				if (!line.TokRead(nuid))
					break;
				nt++;
			}
			if (nt == 4 || nt == 6 || nt == 8)
				ERROREXIT; //linear mesh not allowed
		}
		if (line.isCommentOrBlank(SKIPCONTINATION))
			continue;
		if (line == "end")
			break;
		num++;
	}
}

int SRinput::GetMaterialId(SRstring &name)
{
	//look up the material id with "name"
	//return:
		//material number, -1 if not found

	SRmaterial* mat;
	for (int i = 0; i < model.materials.GetNum(); i++)
	{
		mat = model.GetMaterial(i);
		if (name.Compare(mat->name))
			return i;
	}
	return -1;
}

int SRinput::GetCoordId(SRstring &name)
{
	//look up the coordinate system id with "name"
	//return:
		//coordinate system number, -1 if not found

	SRcoord*     coord;
	for (int i = 0; i < model.Coords.GetNum(); i++)
	{
		coord = model.GetCoord(i);
		if (name == coord->name)
			return i;
	}
	SRASSERT(false);
	return -1;
}
void SRinput::SortNodes()
{
	//sort "nodeUids" array in ascending order of uid
	int n = nodeUids.GetNum();
	qsort(nodeUids.GetVector(), n, sizeof(SRuidData), SRuidDataCompareFunc);
}

void SRinput::SortElems()
{
	//sort "elemUids" array in ascending order of uid
	int n = elemUids.GetNum();
	qsort(elemUids.GetVector(), n, sizeof(SRuidData), SRuidDataCompareFunc);
}

int SRinput::nodalToFaceConstraints(bool countOnly)
{
	//convert nodal to face constraints if all nodes of the face are constrained
	//input:
		//countOnly = true to just count the new face constraints, not add them
	//notes:
		//adds new face constraints to model.constraints
		//use of nodeconstore: duplicate nodal constraints are blended together if they are compatible (bot gcs or both lcs with same coord sys)
		//there will only be duplicate constraints in nodeconstore if they are not compatible.
		//So it's necessary to use nodeconstore to loop over all the duplicate constraints of each node

	int nface = model.GetNumFaces();

	int nfaceCon = 0;
	for (int f = 0; f < nface; f++)
	{
		SRface* face = model.GetFace(f);
		if (face->elementOwners[1] != -1)//not boundary face
			continue;

		int nn = face->GetNumNodes();

		bool allNodesConstrained = true;
		bool anyNodeConstrained = false;
		for (int n = 0; n < nn; n++)
		{
			int nid = face->GetNodeId(n);
			int ncon = nodeConStore.Get(nid).cons.GetNum();
			if (ncon == 0)
			{
				allNodesConstrained = false;
				break;
			}
		}

		if (allNodesConstrained)
		{
			//check midside nodes:
			for (int l = 0; l < nn; l++)
			{

				int mid = face->localEdges.Get(l).GetEdge()->GetMidNodeId();
				if (nodeConStore.Get(mid).cons.GetNum() == 0)
				{
					allNodesConstrained = false;
					break;
				}
			}
		}

		if (!allNodesConstrained)
			continue;

		//convert to face constraint only if dofs constrained are compatible and coordIds are compatible:
		bool constrainedDof0[3];
		int nnodes = nn;
		int nnodesTotal = 2 * nnodes;
		int nid = face->GetNodeId(0);
		int ncon0 = nodeConStore.Get(nid).cons.GetNum();
		SRconstraint* newFaceCon = NULL;

		for (int icon0 = 0; icon0 < ncon0; icon0++)
		{
			bool compatible = true;
			SRconstraint* nodeCon0 = nodeConStore.Get(nid).GetCon(icon0);
			int coord0 = nodeCon0->coordId;
			for (int dof = 0; dof < 3; dof++)
				constrainedDof0[dof] = nodeCon0->IsConstrainedDof(dof);
			bool allNodesConstrainedDofSame[3];
			bool anyNodeEnforcedDof[3];
			for (int dof = 0; dof < 3; dof++)
			{
				allNodesConstrainedDofSame[dof] = true;
				if (fabs(nodeCon0->getDisp(0, dof)) > TINY)
					anyNodeEnforcedDof[dof] = true;
				else
					anyNodeEnforcedDof[dof] = false;
			}
			for (int n = 1; n < nnodesTotal; n++)
			{
				int nid = face->GetNodeOrMidNodeId(n);
				int ncon = nodeConStore.Get(nid).cons.GetNum();
				for (int icon = 0; icon < ncon; icon++)
				{
					compatible = true;
					SRconstraint* nodeCon = nodeConStore.Get(nid).GetCon(icon);
					int coord = nodeCon->coordId;
					if (coord != coord0)
					{
						compatible = false;
						continue;
					}
					for (int dof = 0; dof < 3; dof++)
					{
						if (!constrainedDof0[dof] || (constrainedDof0[dof] && !nodeCon->IsConstrainedDof(dof)))
							allNodesConstrainedDofSame[dof] = false;
						if (fabs(nodeCon->getDisp(0, dof)) > TINY)
							anyNodeEnforcedDof[dof] = true;
					}
					if (!compatible)
						continue;
					bool anydofsame = false;
					for (int dof = 0; dof < 3; dof++)
					{
						if (allNodesConstrainedDofSame[dof])
						{
							anydofsame = true;
							break;
						}
					}

					if (!anydofsame)
						compatible = false;
					if (compatible)
						break;
				}  //endif: for (icon = 0; icon < node->GetNumConstraints(); icon++)
				if (!compatible)
					break;

			} // endif: for (int n = 1; n < nnodesTotal; n++)

			if (compatible)
			{
				//can't have more than one constraint on same face:
				if (newFaceCon != NULL)
					ERROREXIT;
				nfaceCon++;
				if (countOnly)
					break;
				face->constraintId = model.constraints.GetNum();
				newFaceCon = model.constraints.Add();
				newFaceCon->entityId = f;
				for (int dof = 0; dof < 3; dof++)
				{
					if (allNodesConstrainedDofSame[dof])
						newFaceCon->constrainedDof[dof] = 1;
				}
				newFaceCon->type = faceCon;
				newFaceCon->coordId = coord0;
				bool anyEnforcedDisp = false;
				for (int dof = 0; dof < 3; dof++)
				{
					if (anyNodeEnforcedDof[dof])
						anyEnforcedDisp = true;
				}

				if (anyEnforcedDisp)
				{
					newFaceCon->enforcedDisplacementData.Allocate(nnodesTotal, 3);
					for (int dof = 0; dof < 3; dof++)
					{
						if (!newFaceCon->IsConstrainedDof(dof))
							continue;
						for (int n = 0; n < nnodesTotal; n++)
						{
							if (anyNodeEnforcedDof[dof])
							{
								double enfd = nodeCon0->getDisp(0, dof);
								newFaceCon->enforcedDisplacementData.Put(n, dof, enfd);
							}
						}
					}
				}
			} //endif: if (compatible)
		} // endif: for (int icon0 = 0; icon0 < node->GetNumConstraints(); icon0++)
	} //endif: for (int f = 0; f < nface; f++)

	if (!countOnly)
	{
		//turn off constrained dofs of nodes involved in face constraints:
		for (int f = 0; f < nface; f++)
		{
			SRface* face = model.GetFace(f);
			if (face->elementOwners[1] != -1)//not boundary face
				continue;
			SRconstraint* faceCon = face->GetConstraint();
			if (faceCon == NULL)
				continue;

			int nnodesTotal = 2 * face->GetNumNodes();;

			for (int n = 0; n < nnodesTotal; n++)
			{
				int nid = face->GetNodeOrMidNodeId(n);
				int ncon = nodeConStore.Get(nid).cons.GetNum();
				for (int icon = 0; icon < ncon; icon++)
				{
					SRconstraint* nodeCon = nodeConStore.Get(nid).GetCon(icon);
					for (int dof = 0; dof < 3; dof++)
					{
						if (faceCon->IsConstrainedDof(dof))
							nodeCon->constrainedDof[dof] = 0;
					}
				}
			}
		}

		for (int n = 0; n < model.GetNumConstraints(); n++)
		{
			SRconstraint* nodeCon = model.GetConstraint(n);
			if (nodeCon->type == nodalCon)
			{
				bool anydofconstrained = false;
				for (int dof = 0; dof < 3; dof++)
				{
					if (nodeCon->IsConstrainedDof(dof))
					{
						anydofconstrained = true;
						break;
					}
				}
				if (!anydofconstrained)
					nodeCon->type = inactiveCon;
			}
		}
	}

	return nfaceCon;
}

int SRinput::nodalToFaceForces(bool countOnly)
{
	//convert nodal forces to face forces if all the nodes of a face are loaded compatibly
	//and it is a boundary face

	//input:
		//countonly = true to just count the number of new face forces else false
	//note:
		//adds new face forces to model.forces

	int nface = model.GetNumFaces();

	SRintVector compatibleForceNum;
	compatibleForceNum.Allocate(8);
	compatibleForceNum.Set(-1);

	int numFaceForce = 0;
	for (int f = 0; f < nface; f++)
	{
		SRface* face = model.GetFace(f);
		for (int i = 0; i < face->GetNumLocalEdges(); i++)
			int mid = face->localEdges.Get(i).GetEdge()->GetMidNodeId();
		if (face->GetElementOwner(1) != -1) //not boundary face
			continue;

		int nn = face->GetNumNodes();

		bool allNodesLoaded = true;
		for (int n = 0; n < nn; n++)
		{
			int nid = face->GetNodeId(n);
			int nforce = nodeForceStore.Get(nid).forces.GetNum();
			if (nforce == 0)
			{
				allNodesLoaded = false;
				break;
			}
		}
		if (allNodesLoaded)
		{
			//check also for midside nodes
			for (int l = 0; l < nn; l++)
			{

				int mid = face->localEdges.Get(l).GetEdge()->GetMidNodeId();
				int nforce = nodeForceStore.Get(mid).forces.GetNum();
				if (nforce == 0)
				{
					allNodesLoaded = false;
					break;
				}
			}

		}
		if (!allNodesLoaded)
			continue;

		//convert to face force only if coordIds and pressure flags are compatible:
		int nnodes = nn;
		int nnodesTotal = 2 * nnodes;
		int nid0 = face->GetNodeId(0);
		int nforce0 = nodeForceStore.Get(nid0).forces.GetNum();
		SRforce* newFaceForce = NULL;
		for (int iforce0 = 0; iforce0 < nforce0; iforce0++)
		{
			compatibleForceNum.Set(-1);
			compatibleForceNum.Put(0, iforce0);
			SRforce* nodeForce0 = nodeForceStore.Get(nid0).GetForce(iforce0);
			int coord0 = nodeForce0->coordId;
			bool pressure0 = nodeForce0->isPressure();
			double fv0[3];
			int ndof = 3;
			if (pressure0)
				ndof = 1;
			for (int dof = 0; dof < ndof; dof++)
				fv0[dof] = nodeForce0->forceVals.Get(0, dof);

			SRforce* nodeForce;

			for (int n = 1; n < nnodesTotal; n++)
			{
				int nid = face->GetNodeOrMidNodeId(n);
				int nforce = nodeForceStore.Get(nid).forces.GetNum();
				for (int iforce = 0; iforce < nforce; iforce++)
				{
					SRforce* nodeForce = nodeForceStore.Get(nid).GetForce(iforce);
					int coord = nodeForce->coordId;
					if ((coord == coord0) && (nodeForce->isPressure() == pressure0))
					{
						compatibleForceNum.Put(n, iforce);
						break;
					}
				}

			}
			bool compatible = true;
			int nn = nnodesTotal;
			for (int n = 1; n < nn; n++)
			{
				if (compatibleForceNum.Get(n) == -1)
				{
					compatible = false;
					break;
				}
			}
			if (compatible)
			{
				//can't have more than one force on same face:
				if (newFaceForce != NULL)
					ERROREXIT;
				numFaceForce++;
				if (countOnly)
					break;
				face->forceIds.PushBack(model.forces.GetNum());
				if (face->forceIds.GetNum() > 1)
					anyFaceHasMultipleLoads = true;
				newFaceForce = model.forces.Add();
				newFaceForce->type = faceForce;
				newFaceForce->faceFromNodal = true;
				numFaceFromNodalLoad++;
				newFaceForce->coordId = coord0;
				newFaceForce->pressure = pressure0;
				newFaceForce->entityId = f;
				newFaceForce->forceVals.Allocate(nnodesTotal, ndof);
				SRforce nodeForceTmp;
				nodeForceTmp.forceVals.Allocate(1, 3);
				for (int n = 0; n < nnodesTotal; n++)
				{
					int iforce = compatibleForceNum.Get(n);
					int nid;
					if (n == 0)
						nodeForce = nodeForce0;
					else
					{
						nid = face->GetNodeOrMidNodeId(n);
						nodeForce = nodeForceStore.Get(nid).GetForce(iforce);
					}
					for (int dof = 0; dof < ndof; dof++)
					{
						double ft = nodeForce->forceVals.Get(0, dof);
						newFaceForce->forceVals.Put(n, dof, ft);
					}
					//set this node force inactive to flag it should be removed below:
					nodeForce->type = inactiveForce;
				}
			} //endif: if (compatible)
		} //for (int iforce0 = 0; iforce0 < nforce0; iforce0++)
	} //for (int f = 0; f < nface; f++)

	if (countOnly)
		return numFaceForce;

	//for biaxial loading at corner of two shared faces, an adjacent face might pick up the loaded dof at the edge only
	//filter this out:
	for (int f = 0; f < model.GetNumFaces(); f++)
	{
		SRface* face = model.GetFace(f);
		if (face->forceIds.isEmpty())
			continue;
		SRforce* faceForce = model.GetForce(face->forceIds.Get(0));
		int nnodesTotal = face->GetNumNodesTotal();
		int ndof = 3;
		if (faceForce->pressure)
			ndof = 1;

		for (int dof = 0; dof < ndof; dof++)
		{
			int nej = face->localEdges.GetNum();
			bool edgeNode[8];
			for (int lej = 0; lej < nej; lej++)
			{
				for (int n = 0; n < nnodesTotal; n++)
					edgeNode[n] = false;
				bool entireEdgeLoaded = true;
				for (int ejnode = 0; ejnode < 3; ejnode++)
				{
					int ln = model.map.GetFaceEdgeLocalNode(nej, lej, ejnode);
					edgeNode[ln] = true;
					double ft = faceForce->forceVals.Get(ln, dof);
					if (fabs(ft) < TINY)
					{
						entireEdgeLoaded = false;
						break;
					}
				}
				if (!entireEdgeLoaded)
					continue;
				bool anyOtherNodeLoaded = false;
				for (int n = 0; n < nnodesTotal; n++)
				{
					if (edgeNode[n])
						continue;
					double ft = faceForce->forceVals.Get(n, dof);
					if (fabs(ft) > TINY)
					{
						anyOtherNodeLoaded = true;
						break;
					}
				}
				if (!anyOtherNodeLoaded)
				{
					//only one edge is loaded for this dof. filter out by setting all forcevals for this dof to 0:
					for (int n = 0; n < nnodesTotal; n++)
						faceForce->forceVals.Put(n, dof, 0.0);
				}
			}

			//same check for only one node loaded, would happen at corner where 3 faces meet
			int numLoadedNodes = 0;
			for (int n = 0; n < nnodesTotal; n++)
			{
				double ft = faceForce->forceVals.Get(n, dof);
				if (fabs(ft) > TINY)
					numLoadedNodes++;
			}
			if (numLoadedNodes == 1)
			{
				for (int n = 0; n < nnodesTotal; n++)
					faceForce->forceVals.Put(n, dof, 0.0);
			}
		} // end: for (int dof = 0; dof < 3; dof++)
	} // end: for (int f = 0; f < model.GetNumFaces(); f++)

	for (int f = 0; f < model.GetNumFaces(); f++)
	{
		SRface* face = model.GetFace(f);
		if (face->forceIds.isEmpty())
			continue;
		//empty forceIds, they will be reassigned after inactive nodal forces removed
		face->forceIds.Free();
	}

	return numFaceForce;
}

void SRinput::echoConstraints()
{
	OUTPRINT("constraints\n");
	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		SRconstraint* con = model.GetConstraint(i);
		int eid = con->GetEntityId();
		SRstring typestr;
		if (con->GetType() == nodalCon)
			typestr.Copy("Nodal");
		else if (con->GetType() == faceCon)
			typestr.Copy("Face");
		OUTPRINT("%d type, eid, coordid: %s %d %d\n", i, typestr.str, eid, con->GetCoordId());
		if (con->GetType() == faceCon)
		{
			SRface* face = model.GetFace(eid);
			OUTPRINTNORET("nodes: ");
			for (int j = 0; j < face->GetNumNodes(); j++)
				OUTPRINTNORET(" %d", face->GetNodeId(j));
			OUTPRINT("\n");
		}
		OUTPRINT(" constrained dofs: %d %d %d\n", con->IsConstrainedDof(0), con->IsConstrainedDof(1), con->IsConstrainedDof(2));
		OUTPRINT(" enf. disp: \n");
		for (int j = 0; j < 3; j++)
		{
			if (con->IsConstrainedDof(i))
			{
				for (int k = 0; k < con->GetNumEnforcedDisp(); k++)
					OUTPRINTNORET(" %lg", con->GetEnforcedDisp(k, j));
				OUTPRINT("\n");
			}
		}
	}
}

int SRinput::elemFaceFind(SRelement* elem, int nv[], int gno[])
{
	//find the face on an element that matches nodes in nv array
	//input:
		//elem = pointer to element
		//nv = corner nodes of faces. nv[3] = -1 for tri face
	//output:
		//gno = global node orders of the nodes on the face corresponding to nv
	//return
		//id of the global faces in the model with corner nodes nv; -1 if no match

	for (int l = 0; l < elem->localFaces.GetNum(); l++)
	{
		SRface* face = elem->GetFace(l);
		if (SRfaceUtil::FaceMatch(nv, face->nodeIds, gno))
			return elem->GetLocalFaceGlobalId(l);
	}
	return -1;
}

void SRinput::FillAndSortNodeUids()
{
	int nnode = model.GetNumNodes();
	if (nodeUidOffset == -1)
	{
		//fill node uid vector and sort in ascending order of uid for faster
		//node-finding:
		nodeUids.Allocate(nnode);
		SRuidData *nuid;
		for (int i = 0; i < nnode; i++)
		{
			SRnode* node = model.GetNode(i);
			nuid = nodeUids.GetPointer(i);
			nuid->id = i;
			nuid->uid = node->userId;
		}
		SortNodes();
	}
}


void SRinput::breakoutMeshPlot()
{

	// plot mesh to breakout frd
	int nnode = model.GetNumNodes();
	SRstring cgxSave;
	cgxSave.Copy(model.cgxFrdFile.filename);
	SRstring line;
	cgxSave.Left('.', line);
	line += "_bkt.frd";
	model.cgxFrdFile.filename = line;
	if (!model.cgxFrdFile.OutOpenNoFail())
	{
		model.cgxFrdFile.filename = cgxSave;
		return;
	}

	model.post.MeshToFrd();
	model.cgxFrdFile.Close();
	model.cgxFrdFile.filename = cgxSave;
}

void SRinput::MeshOnlyFrdPlot()
{

	//for debugging, plot mesh to frd file right after input
	int nnode = model.GetNumNodes();

	if (!model.cgxFrdFile.OutOpenNoFail())
	{
		OUTPRINT("unable to open cgx frd file for postprocessing");
		OUTPRINT("this model may be open already in Cgx postprocessor");
		REPPRINT("unable to open cgx frd file for postprocessing");
		REPPRINT("this model may be open already in Cgx postprocessor");
		ERROREXIT;
	}
	model.post.MeshToFrd();
	model.cgxFrdFile.Close();
	exit(0);
}

SRconstraint* SRnodeConInputStore::GetCon(int i)
{
	int id = cons.Get(i);
	return model.GetConstraint(id);
}

SRforce* SRnodeForceInputStore::GetForce(int i)
{
	int id = forces.Get(i);
	return model.GetForce(id);
}

void SRinput::GetBinaryFileName(char *ext, SRstring &name)
{
	//get the binary file name, full path, corresponding to a file name and extension
	//input:
	//ext = name of file extension
	//name = file name
	//note:
	//creates working directory ".tmp" if does not already exist
	name = model.wkdir;
	name += "tmp";
	SRfile::CreateDir(name.str);
	name += "\\";
	name += ext;
	name += ".bin";
}

void SRinput::equivMatTest()
{
	//check all materials for equivalence with other materials
	int nm = model.GetNumMaterials();
	for (int i = 0; i < nm; i++)
	{
		SRmaterial* mati = model.GetMaterial(i);
		for (int j = i + 1; j < nm; j++)
		{
			SRmaterial* matj = model.GetMaterial(j);
			if (mati->diffElast(matj))
			{
				model.allMatsEquivElast = false;
				return;
			}
		}

	}
}

int SRinput::GetBreakoutmatId(SRstring& saveBreakoutMat)
{
	//find the matid with name saveBreakoutMat and set model.saveBreakoutMatId

	int matid = -1;
	int nummat = model.GetNumMaterials();
	if (nummat == 1)
	{
		OUTPRINT("Error- saveBreakout by material can only be performed for models with multiple materials");
		LOGPRINT("Error- saveBreakout by material can only be performed for models with multiple materials");
		ERROREXIT;
	}
	for (int m = 0; m < nummat; m++)
	{
		SRmaterial* mat = model.GetMaterial(m);
		if (mat->name.CompareUseLength(saveBreakoutMat))
		{
			matid = m;
			break;
		}
	}
	if (matid == -1)
	{
		OUTPRINT("Error- saveBreakout by material. material with name %s not found", saveBreakoutMat.str);
		LOGPRINT("Error- saveBreakout by material. material with name %s not found", saveBreakoutMat.str);
		ERROREXIT;
	}

	return matid;
}

void SRinput::applyNodalBreakoutCons()
{
	if (numUnSupportedNode == 0)
		return;
	SRvec3 enfd;
	int numnodes = model.GetNumNodes();
	model.post.readPreviousDisplacementForBreakout();

	//nodes attached to unsupported elements need breakout constraints.
	//they were flagged on input as sacrificial:
	int ncon = model.GetNumConstraints();
	ncon += numUnSupportedNode;
	model.constraints.Allocate(ncon);
	SRconstraint* con;
	for (int i = 0; i < model.GetNumNodes(); i++)
	{
		SRnode* node = model.GetNode(i);
		if (node->checkUnsup() && !node->isOrphan())
		{
			node->constraintId = model.constraints.GetNum();
			con = model.constraints.Add();
			con->id = node->constraintId;
			con->type = nodalCon;
			//assign a coord id so the breakout constraint is forced to be processed as penalty constraint, to avoid conflicts:
			con->coordId = 0; //0 is reserved for default gcs coord system
			con->breakout = true;
			con->entityId = node->id;
			enfd.Copy(model.post.nodeDisps.Get(i));
			for (int dof = 0; dof < 3; dof++)
			{
				con->constrainedDof[dof] = true;
				con->PutEnforcedDisplacementData(0, dof, enfd.d[dof]);
			}
		}
	}

	//for any edges with breakout cons on both nodes, need one on midnode also if not already there:
	int nmidbreakout = 0;
	int ncon0 = model.GetNumConstraints();
	int nedge = model.GetNumEdges();
	model.constraints.Allocate(ncon0 + nedge);
	SRintVector nodeNeedsMidNodeCon;
	nodeNeedsMidNodeCon.Allocate(numnodes);
	for (int e = 0; e < nedge; e++)
	{
		SRedge* edge = model.GetEdge(e);
		enfd.Zero();
		bool bothnodesbreakout = true;
		SRnode* node;
		int nuids[2];
		for (int n = 0; n < 2; n++)
		{
			node = edge->GetNode(n);
			if (!node->checkUnsup())
			{
				bothnodesbreakout = false;
				break;
			}
			nuids[n] = node->userId;
		}
		if (bothnodesbreakout)
		{
			int mid = edge->GetMidNodeId();
			node = model.GetNode(mid);
			if (!node->checkUnsup())
			{
				nodeNeedsMidNodeCon.Put(mid, 1);
				node->sacrificialType = unsupported;
				node->constraintId = model.constraints.GetNum();
				SRconstraint* con = model.constraints.Add();
				con->id = node->constraintId;
				con->type = nodalCon;
				con->breakout = true;
				con->entityId = node->id;
				con->coordId = 0;//0 is reserved for default gcs coord system
			}
		}
	}

	SRvec3 enfdmid;
	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		SRconstraint* con = model.GetConstraint(i);
		if (con->type == nodalCon && con->breakout)
		{
			int nid = con->entityId;
			SRnode* node = model.GetNode(nid);
			if (nodeNeedsMidNodeCon.Get(nid) == 1)
			{
				con->allocateEnforcedDisplacementData(1);
				enfdmid = model.post.nodeDisps.Get(nid);
				for (int dof = 0; dof < 3; dof++)
				{
					con->constrainedDof[dof] = true;
					con->PutEnforcedDisplacementData(0, dof, enfdmid.d[dof]);
				}
			}
		}
	}
}

void SRinput::applyShellOrBeamBreakoutCons()
{
	//apply breakout constraints to faces with shell or beam nodes
	int nface = model.GetNumFaces();
	SRintVector faceWithShellNodes;
	faceWithShellNodes.Allocate(nface);
	int numFaceWithShellNodes = 0;

	for (int n = 0; n < model.GetNumNodes(); n++)
	{
		SRnode* node = model.GetNode(n);
		if (node->isOrphan())
			continue;
		if (node->sacrificialType != shellOrBeamNode)
			continue;
		for (int f = 0; f < numNodeFaces.Get(n); f++)
		{
			int fid = nodeFaces.Get(n, f);
			if (faceWithShellNodes.Get(fid) == 0)
			{
				numFaceWithShellNodes++;
				faceWithShellNodes.Put(fid, 1);
			}
		}
	}

	model.post.readPreviousDisplacementForBreakout();

	int ncon = model.GetNumConstraints();
	ncon += numFaceWithShellNodes;
	model.constraints.Allocate(ncon);
	SRconstraint* con;
	SRvec3 enfDisp;
	for (int f = 0; f < nface; f++)
	{
		if (faceWithShellNodes.Get(f) == 1)
		{
			int cid = model.constraints.GetNum();
			con = model.constraints.Add();
			con->id = cid;
			con->type = faceCon;
			//assign a coord id so the breakout constraint is forced to be processed as penalty constraint, to avoid conflicts:
			con->coordId = 0; //0 is reserved for default gcs coord system
			con->breakout = true;
			con->entityId = f;
			SRface* face = model.GetFace(f);
			model.anyUnsupportedFace = true;
			face->unsupported = true;
			face->constraintId = con->id;
			int nn = face->GetNumNodesTotal();
			con->enforcedDisplacementData.Allocate(nn, 3);
			for (int n = 0; n < nn; n++)
			{
				int nid = face->GetNodeOrMidNodeId(n);
				enfDisp.Copy(model.post.nodeDisps.Get(nid));
				for (int dof = 0; dof < 3; dof++)
					con->enforcedDisplacementData.Put(n, dof, enfDisp.d[dof]);
			}
		}
	}
}

void SRinput::InputnodalBreakoutConstraints(int nbr)
{
	int ncon = model.GetNumConstraints();
	ncon += nbr;
	model.constraints.Allocate(ncon);
	SRconstraint* con;
	//input nodal constraints

	//input:
	//lastSet = true if there is only one constraint set in the model or this is the last set

	//input spec: the following are all on one line:
	//node user id
	//enforced disp vals for each constrained dof, 0 if not enforced, "-" if not constrained
	//"coord" coord-name
	//notes:
	//edge and face constraints will be automatically created
	//if all nodes of the edge or face are constrained compatibly
	// after lastSet (in nodalToFaceConstraints)

	SRstring tok, line;
	int nodeuid;
	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (line == "end")
			break;
		line.TokRead(nodeuid);
		int nId = NodeFind(nodeuid);
		SRnode* node = model.GetNode(nId);
		node->constraintId = model.GetNumConstraints();
		node->sacrificialType = breakoutCon;
		SRconstraint* con = model.constraints.Add();
		con->type = nodalCon;
		con->coordId = 0; //0 is reserved for default gcs coord system
		con->breakout = true;
		con->entityId = node->id;
		double enfd;
		for (int dof = 0; dof < 3; dof++)
		{
			con->constrainedDof[dof] = true;
			line.TokRead(enfd);
			if (fabs(enfd) > TINY)
				model.anyEnforcedDisplacement = true;
			if (con->enforcedDisplacementData.isEmpty())
				con->enforcedDisplacementData.Allocate(1, 3);
			con->enforcedDisplacementData.Put(0, dof, enfd);
		}
	}

}

int SRinput::CountMultiFaceConOrForce(int numgroups)
{
	int numnodetotal = 0;
	SRstring line;
	int numnodeingroup;
	for (int i = 0; i < numgroups; i++)
	{
		model.inputFile.GetLine(line); //header
		line.Token();
		line.TokRead(numnodeingroup);
		numnodetotal += numnodeingroup;
		for (int n = 0; n < numnodeingroup; n++)
			model.inputFile.GetLine(line); //skip data line
	}
	model.inputFile.GetLine(line); //"end multi"
	return numnodetotal;
}

void SRinput::fillMultiFaceForceGroups()
{
	if (numFaceFromNodalLoad == 0)
		return;

	if (anyFaceHasMultipleLoads)
		return;

	nodeNeedsNodalForce.Allocate(model.GetNumNodes());
	nodeNeedsNodalForce.Set(-1);

	int nface = model.GetNumFaces();
	model.faceForceGroups.Allocate(nface);
	for (int f = 0; f < nface; f++)
	{
		SRface* face = model.GetFace(f);
		if (!face->hasForce())
			continue;
		if (face->multifaceForceGroupId != -1)
			continue;
		//create new surface for this face:
		int ffid = model.faceForceGroups.GetNum();
		SRFaceForceGroup* ff = model.faceForceGroups.Add();
		face->multifaceForceGroupId = ffid;
		ff->addFace(face);
		checkBoundaryEdgeNeighbors(ffid, face);
		int fid = face->forceIds.Get(0);
		SRforce* force = model.GetForce(fid);
		if (force->faceFromNodal)
			ff->fromNodalForces = true;
		force->faceFromNodal = false;
	}

	SRintVector nodeonff;
	int nnode = model.GetNumNodes();
	nodeonff.Allocate(nnode);
	for (int ffn = 0; ffn < model.faceForceGroups.GetNum(); ffn++)
	{
		nodeonff.Zero();
		SRFaceForceGroup* ff = model.faceForceGroups.GetPointer(ffn);
		for (int f = 0; f < ff->faceIds.GetNum(); f++)
		{
			SRface* face = model.GetFace(ff->faceIds.Get(f));
			for (int n = 0; n < face->GetNumNodesTotal(); n++)
			{
				int nid = face->GetNodeOrMidNodeId(n);
				nodeonff.Put(nid, 1);
			}
		}
		for (int n = 0; n < nnode; n++)
		{
			if (nodeonff.Get(n) == 1)
			{
				ff->nodeIds.PushBack(n);
				nodeNeedsNodalForce.Put(n, 1);
			}
		}
		ff->nodeIds.Sort();
		/*
		OUTPRINT("ffg: %d nodes:", ffn);
		for (int n = 0; n < ff->nodeIds.GetNum(); n++)
		{
			int nid = ff->nodeIds.Get(n);
			int nuid = model.GetNode(nid)->userId;
			OUTPRINT("%d ", nuid);
		}
		*/
	}

	//LOGPRINT("number of faceForceGroups: %d", model.faceForceGroups.GetNum());
	//plot surfaces
	bool dbgplotFaceForceGroups = false;
	if (!dbgplotFaceForceGroups)
		return;
	for (int s = 0; s < model.faceForceGroups.GetNum(); s++)
		PlotFaceForceGroups(s);
}

void SRinput::checkBoundaryEdgeNeighbors(int parentMultifaceForceGroupId, SRface* parentFace)
{
	//loop over this face's local edges. see if face neighbors are on same surface as parentMultifaceForceGroupId (using normal kink test)

	for (int lej = 0; lej < parentFace->GetNumLocalEdges(); lej++)
	{
		SRedge* edge = parentFace->GetEdge(lej);
		for (int ef = 0; ef < edge->boundaryfaceData.GetNum(); ef++)
		{
			int fid2 = edge->boundaryfaceData.Get(ef).faceId;
			if (fid2 == parentFace->id)
				continue;
			int lej2 = edge->boundaryfaceData.Get(ef).localEdgeId;
			SRface* face2 = model.GetFace(fid2);
			if (!face2->hasForce())
				continue;
			if (!edge->checkKinkOK(parentFace, lej, face2, lej2))
				continue;
			if (face2->multifaceForceGroupId != -1)
				continue;
			face2->multifaceForceGroupId = parentMultifaceForceGroupId;
			SRFaceForceGroup* ff = model.faceForceGroups.GetPointer(parentMultifaceForceGroupId);
			ff->addFace(face2);
			checkBoundaryEdgeNeighbors(parentMultifaceForceGroupId, face2);
		}
	}
}



void SRinput::PlotFaceForceGroups(int ffnum)
{
	SRstring name;
	name.Copy(model.outdir);
	name.Cat("\\faceForceGroup");
	char buf[20];
	sprintf_s(buf, "%d", ffnum);
	name.Cat(buf);
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
		frdf.PrintLine("\n");
	}
	frdf.PrintLine(" -3");
	//displacements

	SRvector <SRvec3> nodeDispsTmp;
	SRvec3 allOnes;
	allOnes.Assign(1.0, 1.0, 1.0);
	nodeDispsTmp.Allocate(model.GetNumNodes());
	SRFaceForceGroup* ff = model.faceForceGroups.GetPointer(ffnum);

	for (int n = 0; n < ff->nodeIds.GetNum(); n++)
	{
		int nid = ff->nodeIds.Get(n);
		nodeDispsTmp.Get(nid) = allOnes;
	}

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
		disp = nodeDispsTmp.Get(i);
		frdf.PrintLine(" -1%10d%12.5le%12.5le%12.5le", node->userId, disp.d[0], disp.d[1], disp.d[2]);
	}
	frdf.PrintLine(" -3");

	frdf.PrintLine("9999"); //end of data

	frdf.Close();
}

void SRinput::EnergySmoothNodalForces()
{
	//convert to smooth tractions, energy approach

	if (!model.doEnergySmooth)
		return;

	int nfg = model.faceForceGroups.GetNum();

	bool useSymFull = false;
	if (useSymFull)
	{
		EnergySmoothNodalForcesSymFull();
		return;
	}


	for (int g = 0; g < nfg; g++)
	{
		SRFaceForceGroup* ffg = model.faceForceGroups.GetPointer(g);
		if (!ffg->fromNodalForces)
			continue;
		SRpardiso parSolver;
		ffg->faceFunLoc.Allocate(8,ffg->faceIds.GetNum());
		double N[8];
		double rf, sf, w;
		int nfun = ffg->nodeIds.GetNum();
		ffg->ForceDof.Allocate(3, ffg->nodeIds.GetNum());
		for (int f = 0; f < ffg->faceIds.GetNum(); f++)
		{
			int fid = ffg->faceIds.Get(f);
			SRface* face = model.GetFace(fid);
			int forceid = face->forceIds.Get(0);
			SRforce* faf = model.GetForce(forceid);
			int nfunface = face->GetNumNodesTotal();
			int fafmatlen = nfunface*(nfunface + 1) / 2;
			faf->StiffMat.Allocate(fafmatlen);
			faf->stiffDiag.Allocate(nfunface);
			int nz = 0;
			for (int r = 0; r < nfunface; r++)
			{
				faf->stiffDiag.Put(r, nz);
				nz += (nfunface - r);
			}
			for (int i = 0; i < nfunface; i++)
			{
				int nid = face->GetNodeOrMidNodeId(i);
				int loc = ffg->nodeIds.Find(nid);
				if (loc == -1)
					ERROREXIT;
				ffg->faceFunLoc.Put(i, f, loc);
			}
			int nint = model.math.FillGaussPoints(face);
			for (int gp = 0; gp < nint; gp++)
			{
				model.math.GetGP2d(gp, rf, sf, w);
				w *= (face->Jacobian(rf, sf));
				model.map.FaceShapeFunctions(face, rf, sf, N);
				for (int i = 0; i < nfunface; i++)
				{
					double niw = N[i] * w;
					int iloc = ffg->faceFunLoc.Get(i, f);
					for (int dof = 0; dof < 3; dof++)
						ffg->ForceDof.Put(dof, iloc, faf->forceVals.Get(i, dof));
					for (int j = i; j < nfunface; j++)
					{
						double NN = N[j] * niw;
						int ij = parSolver.getFFGElemStiffLocation(faf, i, j);
						faf->StiffMat.PlusAssign(ij, NN);
					}
				}
			}
		}
		parSolver.solveFFG(ffg);
		SRdoubleVector Fv;
		Fv.Allocate(nfun);
		double *F = Fv.d;

		for (int dof = 0; dof < 3; dof++)
		{
			for (int j = 0; j < nfun; j++)
				F[j] = ffg->ForceDof.Get(dof, j);
			//F was overwritten with tractions:
			for (int f = 0; f < ffg->faceIds.GetNum(); f++)
			{
				int fid = ffg->faceIds.Get(f);
				SRface* face = model.GetFace(fid);
				int forceid = face->forceIds.Get(0);
				SRforce* faf = model.GetForce(forceid);
				faf->coordId = -1;
				int nfunface = face->GetNumNodesTotal();
				for (int i = 0; i < nfunface; i++)
				{
					int iloc = ffg->faceFunLoc.Get(i, f);
					faf->forceVals.Put(i, dof, F[iloc]);
				}
				faf->stiffDiag.Free();
				faf->StiffMat.Free();
			}
		}
		ffg->faceFunLoc.Free();
		ffg->smoothStiff.Free();
	}

	bool echoSmoothTractions = false;
	if (!echoSmoothTractions)
		return;
	OUTPRINT("\nsmoothed face tractions");
	OUTPRINT(" node uid   tx      ty     tz");
	for (int g = 0; g < nfg; g++)
	{
		SRFaceForceGroup* ffg = model.faceForceGroups.GetPointer(g);
		if (!ffg->fromNodalForces)
			continue;
		for (int f = 0; f < ffg->faceIds.GetNum(); f++)
		{
			int fid = ffg->faceIds.d[f];
			SRface* face = model.GetFace(fid);
			int forceid = face->forceIds.Get(0);
			SRforce* faf = model.GetForce(forceid);
			OUTPRINT("local face: %d", f);
			for (int n = 0; n < faf->forceVals.getNumCols(); n++)
			{
				int nuid = face->GetNodeOrMidnode(n)->userId;
				OUTPRINTNORET("%d", nuid);
				for (int dof = 0; dof < 3; dof++)
				{
					double ft = faf->forceVals.Get(n, dof);
					OUTPRINTNORET(" %lg", ft);
				}
				OUTPRINTRET;
			}
		}
	}
}

void SRinput::PlotUnsupportedFaces()
{
	SRstring name;
	name.Copy(model.outdir);
	name.Cat("\\unsupFaces");
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
		frdf.PrintLine("\n");
	}
	frdf.PrintLine(" -3");
	//displacements

	SRintVector nodesOnUnsupFaces;
	nodesOnUnsupFaces.Allocate(numnode);
	for (int f = 0; f < model.faces.GetNum(); f++)
	{
		SRface* face = model.GetFace(f);
		if (face->unsupported)
		{
			for (int i = 0; i < face->GetNumNodesTotal(); i++)
			{
				int nid = face->GetNodeOrMidNodeId(i);
				nodesOnUnsupFaces.Put(nid, 1);
			}
		}
	}

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
		if (nodesOnUnsupFaces.Get(i) == 1)
			disp.Assign(1.0, 1.0, 1.0);
		frdf.PrintLine(" -1%10d%12.5le%12.5le%12.5le", node->userId, disp.d[0], disp.d[1], disp.d[2]);
	}
	frdf.PrintLine(" -3");

	frdf.PrintLine("9999"); //end of data

	frdf.Close();
}

void SRinput::EnergySmoothNodalForcesSymFull()
{
	//convert to smooth tractions, energy approach:

	int nfg = model.faceForceGroups.GetNum();

	for (int g = 0; g < nfg; g++)
	{
		SRFaceForceGroup* ffg = model.faceForceGroups.GetPointer(g);
		if (!ffg->fromNodalForces)
			continue;
		SRintMatrix faceFunLoc;
		faceFunLoc.Allocate(8, ffg->faceIds.GetNum());
		double N[8];
		double rf, sf, w;
		int nfun = ffg->nodeIds.GetNum();
		SRskyline NNsky(nfun);
		NNsky.FillSymDiag();
		SRdoubleMatrix rhs;
		SRdoubleVector F;
		rhs.Allocate(3, nfun);
		F.Allocate(nfun);
		for (int f = 0; f < ffg->faceIds.GetNum(); f++)
		{
			int fid = ffg->faceIds.Get(f);
			SRface* face = model.GetFace(fid);
			int forceid = face->forceIds.Get(0);
			SRforce* faf = model.GetForce(forceid);
			int nfunface = face->GetNumNodesTotal();
			for (int i = 0; i < nfunface; i++)
			{
				int nid = face->GetNodeOrMidNodeId(i);
				int loc = ffg->nodeIds.Find(nid);
				if (loc == -1)
					ERROREXIT;
				faceFunLoc.Put(i, f, loc);
			}
			int nint = model.math.FillGaussPoints(face);
			for (int gp = 0; gp < nint; gp++)
			{
				model.math.GetGP2d(gp, rf, sf, w);
				w *= (face->Jacobian(rf, sf));
				model.map.FaceShapeFunctions(face, rf, sf, N);
				for (int i = 0; i < nfunface; i++)
				{
					double niw = N[i] * w;
					int iloc = faceFunLoc.Get(i, f);
					for (int dof = 0; dof < 3; dof++)
						rhs.Put(dof, iloc, faf->forceVals.Get(i, dof));
					for (int j = 0; j < nfunface; j++)
					{
						double NN = N[j] * niw;
						int jloc = faceFunLoc.Get(j, f);
						if (jloc < iloc)
							continue;
						int ij = NNsky.Location(iloc, jloc);
						NNsky.stiffMat.PlusAssign(ij, NN);
					}
				}
			}
		}
		if (!NNsky.Decomp())
			ERROREXIT;
		for (int dof = 0; dof < 3; dof++)
		{
			for (int j = 0; j < nfun; j++)
				F.d[j] = rhs.Get(dof, j);
			if (!NNsky.BackSolve(F.d))
				ERROREXIT;
			//F was overwritten with tractions:
			for (int f = 0; f < ffg->faceIds.GetNum(); f++)
			{
				int fid = ffg->faceIds.Get(f);
				SRface* face = model.GetFace(fid);
				int forceid = face->forceIds.Get(0);
				SRforce* faf = model.GetForce(forceid);
				int nfunface = face->GetNumNodesTotal();
				for (int i = 0; i < nfunface; i++)
				{
					int iloc = faceFunLoc.Get(i, f);
					faf->forceVals.Put(i, dof, F[iloc]);
				}
			}
		}
	}

	bool echoSmoothTractions = false;
	if (!echoSmoothTractions)
		return;
	OUTPRINT("\nsmoothed face tractions");
	OUTPRINT(" node uid   tx      ty     tz");
	for (int g = 0; g < nfg; g++)
	{
		SRFaceForceGroup* ffg = model.faceForceGroups.GetPointer(g);
		if (!ffg->fromNodalForces)
			continue;
		for (int f = 0; f < ffg->faceIds.GetNum(); f++)
		{
			int fid = ffg->faceIds.d[f];
			SRface* face = model.GetFace(fid);
			int forceid = face->forceIds.Get(0);
			SRforce* faf = model.GetForce(forceid);
			OUTPRINT("local face: %d", f);
			for (int n = 0; n < faf->forceVals.getNumCols(); n++)
			{
				int nuid = face->GetNodeOrMidnode(n)->userId;
				OUTPRINTNORET("%d", nuid);
				for (int dof = 0; dof < 3; dof++)
				{
					double ft = faf->forceVals.Get(n, dof);
					OUTPRINTNORET(" %lg", ft);
				}
				OUTPRINTRET;
			}
		}
	}
}

void SRinput::nodalToFaceUnsup()
{
	if (numUnSupportedNode == 0)
		return;
	for (int f = 0; f < model.faces.GetNum(); f++)
	{
		bool isBsurf = false;
		SRface* face = model.GetFace(f);
		bool allCornersUnsup = true;
		int nc = face->GetNumNodes();
		for (int i = 0; i < nc; i++)
		{
			SRnode* node = face->GetNode(i);
			if (!node->checkUnsup())
			{
				allCornersUnsup = false;
				break;
			}
			if (node->sacrificialType == onBsurf)
				isBsurf = false;
		}
		if (!allCornersUnsup)
			continue;
		model.anyUnsupportedFace = true;
		face->unsupported = true;
		for (int e = 0; e < nc; e++)
		{
			SRnode* node = face->GetNodeOrMidnode(e + nc);
			if (!node->checkUnsup())
				numUnSupportedNode++;
			if (isBsurf)
				node->sacrificialType = onBsurf;
			else
				node->sacrificialType = unsupported;
		}
	}

	//PlotUnsupportedFaces();
}

void SRinput::PlotBktElSacr()
{
	SRintVector ellist;
	for (int e = 0; e < model.post.bktElSacrList.GetNum(); e++)
	{
		int euid = model.post.bktElSacrList.Get(e);
		int eid = elemFind(euid);
		if (eid != -1)
			ellist.PushBack(eid);
	}
	model.post.PlotElems(ellist.GetNum(), ellist.d, "inputPlotBktElSacr");
}



