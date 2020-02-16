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
// SRelement.cpp: implementation of the SRelement class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SRmodel.h"

extern SRmodel model;

SRBTC::SRBTC()
{
	BTC11 = BTC12 = BTC13 = BTC14 = BTC15 = BTC16 = 0.0;
	BTC21 = BTC22 = BTC24 = BTC24 = BTC25 =  BTC26 = 0.0;
	BTC21 = BTC32 = BTC33 = BTC34 = BTC35 = BTC36 = 0.0;
}


SRface * SRlocalFace::GetFace()
{
	//get pointer to global face corresponding to this local face of an element
	//return:
		//pointer to global face

	return model.GetFace(globalFaceId);
}


SRelement::SRelement()
{
	type = tet;
	sacrificial = 0;
	svmMax = 0.0;
	flattened = false;
	flattenedHighStress = false;
	stiffLength = 0;
	constrained = false;
	elVol = 0.0;
	saveForBreakout = false;
	flattenFraction = 0.0;
	approxVol = 0.0;
	allNodesHaveDisp = false;
};

void SRelement::Create(int useridt, int nnodes, int nodest[], SRmaterial* mat)
{
	if (nnodes == 10)
	{
		Create(tet, useridt, nnodes, nodest, mat);
	}
	else if (nnodes == 15)
	{
		Create(wedge, useridt, nnodes, nodest, mat);
	}
	else if (nnodes == 20)
	{
		Create(brick, useridt, nnodes, nodest, mat);
	}

}
void SRelement::Create(SRelementType typet, int useridt, int nnodes, int nodest[], SRmaterial* mat)
{
	//Create element. store nodes, material, edges, faces in element arrays
	//input:
		//typet = element type
		//useridt = user ID
		//nodest = vector of node ids
		//mat = material property

	type = typet;
	userId = useridt;
	elMat = mat;
	int ncorner, nej;
	penaltyStiffness = 0.0;
	if (type == tet)
	{
		ncorner = 4;
		nej = 6;
		localEdges.Allocate(nej);
		localFaces.Allocate(4);

	}
	else if (type == wedge)
	{
		ncorner = 6;
		nej = 9;
		localEdges.Allocate(nej);
		localFaces.Allocate(5);
	}
	else
	{
		ncorner = 8;
		nej = 12;
		localEdges.Allocate(nej);
		localFaces.Allocate(6);
	}

	if (!model.saveBreakout)
	{
		nodeIds.Allocate(ncorner);
		for (int i = 0; i < ncorner; i++)
			nodeIds.Put(i, nodest[i]);
	}
	else
	{
		int nnodesTotal = ncorner + nej;
		nodeIds.Allocate(nnodesTotal);
		for (int i = 0; i < nnodesTotal; i++)
			nodeIds.Put(i, nodest[i]);
	}
	pChanged = true;
	error = 0.0;
}

bool SRelement::DetectThinEdges()
{
	//detect if any of this elements edges are thin.
	//notes:
		//criteria for thin edges:
			//1. the length of all edges parallel to a direction is less than 1/10 of the max-inplane direction,
				//e.g. if direction is r, the inplane directions are s,t
			//2. all the edges parallel to the thin direction must be straight
			//3. all the edges parallel to the thin direction must be normal to the face that touches the edge's 1st node

	SRedge* edge;
	SRvec3 edgeTang, faceNorm;
	double thinTol = 0.101;
	if (type == tet)
		return false;
	else if (type == brick)
	{
		int rej[4] = { 0, 2, 4, 6 };
		int sej[4] = { 1, 3, 5, 7 };
		double rlmin = BIG, rlmax = 0.0, slmin = BIG, slmax = 0.0, tlmin = BIG, tlmax = 0.0;
		bool rejStraight = true;
		bool sejStraight = true;
		bool tejStraight = true;
		for (int i = 0; i < 4; i++)
		{
			edge = GetEdge(rej[i]);
			if (!edge->isStraight())
				rejStraight = false;
			double len = edge->GetSize();
			if (len > rlmax)
				rlmax = len;
			if (len < rlmin)
				rlmin = len;
			edge = GetEdge(sej[i]);
			if (!edge->isStraight())
				sejStraight = false;
			len = edge->GetSize();
			if (len > slmax)
				slmax = len;
			if (len < slmin)
				slmin = len;
			edge = GetEdge(i+8);
			if (!edge->isStraight())
				tejStraight = false;
			len = edge->GetSize();
			if (len > tlmax)
				tlmax = len;
			if (len < tlmin)
				tlmin = len;
		}
		double inplanest = MATHMIN(slmin, tlmin);
		double inplanert = MATHMIN(rlmin, tlmin);
		double inplaners = MATHMIN(rlmin, slmin);
		//r-edges:
		if (rlmax < thinTol*inplanest && rejStraight)
		{
			//thin direction is r. make sure r edges are normal to their faces
			int edgeCorner0[4] = { 0, 3, 4, 7 };
			//element local face corresponding to node 0 of the r edges is 0-3-7-4 = 2
			SRface* face = GetFace(2);
			bool allNormal = true;
			double rf, sf;
			for (int i = 0; i < 4; i++)
			{
				edge = GetEdge(rej[i]);
				model.map.EdgeTangent(edge, 0.0, edgeTang, true);
				int elnodeAtCorner0 = GetNodeId(edgeCorner0[i]);
				if (!face->natCoordsAtCorner(elnodeAtCorner0, rf, sf))
				{
					allNormal = false;
					break;
				}
				face->OutwardNormal(rf, sf, faceNorm, false);
				//face normal should be parallel to edgeTang. ok if opposite sense
				double dot = fabs(faceNorm.Dot(edgeTang));
				if ((1.0 - dot) > SMALL)
				{
					allNormal = false;
					break;
				}
			}
			if (allNormal)
			{
				//ttd only warn once per model
				OUTPRINT("thin element detected: %d, p-order will be set to %d in thin direction", GetUserid(), model.GetThinEdgeMaxPorder());
				for (int i = 0; i < 4; i++)
					GetEdge(rej[i])->SetThin(true);
				return true;
			}
		}
		//s-edges:
		if (slmax < thinTol*inplanert && sejStraight)
		{
			//thin direction is s. make sure s edges are normal to their faces
			int edgeCorner0[4] = { 1, 0, 5, 4 };
			//element local face corresponding to node 0 of the r edges is 0-1-5-4 = 4
			SRface* face = GetFace(4);
			bool allNormal = true;
			double rf, sf;
			for (int i = 0; i < 4; i++)
			{
				edge = GetEdge(rej[i]);
				model.map.EdgeTangent(edge, 0.0, edgeTang, true);
				int elnodeAtCorner0 = GetNodeId(edgeCorner0[i]);
				if (!face->natCoordsAtCorner(elnodeAtCorner0, rf, sf))
				{
					allNormal = false;
					break;
				}
				face->OutwardNormal(rf, sf, faceNorm, false);
				//face normal should be parallel to edgeTang. ok if opposite sense
				double dot = fabs(faceNorm.Dot(edgeTang));
				if ((1.0 - dot) > SMALL)
				{
					allNormal = false;
					break;
				}
			}
			if (allNormal)
			{
				OUTPRINT("thin element detected: %d, p-order will be set to %d in thin direction", GetUserid(), model.GetThinEdgeMaxPorder());
				for (int i = 0; i < 4; i++)
					GetEdge(sej[i])->SetThin(true);
				return true;
			}
		}
		//t-edges:
		if (tlmax < thinTol*inplaners && tejStraight)
		{
			//thin direction is t. make sure t edges are normal to their faces
			int edgeCorner0[4] = { 0, 1, 2, 3 };
			//element local face corresponding to node 0 of the t edges is 0-1-2-3 = 0
			SRface* face = GetFace(0);
			bool allNormal = true;
			double rf, sf;
			for (int i = 0; i < 4; i++)
			{
				//t-edges are 8,9,10,11
				edge = GetEdge(i+8);
				model.map.EdgeTangent(edge, 0.0, edgeTang, true);
				int elnodeAtCorner0 = GetNodeId(edgeCorner0[i]);
				if (!face->natCoordsAtCorner(elnodeAtCorner0, rf, sf))
				{
					allNormal = false;
					break;
				}
				face->OutwardNormal(rf, sf, faceNorm, false);
				//face normal should be parallel to edgeTang. ok if opposite sense
				double dot = fabs(faceNorm.Dot(edgeTang));
				if ((1.0 - dot) > SMALL)
				{
					allNormal = false;
					break;
				}
			}
			if (allNormal)
			{
				OUTPRINT("thin element detected: %d, p-order will be set to %d in thin direction", GetUserid(), model.GetThinEdgeMaxPorder());
				for (int i = 0; i < 4; i++)
					GetEdge(i+8)->SetThin(true);
				return true;
			}
		}
	}
	else if (type == wedge)
	{
		//only t-edges may be thin
		//in-plane edges are 0-5
		//t-edges are > 5
		bool tejStraight = true;

		double inplaneMin = BIG;
		double tMax = 0.0;
		for (int i = 0; i < 6; i++)
		{
			double len = GetEdge(i)->GetSize();
			if (len < inplaneMin)
				inplaneMin = len;
		}
		for (int i = 6; i < 9; i++)
		{
			SRedge* edge = GetEdge(i);
			if (!edge->isStraight())
				tejStraight = false;
			double len = edge->GetSize();
			if (len > tMax)
				tMax = len;
		}
		bool thin = false;
		if (tMax < thinTol*inplaneMin && tejStraight)
		{
			//thin direction is t. make sure t edges are normal to their faces
			int edgeCorner0[3] = { 0, 1, 2 };
			//element local face corresponding to node 0 of the r edges is 0-1-2 = 0
			SRface* face = GetFace(0);
			bool allNormal = true;
			double rf, sf;
			for (int i = 6; i < 9; i++)
			{
				edge = GetEdge(i);
				model.map.EdgeTangent(edge, 0.0, edgeTang, true);
				int elnodeAtCorner0 = GetNodeId(i - 6);
				if (!face->natCoordsAtCorner(elnodeAtCorner0, rf, sf))
				{
					allNormal = false;
					break;
				}
				face->OutwardNormal(rf, sf, faceNorm, false);
				//face normal should be parallel to edgeTang. ok if opposite sense
				double dot = fabs(faceNorm.Dot(edgeTang));
				if ((1.0 - dot) > SMALL)
				{
					allNormal = false;
					break;
				}
			}
			if (allNormal)
			{
				OUTPRINT("thin element detected: %d, p-order will be set to %d in thin direction", GetUserid(), model.GetThinEdgeMaxPorder());
				for (int i = 6; i < 9; i++)
					GetEdge(i)->SetThin(true);
				return true;
			}
		}
	}
	return false;
}

void SRelement::PutRawStrain(int i, double *strain)
{
	//store raw strains at a gauss point of this elements
	//input:
		//i = gauss point number
		//strain = strain tensor
	//note:
		//stores strain in element rawStrains matrix

	for (int e = 0; e < 6; e++)
		rawStrains.Put(i, e, strain[e]);
}


void SRelement::GetEdgeLocalNodeIds(int lej, int& ln0, int& ln1)
{
	//find the element local node ids of the two ends of a local edge
	//input:
		//lej = local edge number
	//output:
		//ln0, ln1 = element local node ids of the ends of the local edge
	if (type == tet)
	{
		ln0 = model.map.GetTetEdgeLocalNode(lej, 0);
		ln1 = model.map.GetTetEdgeLocalNode(lej, 1);
	}
	else if (type == wedge)
	{
		ln0 = model.map.GetWedgeEdgeLocalNode(lej, 0);
		ln1 = model.map.GetWedgeEdgeLocalNode(lej, 1);
	}
	else
	{
		ln0 = model.map.GetBrickEdgeLocalNode(lej, 0);
		ln1 = model.map.GetBrickEdgeLocalNode(lej, 1);
	}
}


void SRelement::GetEdgeNodeIds(int lej, int &n1, int &n2)
{
	//get the node ids at the ends of a local edge of an element.
	//this is for use before global edge assignment, so this cannot be looked
	//up with localEdges.
	//input:
		//lej = local edge number
	//output:
		//n1,n2 = nodes at ends of edge
	int ln0, ln1;
	GetEdgeLocalNodeIds(lej, ln0, ln1);
	n1 = nodeIds.Get(ln0);
	n2 = nodeIds.Get(ln1);
}

int SRelement::GetLocalEdgeNodeId(int lej, int lnode)
{ 
	//look up the global node id corresponding to a local node of a local edge

	return localEdges.Get(lej).GetEdge()->GetNodeId(lnode);
}

void SRelement::GetFaceNodes(int lface, int &n1, int &n2, int &n3, int &n4)
{
	//get the node numbers at the corners of a local face of an element.
	//this is for use before global face assignment, so this cannot be looked
	//up with localFaces.
	//input:
		//lface = local face number
	//output:
		//n1,n2,n3,n4 = nodes at corners of face; n4 = -1 for tri face

	if (type == tet)
	{
		if (lface == 0)
		{
			//local face 1 = 2-3-4 (face for which L1=0)
			n1 = nodeIds.Get(1);
			n2 = nodeIds.Get(2);
			n3 = nodeIds.Get(3);
		}
		else if (lface == 1)
		{
			//local face 2 = 1-3-4 (face for which L2=0)
			n1 = nodeIds.Get(0);
			n2 = nodeIds.Get(2);
			n3 = nodeIds.Get(3);
		}
		else if (lface == 2)
		{
			//local face 3 = 1-2-4 (face for which L3=0)
			n1 = nodeIds.Get(0);
			n2 = nodeIds.Get(1);
			n3 = nodeIds.Get(3);
		}
		else if (lface == 3)
		{
			//local face 4 = 1-2-3 (face for which L4=0)
			n1 = nodeIds.Get(0);
			n2 = nodeIds.Get(1);
			n3 = nodeIds.Get(2);
		}
		n4 = -1;
	}
	else if (type == wedge)
	{
		GetWedgeFaceNodes(lface, n1, n2, n3, n4);
	}
	else
	{
		GetBrickFaceNodes(lface, n1, n2, n3, n4);
	}
}

SRedge* SRelement::GetEdge(int localedgenum)
{
	//look up global edge corresponding to element local edge
	//input:
		//localedgenum = element local edge number
	//return
		//pointer to global edge

	int gej = localEdges.GetPointer(localedgenum)->GetGlobalEdgeId();
	return model.GetEdge(gej);
}

SRface* SRelement::GetFace(int localfacenum)
{
	//look up global face corresponding to element local face
	//input:
		//localfacenum = element local face number
	//return
		//pointer to global face

	return localFaces.GetPointer(localfacenum)->GetFace();
}


SRnode* SRelement::GetNode(int localnodenum)
{
	//get the pointer to the global node corresponding to local node number
	//input:
		//localnodenum = local number in element
	//return:
		//pointer to the global node

	int nodeid = nodeIds.Get(localnodenum);
	return model.GetNode(nodeid);
}

double* SRelement::CalculateStiffnessMatrix(int processorNum, int& len)
{
	//Calculate Element Stiffness Matrix for a three dimensional element
	//input:
		//processorNum =  else processor calling this element, 0 for single-thread
	//output:
		//length of stiffness matrix
	//return:
		//Element Stiffness Matrix stored symmetrically
	thread = processorNum;

	double *stiff = NULL;

	int i, neq;
	int nfun = globalFunctionNumbers.GetNum();
	int nint, gp, coleq0;
	SRBTC btc;
	SRmap* map = &model.map;
	double w, dbdx, dbdy, dbdz;
	nint = model.math.FillGaussPoints(this);

	neq = 3 * nfun;
	FillStiffDiag(neq);
	len = neq*(neq + 1) / 2;

	if (model.areElementsInMemory())
	{
		if (len > stiffnessMatrix.GetNum())
			stiffnessMatrix.Allocate(len);
		//don't need to recalculate stiffness if p not changed:
		stiff = stiffnessMatrix.GetVector();
		if (!pChanged)
			return stiff;
	}
	else
	{
		if (!pChanged)
		{
			//try to read in previous p:
			stiff = model.ReadElementStiffness(this, true, processorNum);
			if (stiff != NULL)
				return stiff;
		}
		else
			stiff = model.GetElementStiffnessVector(processorNum);
	}
	for (i = 0; i < len; i++)
		stiff[i] = 0.0;

	SRElementData *eldata = model.GetElementData(processorNum);

	SRdoubleMatrix& dbdxm = eldata->Getdbdx();
	SRdoubleMatrix& dbdym = eldata->Getdbdy();
	SRdoubleMatrix& dbdzm = eldata->Getdbdz();
	double *wv = eldata->GetIntwt();
	basisVec = eldata->GetBasisVec();
	dbasisdr = eldata->Getdbasisdr();
	dbasisds = eldata->Getdbasisds();
	dbasisdt = eldata->Getdbasisdt();

	elVol = 0.0;
	for (gp = 0; gp < nint; gp++)
	{
		double r, s, t;
		model.math.GetGP3d(gp, r, s, t, w, thread);
		double detj = FillMapping(r, s, t);
		if (detj < approxVol*SMALL)
		{

			OUTPRINT(" bad Jacobian. Element: %d", GetUserid());
			ERROREXIT;
		}
		w *= detj;
		elVol += w;
		wv[gp] = w;
		FillBasisFuncs(r, s, t, derivonly);

		for (int fun = 0; fun < nfun; fun++)
		{
			double dbdr = dbasisdr[fun];
			double dbds = dbasisds[fun];
			double dbdt = dbasisdt[fun];
			XyzDerivatives(dbdr, dbds, dbdt, dbdx, dbdy, dbdz);
			dbdxm.Put(gp, fun, dbdx);
			dbdym.Put(gp, fun, dbdy);
			dbdzm.Put(gp, fun, dbdz);
		}
	}

	bool isGenAniso = (elMat->type == genAniso) || model.UseSimpleElements();

	for (gp = 0; gp < nint; gp++)
	{
		double *dbdxv = dbdxm.GetRow(gp);
		double *dbdyv = dbdym.GetRow(gp);
		double *dbdzv = dbdzm.GetRow(gp);
		w = wv[gp];
		if (!isGenAniso)
		{
			for (int rowfun = 0; rowfun < nfun; rowfun++)
				colFunLoopIsoOrtho(dbdxv, dbdyv, dbdzv, w, rowfun, btc, stiff);
		}
		else
		{
			for (int rowfun = 0; rowfun < nfun; rowfun++)
				colFunLoopGenAniso(dbdxv, dbdyv, dbdzv, w, rowfun, btc, stiff);
		}
	}




	if ( hasLcsConstraint())
		AddPenaltyToStiffnessandEnfDisp(stiff);


	stiffLength = len;

	return stiff;
}

void SRelement::AddPenaltyToStiffnessandEnfDisp(double *stiff)
{
	//handle constraints in non-gcs coordinates with penalty method
	//output:
		//stiff updated with contribution of penalty
	//notes:
		//if there are any enforced displacements, global force vector is modified
		//this routine if called if there is a non-gcs edge or face constraint on one of the boundary faces of this element

	int faceToElFuns[MAXFACEBASISFUNCTIONS];
	int condof, rowfun, colfun, rowdof, coldof, roweq, coleq, elrowfun, elcolfun;
	SRface *face;
	SRconstraint *con;
	double rf, sf, w, bi, bj, enfdisp, detj, bibj, bibjelB;
	int i, nint, nfun, symloc, nfun2;
	double *basisv = basisVec;
	SRvec3 elocal;
	double elA, elB;
	SRvec3 pos;

	double *forceVec = model.GetElementData(thread)->GetEnfdForceVec();
	for (int facon = 0; facon < faceLCSConstraints.GetNum(); facon++)
	{
		int lface = faceLCSConstraints.Get(facon);
		face = localFaces.GetPointer(lface)->GetFace();
		con = face->GetConstraint();
		if (con == NULL)
			ERROREXIT;
		SRcoord* coord = con->GetCoord();
		//map face basis function numbers to element:
		nfun = MapFaceToElFuns(lface, faceToElFuns);
		nint = model.math.FillGaussPoints(face, -1, thread);// -1 means p not being provided, use face pmax
		int gfun, geq;
		for (condof = 0; condof < 3; condof++)
		{
			if (!con->IsConstrainedDof(condof))
				continue;
			for (i = 0; i < nint; i++)
			{
				model.math.GetGP2d(i, rf, sf, w, thread);
				model.basis.FaceBasisFuncs(rf, sf, face, basisv);
				detj = face->Jacobian(rf, sf);
				w *= detj;
				double pw = penaltyStiffness *w;
				face->Position(rf, sf, pos);
				coord->CalculateLocalDirection(pos, condof, elocal);
				if (con->hasEnforcedDisp())
					enfdisp = con->GetFaceEnforcedDisp(rf, sf, condof);
				for (rowfun = 0; rowfun< nfun; rowfun++)
				{
					gfun = face->GetGlobalFunctionNumber(rowfun);
					bj = basisv[rowfun];
					bj *= pw;
					elrowfun = faceToElFuns[rowfun];
					//contribution to enforced displacement:
					if (con->hasEnforcedDisp())
					{
						for (rowdof = 0; rowdof < 3; rowdof++)
						{
							elB = elocal.d[rowdof];
							geq = model.GetFunctionEquation(gfun, rowdof);
							if (geq >= 0)
								forceVec[geq] += (enfdisp*bj*elB);
						}
					}
					for (colfun = 0; colfun < nfun; colfun++)
					{
						bi = basisv[colfun];
						elcolfun = faceToElFuns[colfun];
						bibj = bj*bi;
						for (rowdof = 0; rowdof < 3; rowdof++)
						{
							roweq = elrowfun * 3 + rowdof;
							elB = elocal.d[rowdof];
							bibjelB = bibj*elB;
							for (coldof = 0; coldof<3; coldof++)
							{
								elA = elocal.d[coldof];
								coleq = elcolfun * 3 + coldof;
								if (roweq > coleq)
									continue;
								symloc = GetStiffnessLocationAboveDiag(roweq, coleq);
								stiff[symloc] += (bibjelB*elA);
							}
						}
					}
				}
			}
		}
	}
	for (int ncon = 0; ncon < nodeLcSConstraints.GetNum(); ncon++)
	{
		int lnode = nodeLcSConstraints.Get(ncon);
		SRnode* node = getNodeOrMidNode(lnode);
		con = node->GetConstraint();
		SRcoord* coord = con->GetCoord();
		int elfun = lnode;
		if (node->isMidSide())
		{
			int eid = node->midSideEdgeOwner;
			int lej = localEdgeMatch(eid);
			elfun = MapEdgeP2FnToElFun(lej);
		}
		for (condof = 0; condof < 3; condof++)
		{
			if (!con->IsConstrainedDof(condof))
				continue;
			coord->CalculateLocalDirection(node->Position(), condof, elocal);
			if (con->hasEnforcedDisp())
				enfdisp = con->getDisp(0, condof);
			int gfun;
			if (!node->isMidSide())
				gfun = node->globalFunctionNumber;
			else
			{
				int eid = node->midSideEdgeOwner;
				SRedge* edge = model.GetEdge(eid);
				gfun = edge->GetGlobalFunctionNumber(2);
			}
			for (rowdof = 0; rowdof < 3; rowdof++)
			{
				elB = elocal.d[rowdof] * penaltyStiffness;
				//contribution to enforced displacement:
				if (con->hasEnforcedDisp())
				{
					int geq = model.GetFunctionEquation(gfun, rowdof);
					if (geq >= 0)
						forceVec[geq] += (enfdisp*elB);
				}

				roweq = elfun * 3 + rowdof;
				for (coldof = 0; coldof < 3; coldof++)
				{
					elA = elocal.d[coldof];
					coleq = elfun * 3 + coldof;
					if (roweq > coleq)
						continue;
					symloc = GetStiffnessLocationAboveDiag(roweq, coleq);
					stiff[symloc] += (elA*elB);
				}
			}
		}
	}
}

void SRelement::CalibratePenalty()
{
	//calibrate penalty stiffness for constraints in curvilinear coordinates
	//note:
		//stores class variable penaltyStiffness = stiffness of a thin layer of same material as element
		//duplicate calls ignored

	if (penaltyStiffness > TINY)
		return;
	penaltyStiffness = 0.0;
	double E = elMat->MatScale();

	//penalty multiplier. this makes penalty stiffness equivalent to
	//that of a layer, same material as element,
	//thickness = (element size)/penaltyMult
	double penaltyMult = 1.0e4;

	penaltyStiffness = (penaltyMult*E / size);
}

void SRelement::FillKel33GenAnisoRowWise(SRBTC& btc, int rowfun, int colfun, double dbdx, double dbdy, double dbdz, double kel33[3][3])
{
	//Fill 3x3 portion of elemental stiffness matrix for a general Anisotropic material,
	//compatible with rowwise storage of elements
	//input:
		//btc = strain-displacement matrixe times elasticity matrix
		//rowfun = row function number
		//colfun = column function number
		//dbdx, dbdy, dbdz = derivatives of basis function
	//output:
		//fills kel33[3][3]


	kel33[0][0] = btc.BTC11*dbdx + btc.BTC14*dbdy + btc.BTC15*dbdz;
	kel33[0][1] = btc.BTC12*dbdy + btc.BTC14*dbdx + btc.BTC16*dbdz;
	kel33[0][2] = btc.BTC13*dbdz + btc.BTC15*dbdx + btc.BTC16*dbdy;
	kel33[1][1] = btc.BTC22*dbdy + btc.BTC24*dbdx + btc.BTC26*dbdz;
	kel33[1][2] = btc.BTC23*dbdz + btc.BTC25*dbdx + btc.BTC26*dbdy;
	kel33[2][2] = btc.BTC33*dbdz + btc.BTC35*dbdx + btc.BTC36*dbdy;
	if (rowfun != colfun)
	{
		kel33[1][0] = btc.BTC21*dbdx + btc.BTC24*dbdy + btc.BTC25*dbdz;
		kel33[2][0] = btc.BTC31*dbdx + btc.BTC34*dbdy + btc.BTC35*dbdz;
		kel33[2][1] = btc.BTC32*dbdy + btc.BTC34*dbdx + btc.BTC36*dbdz;
	}
}

void SRelement::StraintoStress(double r, double s, double t, double strain[], double stress[])
{
	//calculate stress given strain at a point of this elements
	//input:
		//r,s,t = natural coordinate at this point
		//strain = strain tensor
	//output:
		//stress = stress tensor

	SRthermalForce *tf = model.GetThermalForce();
	double mechStrain[6];
	int i;
	for(i = 0; i < 6; i++)
		mechStrain[i] = strain[i];
	if(tf != NULL)
	{
		double deltaTemp, etx, ety, etz;
		deltaTemp = tf->GetTemp(this, r, s, t) - elMat->tref;
		etx = deltaTemp*elMat->alphax;
		ety = deltaTemp*elMat->alphay;
		etz = deltaTemp*elMat->alphaz;
		mechStrain[0] -= etx;
		mechStrain[1] -= ety;
		mechStrain[2] -= etz;
	}
	if (elMat->type == iso)
	{
		double c11 = elMat->c11;
		double lambda = elMat->lambda;
		double G = elMat->G;
		stress[0] = c11*mechStrain[0] + lambda*(mechStrain[1] + mechStrain[2]);
		stress[1] = c11*mechStrain[1] + lambda*(mechStrain[0] + mechStrain[2]);
		stress[2] = c11*mechStrain[2] + lambda*(mechStrain[0] + mechStrain[1]);
		stress[3] = G*mechStrain[3];
		stress[4] = G*mechStrain[4];
		stress[5] = G*mechStrain[5];
	}
	else if (elMat->type == ortho)
	{
		SRcij& cij = elMat->orthoCij;
		stress[0] = cij.c11*mechStrain[0] + cij.c12*mechStrain[1] + cij.c13*mechStrain[2];
		stress[1] = cij.c12*mechStrain[0] + cij.c22*mechStrain[1] + cij.c23*mechStrain[2];
		stress[2] = cij.c13*mechStrain[0] + cij.c23*mechStrain[1] + cij.c33*mechStrain[2];
		stress[3] = cij.c44*mechStrain[3];
		stress[4] = cij.c55*mechStrain[4];
		stress[5] = cij.c66*mechStrain[5];
	}
	else if (elMat->type == genAniso)
	{
		SRgenAnisoCij& gcij = elMat->genAnisoCij;
		for(int i = 0; i < 6; i++)
		{
			stress[i] = 0.0;
			for(int j = 0; j < 6; j++)
				stress[i] += (gcij.c[i][j] * mechStrain[j]);
		}
	}

	for (int i = 0; i < 6; i++)
	{
		double ae = fabs(stress[i]);
		if (ae > maxStress)
			maxStress = ae;
	}

}

void SRelement::NodeNaturalCoords(int lnode, double &r, double &s, double &t)
{
	//look up natural coordinate position of a local node of an element (corner or midside node)
	//input:
		//lnode = local node number
	//output:
		//r,s,t = natural coordinates in element
	if (lnode >= nodeIds.GetNum())
	{
		//midside node
		int lej = lnode - nodeIds.GetNum();
		model.map.ElementNaturalCoordsAtMidedge(this, lej, r, s, t);
		return;
	}

	if (type == tet)
	{
		if (lnode == 0)
		{
			r = -1.0;
			s = 0.0;
			t = 0.0;
		}
		else if (lnode == 1)
		{
			r = 1.0;
			s = 0.0;
			t = 0.0;
		}
		else if (lnode == 2)
		{
			r = 0.0;
			s = SQRT3;
			t = 0.0;
		}
		else if (lnode == 3)
		{
			r = 0.0;
			s = SQRT3OVER3;
			t = TWOSQRTTWOTHIRDS;
		}
		else
			ERROREXIT;
	}
	else if (type == wedge)
		GetWedgeNodePos(lnode, r, s, t);
	else
		GetBrickNodePos(lnode, r, s, t);
}

void SRelement::Clear()
{
	//free memory of an element that is no longer needed
	localEdges.Free();
	localFaces.Free();
	globalFunctionNumbers.Free();
}

void SRelement::DownloadDisplacement()
{
	//download from global solution to local displacements of this element
	//note:
		//fills up class variable dispEl
	int elfun;
	int dof, gfun, leq;
	int nfun = globalFunctionNumbers.GetNum();
	dispEl.Allocate(3 * nfun);
	double *dispElVec = dispEl.GetVector();
	for (elfun = 0; elfun < nfun; elfun++)
	{
		gfun = globalFunctionNumbers.Get(elfun);
		for (dof = 0; dof < 3; dof++)
		{
			leq = 3 * elfun + dof;
			double disp = model.GetDisplacementCoeff(gfun, dof);
			dispElVec[leq] = disp;
		}
	}
}

int SRelement::GetMaxPorder()
{
	//get max p order of this element
	int i, p;
	maxPorder = 0;
	for(i = 0; i < localEdges.GetNum(); i++)
	{
		p = localEdges.Get(i).GetPOrder();
		if(p > maxPorder)
			maxPorder = p;
	}
	return maxPorder;
}

double SRelement::CalculateRawStrain(double r, double s, double t, double strain[], double& etx, double& ety, double& etz)
{
	//calculate raw strain in an element at a point
	//input:
		//r,s,t = natural coordinates in element
		//dbasisdr, dbasisds, dbasisdt = element basis derivative vectors if already calculated, else NULL
	//output:
		//strain[6] = vector of strains at this location
		//etx, ety, etz = thermal strain at this location
	//return:
		//maximum thermal strain
	//Note: if thermal loading is present, strain is total strain
	int fun, nfun, eq;
	SRmap* map = &model.map;
	double dbdx, dbdy, dbdz, dbdr, dbds, dbdt, uel, vel, wel, ex, ey, ez, gamxy, gamxz, gamyz;

	nfun = globalFunctionNumbers.GetNum();

	ex = 0.0;
	ey = 0.0;
	ez = 0.0;
	gamxy = 0.0;
	gamxz = 0.0;
	gamyz = 0.0;
	double *dispElVec = dispEl.GetVector();
	for (fun = 0; fun < nfun; fun++)
	{
		dbdr = dbasisdr[fun];
		dbds = dbasisds[fun];
		dbdt = dbasisdt[fun];
		XyzDerivatives(dbdr, dbds, dbdt, dbdx, dbdy, dbdz);
		eq = 3 * fun;
		uel = dispElVec[eq];
		vel = dispElVec[eq + 1];
		wel = dispElVec[eq + 2];
		ex += dbdx*uel;
		ey += dbdy*vel;
		ez += dbdz*wel;
		gamxy += (dbdy*uel + dbdx*vel);
		gamxz += (dbdz*uel + dbdx*wel);
		gamyz += (dbdz*vel + dbdy*wel);
	}

	strain[0] = ex;
	strain[1] = ey;
	strain[2] = ez;
	strain[3] = gamxy;
	strain[4] = gamxz;
	strain[5] = gamyz;

	//max thermal strain component at this point:
	double eT = 0.0;
	etx = ety = etz = 0.0;
	SRthermalForce* tf = model.GetThermalForce();
	if(tf != NULL)
	{
		double deltaTemp, etx, ety, etz;
		deltaTemp = tf->GetTemp(this, r, s, t) - elMat->tref;
		double alfx, alfy, alfz;
		etx = deltaTemp*elMat->alphax;
		ety = deltaTemp*elMat->alphay;
		etz = deltaTemp*elMat->alphaz;
		eT = fabs(etx);
		ety = fabs(ety);
		etz = fabs(etz);
		if(ety > eT)
			eT = ety;
		if(etz > eT)
			eT = etz;
	}

	return eT;
}

int SRelement::MapFaceToElFuns(int lface, int FaceToElFuns[])
{
	//get basis function numbers on a face to corresponding local element basis function numbers
	//for this element that owns the face
	//input:
	//lface = local face number
	//output:
	//FaceToElFuns[face function number] = element function number
	//return:
	//number of functions on the face

	int i, lej, lfaceej, elfun, facefun = 0, gf, ge, n, elf, elfun0;
	SRface* face = GetFace(lface);
	for (i = 0; i < face->GetNumNodes(); i++)
	{
		gf = face->GetNodeId(i);
		for (elfun = 0; elfun < nodeIds.GetNum(); elfun++)
		{
			ge = nodeIds.Get(elfun);
			if (ge == gf)
			{
				FaceToElFuns[facefun] = elfun;
				facefun++;
				break;
			}
		}
	}

	elfun0 = nodeIds.GetNum();

	//p2 edge funs of the face:
	for (lfaceej = 0; lfaceej < face->GetNumLocalEdges(); lfaceej++)
	{
		gf = face->GetLocalEdgeGlobalId(lfaceej);
		int nf = face->GetLocalEdgePOrder(lfaceej) + 1;
		elfun = elfun0;
		for (lej = 0; lej < localEdges.GetNum(); lej++)
		{
			ge = localEdges.GetPointer(lej)->globalEdgeId;
			if (ge == gf)
			{
				FaceToElFuns[facefun] = elfun;
				facefun++;
			}
			elfun++;
		}
	}

	elfun0 = elfun;

	//higher edge funs of the face:
	for (lfaceej = 0; lfaceej < face->GetNumLocalEdges(); lfaceej++)
	{
		gf = face->GetLocalEdgeGlobalId(lfaceej);
		int nf = face->GetLocalEdgePOrder(lfaceej) + 1;
		elfun = elfun0;
		for (lej = 0; lej < localEdges.GetNum(); lej++)
		{
			ge = localEdges.GetPointer(lej)->globalEdgeId;
			if (ge == gf)
			{
				for (i = 3; i < nf; i++)
				{
					FaceToElFuns[facefun] = elfun;
					facefun++;
					elfun++;
				}
			}
			else
			{
				int nfl = localEdges.GetPointer(lej)->GetEdge()->GetNumGlobalFunctions() - 3;
				if (nfl > 0)
					elfun += nfl;
			}
		}
	}


	for (elf = 0; elf < lface; elf++)
	{
		SRface* elface = GetFace(elf);
		n = model.basis.CountFaceInternalFunctions(elface);
		elfun += n;
	}

	n = model.basis.CountFaceInternalFunctions(face);
	for (i = 0; i < n; i++)
	{
		FaceToElFuns[facefun] = elfun;
		facefun++;
		elfun++;
	}

	n = model.basis.CountFaceTotalFunctions(face, i);
	SRASSERT(n == facefun);

	return facefun;
}

int SRelement::MapEdgeToElFuns(int ledge, int EdgeToElFuns[])
{
	//get basis function numbers on an edge to corresponding local element basis function numbers
	//for this element that owns the edge
	//input:
	//ledge = local edge number
	//output:
	//EdgeToElFuns[edge function number] = element function number

	SRedge* edge = GetEdge(ledge);
	int edgefun = 0;
	for (int i = 0; i < 2; i++)
	{
		int gf = edge->GetNodeId(i);
		for (int elfun = 0; elfun < nodeIds.GetNum(); elfun++)
		{
			int ge = nodeIds.Get(elfun);
			if (ge == gf)
			{
				EdgeToElFuns[edgefun] = elfun;
				edgefun++;
				break;
			}
		}
	}

	int elfun = nodeIds.GetNum();
	//midedge fn:
	for (int e = 0; e < ledge; e++)
		elfun++;
	EdgeToElFuns[edgefun] = elfun;
	elfun++;
	edgefun++;

	elfun = nodeIds.GetNum() + localEdges.GetNum();

	for (int e = 0; e < ledge; e++)
	{
		int nf = localEdges.Get(e).GetEdge()->GetNumGlobalFunctions();
		if (nf > 2)
			elfun += (nf - 2);
	}
	for (int f = 2; f < edge->GetNumGlobalFunctions(); f++)
	{
		EdgeToElFuns[edgefun] = elfun;
		elfun++;
		edgefun++;
	}

	return edge->GetNumGlobalFunctions();
}

int SRelement::MapEdgeP2FnToElFun(int ledge)
{
	//get element's local edge p2 function numbers on an local edge
	//input:
	//ledge = local edge number
	//return:
	//local element function number for the edge's p2 function

	int elfun = nodeIds.GetNum() + ledge;
	return elfun;
}

void SRelement::GetSmoothedStrain(double r, double s, double t, double strain[])
{
	//look up smoothed strain in an element at a point
	//input:
		//r,s,t = natural coordinates
	//output:
		//strain[6] = vector of strains at this location
	int i, n, c;
	double* basis = basisVec;
	
	for (c = 0; c < 6; c++)
	{
		strain[c] = 0.0;
		n = model.basis.ElementBasisFuncs(r, s, t, this, basis);
		for (i = 0; i < n; i++)
			strain[c] += (smoothedStrains[c].Get(i)*basis[i]);
	}
}

void SRelement::GetStress(double r, double s, double t, double stress[], bool checkMaxStrain)
{
	//look up smoothed stress in an element at a point
	//input:
		//r,s,t = natural coordinates
	//output:
		//stress[6] = vector of strains at this location
	double strain[6];
	GetSmoothedStrain(r, s, t, strain);
	StraintoStress(r, s, t, strain, stress);
	if (checkMaxStrain)
	{
		double evm = model.math.GetSvm(strain);
		if (evm > model.GetStrainMax())
			model.SetStrainMax(evm);
	}
}

void SRelement::DownloadSmoothedStrains()
{
	//download smoothed strains from global solution to local strains of this element
	//note:
		//fills up class variable smoothedStrains

	int n = globalFunctionNumbers.GetNum();
	int c, i, f, eq;
	for (i = 0; i < 6; i++)
		smoothedStrains[i].Allocate(n);
	double *strainVec;
	for(c=0;c<6;c++)
	{
		strainVec = model.post.getSmoothedStrainVec(c);
		if(strainVec == NULL)
			return;
		for(i = 0; i < n; i++)
		{
			f = globalFunctionNumbers.Get(i);
			eq = model.GetSmoothFunctionEquation(f);
			if(eq >= 0)
				smoothedStrains[c].Put(i, strainVec[eq]);
		}
	}
}

void SRelement::FillMappingNodesBreakoutOnly(bool needShapeDerivs)
{
	//fill in mapping node vectors xnode, ynode, znode
	//from nodes of element
	//breakout version:

	//note:
	//fills class variables xnode, ynode, znode

	numNodesTotal = GetNumNodes();

	xnodev.Allocate(numNodesTotal);
	ynodev.Allocate(numNodesTotal);
	znodev.Allocate(numNodesTotal);
	xnode = xnodev.d;
	ynode = ynodev.d;
	znode = znodev.d;
	for (int i = 0; i < numNodesTotal; i++)
		NodeLocalCopy(i, GetNodeId(i));
	if (needShapeDerivs)
	{
		dNdrv.Allocate(numNodesTotal);
		dNdsv.Allocate(numNodesTotal);
		dNdtv.Allocate(numNodesTotal);
		dNdr = dNdrv.d;
		dNds = dNdsv.d;
		dNdt = dNdtv.d;
	}
}

void SRelement::FillMappingNodes()
{
	//fill in mapping node vectors xnode, ynode, znode
	//from nodes and midnodes of element

	//note:
		//fills class variables xnode, ynode, znode

	int numcorner, nummidside;
	numcorner = GetNumNodes();
	nummidside = GetNumLocalEdges();
	if (model.saveBreakout)
		numcorner -= nummidside;
	numNodesTotal = numcorner + nummidside;

	xnodev.Allocate(numNodesTotal);
	ynodev.Allocate(numNodesTotal);
	znodev.Allocate(numNodesTotal);
	xnode = xnodev.d;
	ynode = ynodev.d;
	znode = znodev.d;
	dNdrv.Allocate(numNodesTotal);
	dNdsv.Allocate(numNodesTotal);
	dNdtv.Allocate(numNodesTotal);
	dNdr = dNdrv.d;
	dNds = dNdsv.d;
	dNdt = dNdtv.d;

	int i;
	for (i = 0; i < numcorner; i++)
		NodeLocalCopy(i, GetNodeId(i));
	for (i = 0; i < nummidside; i++)
		NodeLocalCopy(i + numcorner, GetLocalEdgeMidNodeId(i));

	//Element size; approximate as average length of all element
	//edges:
	size = 0.0;
	SRnode* n0;
	SRnode* n1;
	for (i = 0; i < nummidside; i++)
	{
		n0 = GetLocalEdgeNode(i, 0);
		n1 = GetLocalEdgeNode(i, 1);
		size += model.math.Distance(n0->Position(), n1->Position());
	}
	size /= nummidside;

	//approx volume is volume of BB. use as scale in Jacobian checks:
	double x, y, z, xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = BIG;
	xmax = -BIG;
	ymin = BIG;
	ymax = -BIG;
	zmin = BIG;
	zmax = -BIG;
	for (i = 0; i < numcorner; i++)
	{
		x = xnode[i];
		y = ynode[i];
		z = znode[i];
		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
		if (z < zmin)
			zmin = z;
		if (z > zmax)
			zmax = z;
	}
	approxVol = (xmax - xmin)*(ymax - ymin)*(zmax - zmin);
	if (approxVol < TINY)
		ERROREXIT;
}

void SRelement::approxCentroid(SRvec3& elcen)
{
	//approx centroid is average corner position
	//use this when mapping not setup
	int nn = GetNumCorners();
	for (int i = 0; i < nn; i++)
	{
		SRnode* node = GetNode(i);
		elcen.PlusAssign(node->pos);
	}
	double scale = 1.0 / nn;
	elcen.Scale(scale);
}


void SRelement::NodeLocalCopy(int lnode, int gnode)
{
	//copy nodal coordinates into local xnode, ynode, znode vectors
	//input:
		//lnode = element local node number
		//gnode = global node id
	//note:
		//fills class variables xnode, ynode, znode

	SRnode* node;
	node = model.GetNode(gnode);
	xnode[lnode] = node->GetXyz(0);
	ynode[lnode] = node->GetXyz(1);
	znode[lnode] = node->GetXyz(2);
}

double SRelement::FillMapping(double r, double s, double t, bool detJonly)
{
	//Calculate mapping at r,s,t for an element
	//input:
		//r,s,t = natural coordinates
		//detJonly = true if only need detJ else false
	//return:
		//determinant of Jacobian of mapping
	//note:
		//fill class variable elJac.

	FillJacobian(r, s, t);

	if (detJonly)
		return elJac.Determinant();
	//invert mapping:
	double detj = elJac.Invert();
	return detj;
}

void SRelement::FillJacobian(double r, double s, double t, bool fillElJac)
{
	//fill jacobian for an element; support routine for FillMapping
	//input:
		//r,s,t = natural coordinates in mapElem
		//fillEljac = true to fill jacobian (eljac), false to just fill mapping derivatives
	//note:
		//fills class variable elJac = jacobian

	if (GetType() == tet)
		FillTetJacobian(r, s, t, fillElJac);
	else if (GetType() == wedge)
		FillWedgeJacobian(r, s, t, fillElJac);
	else if (GetType() == brick)
		FillBrickJacobian(r, s, t, fillElJac);
}

void SRelement::FillTetJacobian(double r, double s, double t, bool fillElJac)
{
	//calculate the jacobian of a tet
	//input:
		//r,s,t = natural coordinate
		//fillEljac = true to fill jacobian (eljac), false to just fill mapping derivatives
	//note:
		//fill class variables dNdr, dNds, dNdt, and elJac.

	double L[4];
	int i;
	model.map.TetVolumeCoords(r, s, t, L);

	//corner mapping functions:
	for (i = 0; i < 4; i++)
	{
		dNdr[i] = model.map.GetDLtdr(i) * (4.0*L[i] - 1.0);
		dNds[i] = model.map.GetDLtds(i) * (4.0*L[i] - 1.0);
		dNdt[i] = model.map.GetDLtdt(i) * (4.0*L[i] - 1.0);
	}


	//mid-edge mapping functions:
	double *dLtdr = model.map.GetDLtdrVec();
	double *dLtds = model.map.GetDLtdsVec();
	double *dLtdt = model.map.GetDLtdtVec();

	dNdr[4] = model.map.TetMidEdgeMappingDeriv(0, 1, dLtdr, L);
	dNds[4] = model.map.TetMidEdgeMappingDeriv(0, 1, dLtds, L);
	dNdt[4] = model.map.TetMidEdgeMappingDeriv(0, 1, dLtdt, L);
	dNdr[5] = model.map.TetMidEdgeMappingDeriv(1, 2, dLtdr, L);
	dNds[5] = model.map.TetMidEdgeMappingDeriv(1, 2, dLtds, L);
	dNdt[5] = model.map.TetMidEdgeMappingDeriv(1, 2, dLtdt, L);
	dNdr[6] = model.map.TetMidEdgeMappingDeriv(0, 2, dLtdr, L);
	dNds[6] = model.map.TetMidEdgeMappingDeriv(0, 2, dLtds, L);
	dNdt[6] = model.map.TetMidEdgeMappingDeriv(0, 2, dLtdt, L);
	dNdr[7] = model.map.TetMidEdgeMappingDeriv(0, 3, dLtdr, L);
	dNds[7] = model.map.TetMidEdgeMappingDeriv(0, 3, dLtds, L);
	dNdt[7] = model.map.TetMidEdgeMappingDeriv(0, 3, dLtdt, L);
	dNdr[8] = model.map.TetMidEdgeMappingDeriv(1, 3, dLtdr, L);
	dNds[8] = model.map.TetMidEdgeMappingDeriv(1, 3, dLtds, L);
	dNdt[8] = model.map.TetMidEdgeMappingDeriv(1, 3, dLtdt, L);
	dNdr[9] = model.map.TetMidEdgeMappingDeriv(2, 3, dLtdr, L);
	dNds[9] = model.map.TetMidEdgeMappingDeriv(2, 3, dLtds, L);
	dNdt[9] = model.map.TetMidEdgeMappingDeriv(2, 3, dLtdt, L);

	if (fillElJac)
		FillJac();
}

void SRelement::FillJac()
{
	//fill the jacobian at a natural coordinate point in an element
	//note:
	//class variables xnode, ynode, znode and dNdr, dNds, dNdr must be filled first
	//fills the class variable elJac

	elJac.Zero();
	int i;
	for (i = 0; i< numNodesTotal; i++)
	{
		elJac.rows[0].d[0] += xnode[i] * dNdr[i]; //dxdr
		elJac.rows[1].d[0] += xnode[i] * dNds[i]; //dxds
		elJac.rows[2].d[0] += xnode[i] * dNdt[i]; //dxdt
		elJac.rows[0].d[1] += ynode[i] * dNdr[i]; //dydr
		elJac.rows[1].d[1] += ynode[i] * dNds[i]; //dyds
		elJac.rows[2].d[1] += ynode[i] * dNdt[i]; //dydt
		elJac.rows[0].d[2] += znode[i] * dNdr[i]; //dzdr
		elJac.rows[1].d[2] += znode[i] * dNds[i]; //dzds
		elJac.rows[2].d[2] += znode[i] * dNdt[i]; //dzdt
	}
}

int SRelement::ShapeFunctionsLinear(double r, double s, double t, double N[])
{
	if (type == tet)
	{
		model.map.TetVolumeCoords(r, s, t, N);
		return 4;
	}
	else if (type == wedge)
	{
		model.map.WedgeLinearShapeFunctions(r, s, t, N);
		return 6;
	}
	else
	{
		model.map.BrickLinearShapeFunctions(r, s, t, N);
		return 8;
	}
}

void SRelement::PositionLinear(double r, double s, double t, SRvec3 &p)
{
	//determine element position at r,s,t using quadratic mapping
	//input:
	//r,s,t = natural coordinates in element
	//output:
	//p = position vector

	int i;
	double N[8];
	int nn = ShapeFunctionsLinear(r, s, t, N);
	p.Zero();
	for (i = 0; i < nn; i++)
	{
		p.d[0] += N[i] * xnode[i];
		p.d[1] += N[i] * ynode[i];
		p.d[2] += N[i] * znode[i];
	}
}

void SRelement::Position(double r, double s, double t, SRvec3 &p)
{
	//determine element position at r,s,t using quadratic mapping
	//input:
		//r,s,t = natural coordinates in element
	//output:
		//p = position vector

	int i;
	double N[20];
	ShapeFunctions(r, s, t, N);
	p.Zero();
	for (i = 0; i < numNodesTotal; i++)
	{
		p.d[0] += N[i] * xnode[i];
		p.d[1] += N[i] * ynode[i];
		p.d[2] += N[i] * znode[i];
	}
}

int SRelement::ShapeFunctions(double r, double s, double t, double N[])
{
	//calculate the shape functions for an element at a natural coordinate point
	//input:
		//r, s, t = natural coordinates
	//output:
	//N = shape functions
	//return:
		//number of functions

	double L[4], l04, l14;
	int i;
	if (GetType() == tet)
	{
		model.map.TetVolumeCoords(r, s, t, L);
		for (i = 0; i < 4; i++)
			N[i] = L[i] * (2.0*L[i] - 1.0);
		l04 = 4.0*L[0];
		l14 = 4.0*L[1];
		N[4] = l04*L[1];
		N[5] = l14*L[2];
		N[6] = l04*L[2];
		N[7] = l04*L[3];
		N[8] = l14*L[3];
		N[9] = 4.0*L[2] * L[3];
		return 10;
	}
	else if (GetType() == wedge)
		return WedgeShapeFns(r, s, t, N);
	else
		return BrickShapeFns(r, s, t, N);
}

void SRelement::XyzDerivatives(double dfdr, double dfds, double dfdt, double &dfdx, double &dfdy, double &dfdz)
{
	//Calculate x,y,z derivatives of a function given r,s,t derivatives and mapping
	//input:
		//dfdr,dfds,dfdt = derivatives of function in natural coordinates
	//output:
		//dfdx,dfdy,dfdz = derivatives of function in global coordinates
	//note:
		//must first have called FillMapping to
		//store  jacobian inverse j11, etc, in "elJac"

	dfdx = (elJac.rows[0].d[0] * dfdr) + (elJac.rows[0].d[1] * dfds) + (elJac.rows[0].d[2] * dfdt);
	dfdy = (elJac.rows[1].d[0] * dfdr) + (elJac.rows[1].d[1] * dfds) + (elJac.rows[1].d[2] * dfdt);
	dfdz = (elJac.rows[2].d[0] * dfdr) + (elJac.rows[2].d[1] * dfds) + (elJac.rows[2].d[2] * dfdt);
}

double* SRelement::FillBasisFuncs(double r, double s, double t, SRBasisCallType calltype)
{
	//calculate basis functions and derivatives
	//input:
		//r,s,t = natural coordinates
	//return:
		//basv = pointer to basisVec vector
	//note: fills class variables basisVec, etc
	model.basis.ElementBasisFuncs(r, s, t, this, basisVec, dbasisdr, dbasisds, dbasisdt, calltype);
	return basisVec;
}

void SRelement::FreeRawStrains()
{
	rawStrains.Free();
	dispEl.Free();
};


void SRelement::Cleanup()
{
	nodeIds.Free();
	localEdges.Free();
	localFaces.Free();
	globalFunctionNumbers.Free();
	rawStrains.Free();
	for (int i = 0; i < 6; i++)
		smoothedStrains[i].Free();
	nodeLcSConstraints.Free();
	faceLCSConstraints.Free();
	stiffnessMatrix.Free();
	stiffDiag.Free();
	dispEl.Free();
}

void SRelement::SetSacrificial(int pf)
{
	//set the sacrificial p order of this element to pf
	//note:
		//if this is not a sacrificial element, pf is left at 0, else
		//it is the p order the element will be frozen to
	sacrificial = pf;
}

void SRelement::FillStiffDiag(int neq)
{
	//fill the element storage ids of the diagonal of this element
	//note:
		//fills class variable stiffDiag
	if (stiffDiag.num == neq)
		return;
	stiffDiag.Allocate(neq);
	int nz = 0;
	for (int r = 0; r < neq; r++)
	{
		stiffDiag.Put(r, nz);
		nz += (neq - r);
	}
}

int SRelement::GetStiffnessLocationAboveDiag(int row, int col)
{
	//given the row, column of a stiffness matrix, return the storage location
	//for a symmetric element stored as upper triangle
	//note:
		//only call this if row, col are on or above the diagonal of the element, 
		//else call the more general function GetStiffnessLocation
	return stiffDiag.Get(row) + col - row;
}

int SRelement::GetStiffnessLocation(int row, int col)
{
	//given the row, column of a stiffness matrix, return the storage location
	//for a symmetric element stored as upper triangle
	//note:
		//to avoid if statement in an inner loop, call GetStiffnessLocationAboveDiag
		//if row, col are on or above the diagonal of the element
	if (col >= row)
		return stiffDiag.Get(row) + col - row;
	else
		return stiffDiag.Get(col) + row - col;
}

void SRelement::colFunLoopIsoOrtho(double* dbdxv, double* dbdyv, double* dbdzv, double w, int rowfun, SRBTC& btc, double *stiff)
{
	//implement code inside "colfun" loop for calculating stiffness matrix
	//isotropic/orthotropic material version
	//input:
		// dbdxv, dbdyv, dbdzv: vectors containing derivatives of all basis functions w.r.t. x,y,z at current gauss point
		// w = integration weight at current gauss point
		// rowfun = current row function number
	//modified:
		// stiff = element stiffness matrix
	//scratch:
		//btc = "BT*C" = transpose of strain-displacement times elasticity matrix for current row function
	//note:
		//this uses "vanilla" colfun loops.
		//I tried using cblas_daxpy instead but it is slower- have to use inc 3 for elstiff updates 
		//separate tuning tests showed daxpy for inc 3 becomes faster for vector length > 300
		//for p8 brick, nfun = 192 = max vector length, avg length is 96. so daxpy is slower

	double dbdxi = dbdxv[rowfun];
	double dbdyi = dbdyv[rowfun];
	double dbdzi = dbdzv[rowfun];
	FillBTC(w, dbdxi, dbdyi, dbdzi, btc);

	int nfun = GetNumFunctions();
	int colfun;
	int rowcolloc;
	double a, b, c;

	//rowdof 0:
	int roweq = rowfun * 3;
	int rowcolloc0 = stiffDiag.Get(roweq);
	//coldof 0
	rowcolloc = rowcolloc0;
	a = btc.BTC11;
	b = btc.BTC14;
	c = btc.BTC15;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a*dbdxv[colfun] + b*dbdyv[colfun] + c*dbdzv[colfun]);
		rowcolloc += 3;
	}
	//coldof 1
	rowcolloc = rowcolloc0 + 1;
	a = btc.BTC12;
	b = btc.BTC14;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a*dbdyv[colfun] + b*dbdxv[colfun]);
		rowcolloc += 3;
	}
	//coldof 2
	rowcolloc = rowcolloc0 + 2;
	a = btc.BTC13;
	b = btc.BTC15;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a*dbdzv[colfun] + b*dbdxv[colfun]);
		rowcolloc += 3;
	}

	//rowdof 1:
	//note: 1st coldof is skipped for colfun = rowfun, rowdof 1 because it is below diagonal
	roweq++;
	rowcolloc0 = stiffDiag.Get(roweq);
	//coldof 0
	//stiffdiag points to coldof 1 for this fun, so skip 2 to get to coldof 0 for rowfun + 1:
	rowcolloc = rowcolloc0 + 2;
	a = btc.BTC21;
	b = btc.BTC24;
	for (colfun = rowfun + 1; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a*dbdxv[colfun] + b*dbdyv[colfun]);
		rowcolloc += 3;
	}
	//coldof 1
	rowcolloc = rowcolloc0;
	a = btc.BTC22;
	b = btc.BTC24;
	c = btc.BTC26;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a*dbdyv[colfun] + b*dbdxv[colfun] + c*dbdzv[colfun]);
		rowcolloc += 3;
	}
	//coldof 2
	rowcolloc = rowcolloc0 + 1;
	a = btc.BTC23;
	b = btc.BTC26;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a*dbdzv[colfun] + b*dbdyv[colfun]);
		rowcolloc += 3;
	}

	//rowdof 2:
	//note: 1st and 2nd coldofs are skipped for colfun = rowfun, rowdof 2 because they are below diagonal
	roweq++;
	rowcolloc0 = stiffDiag.Get(roweq);
	//coldof 0
	//stiffdiag points to coldof 2 for this fun, so skip 1 to get to coldof 0 for rowfun + 1:
	rowcolloc = rowcolloc0 + 1;
	a = btc.BTC31;
	b = btc.BTC35;
	for (colfun = rowfun + 1; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a*dbdxv[colfun] + b*dbdzv[colfun]);
		rowcolloc += 3;
	}
	//coldof 1
	//stiffdiag points to coldof 2 for this fun, so skip 2 to get to coldof 1 for rowfun + 1:
	rowcolloc = rowcolloc0 + 2;
	a = btc.BTC32;
	b = btc.BTC36;
	for (colfun = rowfun + 1; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a*dbdyv[colfun] + b*dbdzv[colfun]);
		rowcolloc += 3;
	}
	//coldof 2
	rowcolloc = rowcolloc0;
	a = btc.BTC33;
	b = btc.BTC35;
	c = btc.BTC36;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a*dbdzv[colfun] + b*dbdxv[colfun] + c*dbdyv[colfun]);
		rowcolloc += 3;
	}
}

void SRelement::colFunLoopGenAniso(double* dbdxv, double* dbdyv, double* dbdzv, double w, int rowfun, SRBTC& btc, double *stiff)
{
	//implement code inside "colfun" loop for calculating stiffness matrix
	//isotropic material version
	//input:
		// dbdxv, dbdyv, dbdzv: vectors containing derivatives of all basis functions w.r.t. x,y,z at current gauss point
		// w = integration weight at current gauss point
		// rowfun = current row function number
	//modified:
		// stiff = element stiffness matrix
	//scratch:
		//btc = "BT*C" = transpose of strain-displacement times elasticity matrix for current row function
	double kel33[3][3];

	double dbdxi = dbdxv[rowfun];
	double dbdyi = dbdyv[rowfun];
	double dbdzi = dbdzv[rowfun];
	FillBTC(w, dbdxi, dbdyi, dbdzi, btc);
	int nfun = GetNumFunctions();
	int colfun;
	int rowcolloc;
	for (int colfun = rowfun; colfun < nfun; colfun++)
	{
		double dbdxj = dbdxv[colfun];
		double dbdyj = dbdyv[colfun];
		double dbdzj = dbdzv[colfun];
		FillKel33GenAnisoRowWise(btc, rowfun, colfun, dbdxj, dbdyj, dbdzj, kel33);

		int roweq = rowfun * 3;
		for (int rowdof = 0; rowdof < 3; rowdof++)
		{
			int coldof0 = 0;
			//colfun = rowfun, special case because LT of Kel33 not stored:
			if (colfun == rowfun)
				coldof0 = rowdof;
			int rowcolloc = stiffDiag.Get(roweq);
			for (int coldof = coldof0; coldof < 3; coldof++)
			{
				int coleq = colfun * 3 + coldof;
				rowcolloc = GetStiffnessLocation(roweq, coleq);
				stiff[rowcolloc] += kel33[rowdof][coldof];
			}
			roweq++;
		}
	}
}

void SRelement::FillBTC(double intwt, double dbdx, double dbdy, double dbdz, SRBTC& btc)
{
	//fill matrix "BTC" = strain-displacement matrixe times elasticity matrix
	//input:
		//intwt = integration weight at this gauss point
		//dbdx,dbdy,dbdz = derivatives of basis functions at this element function and gauss point
	//output:
		//fills BTC = B-T*C
	if (elMat->type == iso)
	{
		double c11 = intwt * elMat->c11;
		double lambda = intwt * elMat->lambda;
		double G = intwt * elMat->G;
		btc.BTC11 = c11*dbdx;
		btc.BTC12 = btc.BTC13 = lambda*dbdx;
		btc.BTC14 = G*dbdy;
		btc.BTC15 = G*dbdz;
		btc.BTC21 = btc.BTC23 = lambda*dbdy;
		btc.BTC22 = c11*dbdy;
		btc.BTC24 = G*dbdx;
		btc.BTC26 = G*dbdz;
		btc.BTC31 = btc.BTC32 = lambda*dbdz;
		btc.BTC33 = c11*dbdz;
		btc.BTC35 = G*dbdx;
		btc.BTC36 = G*dbdy;
	}
	else if (elMat->type == ortho)
	{
		SRcij& cij = elMat->orthoCij;
		btc.BTC11 = intwt * cij.c11*dbdx;
		btc.BTC12 = intwt * cij.c12*dbdx;
		btc.BTC13 = intwt *  cij.c13*dbdx;
		btc.BTC14 = intwt *  cij.c44*dbdy;
		btc.BTC15 = intwt *  cij.c55*dbdz;
		btc.BTC21 = intwt *  cij.c12*dbdy;
		btc.BTC22 = intwt *  cij.c22*dbdy;
		btc.BTC23 = intwt *  cij.c23*dbdy;
		btc.BTC24 = intwt *  cij.c44*dbdx;
		btc.BTC26 = intwt *  cij.c66*dbdz;
		btc.BTC31 = intwt *  cij.c13*dbdz;
		btc.BTC32 = intwt *  cij.c23*dbdz;
		btc.BTC33 = intwt *  cij.c33*dbdz;
		btc.BTC35 = intwt *  cij.c55*dbdx;
		btc.BTC36 = intwt *  cij.c66*dbdy;
	}
	else
	{
		SRgenAnisoCij& gcij = elMat->genAnisoCij;
		btc.BTC11 = intwt * (gcij.c[0][0] * dbdx + gcij.c[3][0] * dbdy + gcij.c[4][0] * dbdz);
		btc.BTC12 = intwt * (gcij.c[0][1] * dbdx + gcij.c[3][1] * dbdy + gcij.c[4][1] * dbdz);
		btc.BTC13 = intwt * (gcij.c[0][2] * dbdx + gcij.c[3][2] * dbdy + gcij.c[4][2] * dbdz);
		btc.BTC14 = intwt * (gcij.c[0][3] * dbdx + gcij.c[3][3] * dbdy + gcij.c[4][3] * dbdz);
		btc.BTC15 = intwt * (gcij.c[0][4] * dbdx + gcij.c[3][4] * dbdy + gcij.c[4][4] * dbdz);
		btc.BTC16 = intwt * (gcij.c[0][5] * dbdx + gcij.c[3][5] * dbdy + gcij.c[4][5] * dbdz);
		btc.BTC21 = intwt * (gcij.c[1][0] * dbdy + gcij.c[3][0] * dbdx + gcij.c[5][0] * dbdz);
		btc.BTC22 = intwt * (gcij.c[1][1] * dbdy + gcij.c[3][1] * dbdx + gcij.c[5][1] * dbdz);
		btc.BTC23 = intwt * (gcij.c[1][2] * dbdy + gcij.c[3][2] * dbdx + gcij.c[5][2] * dbdz);
		btc.BTC24 = intwt * (gcij.c[1][3] * dbdy + gcij.c[3][3] * dbdx + gcij.c[5][3] * dbdx);
		btc.BTC25 = intwt * (gcij.c[1][4] * dbdy + gcij.c[3][4] * dbdx + gcij.c[5][4] * dbdx);
		btc.BTC26 = intwt * (gcij.c[1][5] * dbdy + gcij.c[3][5] * dbdz + gcij.c[5][5] * dbdz);
		btc.BTC31 = intwt * (gcij.c[2][0] * dbdz + gcij.c[4][0] * dbdx + gcij.c[5][0] * dbdy);
		btc.BTC32 = intwt * (gcij.c[2][1] * dbdz + gcij.c[4][1] * dbdx + gcij.c[5][1] * dbdy);
		btc.BTC33 = intwt * (gcij.c[2][2] * dbdz + gcij.c[4][2] * dbdx + gcij.c[5][2] * dbdy);
		btc.BTC34 = intwt * (gcij.c[2][3] * dbdz + gcij.c[4][3] * dbdx + gcij.c[5][3] * dbdy);
		btc.BTC35 = intwt * (gcij.c[2][4] * dbdz + gcij.c[4][4] * dbdx + gcij.c[5][4] * dbdy);
		btc.BTC36 = intwt * (gcij.c[2][5] * dbdz + gcij.c[4][5] * dbdx + gcij.c[5][5] * dbdy);
	}
}

void SRelement::rescaleFlattenFraction()
{
	if (!flattened)
		return;
	flattenFraction = 0.0;
	for (int l = 0; l < localEdges.GetNum(); l++)
	{
		SRedge* edge = GetEdge(l);
		if (edge->straightenFraction > flattenFraction)
			flattenFraction = edge->straightenFraction;
	}
}

bool SRelement::checkForStraightenedEdges()
{
	//see if this element has any edges that were straightened to fix bad mapping
	//return:
	//true if any edges were straightened else false
	//note:
	//if any edges were straightened, the element is treated as sacrificialPorder and frozen at p3
	//since the geometry is not being followed properly, solution is inaccurate in this element regardless of p order
	for (int l = 0; l < localEdges.GetNum(); l++)
	{
		SRedge* edge = GetEdge(l);
		if (edge->straightenFraction > approxVol*SMALL)
		{
			flattened = true;
			if (edge->straightenFraction > flattenFraction)
				flattenFraction = edge->straightenFraction;
		}
		if (edge->straightenFraction > 0.2)
		{
			FillMappingNodes();
			sacrificial = 3; //this is set to freeze order of 3 not 2 in case these edges might be near a hot spot
			return true;
		}
	}

	if (flattenFraction > model.maxFlattened)
	{
		model.maxFlattened = flattenFraction;
		model.maxFlattenedElUId = userId;
	}
	//also refill the mapping nodes of the element's faces:
	for (int i = 0; i < localFaces.GetNum(); i++)
	{
		SRface* face = GetFace(i);
		face->FillMappingNodes();
	}

	return false;
}

bool SRelement::testMapping(bool okToFail)
{
	minJac = BIG;
	double fractions[7] = { 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 };
	int nint = model.math.FillGaussPoints(this);
	double maxStraightenFraction = 0.0;
	bool flattenedThisPass = false;
	int maxTrys = 7;
	bool mappingOk = true;
	if (!model.partialFlatten)
	{
		maxTrys = 1;
		fractions[0] = 1.0;
	}
	for (int gp = 0; gp < nint; gp++)
	{
		double r, s, t, w;
		model.math.GetGP3d(gp, r, s, t, w);
		double detj = FillMapping(r, s, t, true);//detj only is true
		if (detj < minJac)
			minJac = detj;
		if (detj < approxVol*SMALL)
		{
			mappingOk = false;
			flattened = true;
			flattenedThisPass = true;
			if (!okToFail)
			{
				OUTPRINT(" bad Jacobian. Element: %d", userId);
				ERROREXIT;
			}

			double straightenFraction;
			for (int trynum = 0; trynum < maxTrys; trynum++)
			{
				straightenFraction = fractions[trynum];
				if (straightenFraction > maxStraightenFraction)
					maxStraightenFraction = straightenFraction;
				for (int l = 0; l < localEdges.GetNum(); l++)
				{
					SRedge* edge = GetEdge(l);
					edge->Straighten(straightenFraction);
				}
				FillMappingNodes();
				detj = FillMapping(r, s, t, true);//detj only is true
				if (detj > approxVol*SMALL)
					break;
			}
		}
	}
	if (flattenedThisPass)
	{
		//OUTPRINT(" straightened Element: %d %lg", userId, maxStraightenFraction);
		flattenFraction = maxStraightenFraction;
	}
	return mappingOk;
}

void SRelement::SetBasisData(int processorNum )
{
	//set pointers to store basis function data for this element
	//input:
		//processerNum = processor that is calling the element stiffness routing, 0 if not multi-threaded
	SRElementData *eldata = model.GetElementData(processorNum);
	basisVec = eldata->GetBasisVec();
	dbasisdr = eldata->Getdbasisdr();
	dbasisds = eldata->Getdbasisds();
	dbasisdt = eldata->Getdbasisdt();
}

SRvec3 SRelement::GetDisp(double r, double s, double t)
{
	//calculate the displacement at a natural coordinate in an element
	//input:
		//r,s,t = natural coordinates
	//return
		//disp = 3 dof vector of displacement at r,s,t
	//note:
		//DownloadDisplacement must be called before entering loop that calls this routine
		
	SRElementData *eldata = model.GetElementData(0);
	basisVec = eldata->GetBasisVec();
	double *dispElVec = dispEl.GetVector();
	int nfun = globalFunctionNumbers.GetNum();
	FillBasisFuncs(r, s, t, basisonly);
	SRvec3 disp;
	for (int fun = 0; fun < nfun; fun++)
	{
		double bf = basisVec[fun];
		int gfun = globalFunctionNumbers.Get(fun);
		for (int dof = 0; dof < 3; dof++)
		{
			double u = model.GetDisplacementCoeff(gfun, dof);
			disp.d[dof] += u*bf;
		}
	}
	return disp;
}

bool SRelement::nodeDistCheck(SRvec3& pos, double radius)
{
	//check the distance between all nodes of the element and a position.
	//input:
		//pos = position
		//radius = test radius
	//return
		//true if any nodes are within radius of pos, else false
	bool anYNodeInside = false;
	for (int n = 0; n < numNodesTotal; n++)
	{
		double d = pos.Distance(xnode[n], ynode[n], znode[n]);
		if (d <= radius)
		{
			anYNodeInside = true;
			break;
		}
	}
	return anYNodeInside;
}

bool SRelement::InsideBoundingBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
	//check whether all nodes of the element are withn a bounding bix
	//input:
		//xmin, xmax, ymin, ymax, zmin, zmax = bounding box
	//return
		//true if all nodes are within bounding box, else false

	for (int i = 0; i < numNodesTotal; i++)
	{
		double x = xnode[i];
		double y = ynode[i];
		double z = znode[i];
		if (x < xmin || x > xmax || y < ymin || y > ymax || z < zmin || z > zmax)
			return false;
	}
	return true;
}

void SRelement::dumpData()
{
	SRfile& f = model.dumpFile;
	f.Print("%d %d nodes: ", id, userId);
	for (int i = 0; i < nodeIds.GetNum(); i++)
		f.Print(" %d", nodeIds.Get(i));
	for (int i = 0; i < localEdges.GetNum(); i++)
		f.Print(" %d", GetEdge(i)->GetMidNodeId());
	f.PrintReturn();
}

void SRelement::checkMaterialStressMax(double s)
{
	if (isSacrificial())
		return;
	if (s > elMat->GetMaxStress())
		elMat->PutMaxStress(s);
}

void SRelement::checkMaterialStrainMax(double e)
{
	if (sacrificial)
		return;
	if (e > elMat->GetMaxStrain())
		elMat->PutMaxStrain(e);
}

int SRelement::GetNumCorners()
{
	if (type == brick)
		return 8;
	else if (type == tet)
		return 4;
	else
		return 6;
}

int SRelement::GetLocalCornerIdFromGlobalFaceNodeId(int faceCornerId)
{
	for (int i = 0; i < GetNumCorners(); i++)
	{
		if (GetNodeId(i) == faceCornerId)
			return i;
	}
	return -1;
}

double SRelement::minCornerDistance(SRvec3& pos)
{
	double minDist = BIG;
	for (int i = 0; i < GetNumCorners(); i++)
	{
		double dist = pos.Distance(GetNode(i)->pos);
		if (dist < minDist)
			minDist = dist;
	}
	return minDist;
}

double SRelement::centroidDistance(SRvec3& pos)
{
	SRvec3 elcen;
	approxCentroid(elcen);
	return elcen.Distance(pos);
}

int SRelement::getNodeOrMidNodeId(int n)
{
	int ncorner = nodeIds.GetNum();
	if (n < ncorner)
		return nodeIds.Get(n);
	else
	{
		int lej = n - ncorner;
		return GetLocalEdgeMidNodeId(lej);
	}
}

SRnode* SRelement::getNodeOrMidNode(int n)
{
	int nid = getNodeOrMidNodeId(n);
	return model.GetNode(nid);
}

int SRelement::localEdgeMatch(int gId)
{
	for (int i = 0; i < localEdges.GetNum(); i++)
	{
		if (GetLocalEdgeGlobalId(i) == gId)
			return i;
	}
	ERROREXIT;
	return -1;
}

bool SRelement::checkSlopeKink()
{
	double dotmin = 1.0;
	SRface* facemin = NULL;
	for (int f = 0; f < localFaces.GetNum(); f++)
	{
		SRface* face = GetFace(f);
		if (!face->IsBoundaryFace())
			continue;
		for (int e = 0; e < face->localEdges.GetNum(); e++)
		{

			SRedge* edge = face->GetEdge(e);
			SRvec3 norm, pos1, v21;
			double rf, sf;
			face->NaturalCoordinatesNearMidedge(e, rf, sf);
			face->Position(rf, sf, pos1);
			face->OutwardNormal(rf, sf, norm);
			for (int bf = 0; bf < edge->boundaryfaceData.GetNum(); bf++)
			{
				int fid = edge->boundaryfaceData.Get(bf).faceId;
				if (fid == face->id)
					continue;
				int lej = edge->boundaryfaceData.Get(bf).localEdgeId;
				SRface* face2 = model.GetFace(fid);
				SRvec3 norm2, pos2;
				face2->NaturalCoordinatesNearMidedge(lej, rf, sf);
				face2->Position(rf, sf, pos2);
				face2->OutwardNormal(rf, sf, norm2);
				pos1.Subtract(pos2, v21);
				v21.Normalize();
				double normdot = norm.Dot(norm2);
				if (normdot < 0.0)
					continue;
				if ((norm2.Dot(v21) > 0.0) && (norm.Dot(v21) < 0.0))
				{
					//reentrant if vector from point on face2 to point on face 1
					//points in same direction as normal to face2 and in opposite direction to normal to face1
					if (normdot < dotmin)
					{
						dotmin = normdot;
						facemin = face;
					}
				}
			}
		}
	}

	if (dotmin < 0.984)//tol ~10 degrees
	{
		bool plotmaxface = true;
		SCREENPRINT("slopekinkatmax check. eluid %d dotmin: %lg\n", this->userId, dotmin);
		OUTPRINT("slopekinkatmax check. eluid %d dotmin: %lg\n", this->userId, dotmin);
		if (plotmaxface)
		{
			model.post.PlotFaces(1, &(facemin->id));
			model.post.PlotElems(1, &(this->id));
		}
		return true;
	}
	else
		return false;
}

void SRelement::checkAllNodesHaveDisp()
{
	allNodesHaveDisp = true;
	for (int i = 0; i < GetNumNodes(); i++)
	{
		SRnode* node = this->GetNode(i);
		if (!node->hasDisp)
		{
			allNodesHaveDisp = false;
			break;
		}
	}
}
bool SRelement::checkNearBktBdry()
{
	//for breakout model, see if element is near the boundary for purposes
	//of filtering stresses.
	//"near" is defined as: element, or any of it's neighbors, has a node on bdry
	for (int n = 0; n < numNodesTotal; n++)
	{
		SRnode* node = getNodeOrMidNode(n);
		if (node->checkUnsup())
			return true;
	}

	int adjelemlist[8];
	int nface = GetNumLocalFaces();
	int nadjelem = 0;
	for (int f = 0; f < nface; f++)
	{
		SRface* face = GetFace(f);
		int elemid2 = face->elementOwners[1];
		if (elemid2 == -1)
			continue;
		int elemid1 = face->elementOwners[0];
		if (elemid1 != id)
			elemid2 = elemid1;
		SRelement *elem2 = model.GetElement(elemid2);
		for (int n = 0; n < numNodesTotal; n++)
		{
			SRnode* node = elem2->getNodeOrMidNode(n);
			if (node->checkUnsup())
				return true;
		}
		adjelemlist[nadjelem] = elem2->id;
	}

	//repeat 1 more level, face neighbor of the adjacent elems:
	for (int a = 0; a < nadjelem; a++)
	{
		SRelement* elem = model.GetElement(adjelemlist[a]);
		for (int f = 0; f < elem->GetNumLocalFaces(); f++)
		{
			SRface* face = GetFace(f);
			int elemid2 = face->elementOwners[1];
			if (elemid2 == -1)
				continue;
			int elemid1 = face->elementOwners[0];
			if (elemid1 != id)
				elemid2 = elemid1;
			SRelement *elem2 = model.GetElement(elemid2);
			for (int n = 0; n < numNodesTotal; n++)
			{
				SRnode* node = elem2->getNodeOrMidNode(n);
				if (node->checkUnsup())
					return true;
			}
		}
	}
	return false;
}

void SRelement::scaleStress(double stress[], double& svm, double scaleClip)
{
	//scale stress in element so svm matches min of Svmmax in element,scaleClip*svmMax in model
	double target = svmMax;
	double target2 = scaleClip*model.stressMax;
	if (target2 < target)
		target = target2;
	double rat = target / svm;
	svm *= rat;
	for (int c = 0; c < 6; c++)
		stress[c] *= rat;
}

bool SRelement::checkAnyNodeBreakout()
{
	for (int n = 0; n < GetNumNodesTotal(); n++)
	{
		if (getNodeOrMidNode(n)->checkSacr())
			return true;
	}
	return false;
}

bool SRelement::checkForNodalHotSpot()
{
	//check for nodal hot spot: loop node, see if max stress is much higher than minStress at all other nodes
	double highTol = .0;
	double maxNodalStress = 0.0;
	double minNodalStress = BIG;
	double stress[6];
	int nn = GetNumNodesTotal();
	for (int n = 0; n < nn; n++)
	{
		int nid = getNodeOrMidNodeId(n);
		for (int c = 0; c < 6; c++)
			stress[c] = model.post.nodalStress.Get(nid,c);
		double svm = model.math.GetSvm(stress);
		if (svm > maxNodalStress)
			maxNodalStress = svm;
		if (svm < minNodalStress)
			minNodalStress = svm;
	}

	if (minNodalStress < highTol*minNodalStress)
		return false;

	//node at max is hot spot. mark this element, and any elements that own any of it's corners,
	//as sacrificial
	this->sacrificial = 1;
	for (int n = 0; n < nn; n++)
	{
		int nidElMax = getNodeOrMidNodeId(n);
		for (int e = 0; e < model.GetNumElements(); e++)
		{
			SRelement* elem = model.GetElement(e);
			if (elem->sacrificial != 0)
				continue;
			for (int n = 0; n < elem->GetNumNodesTotal(); n++)
			{
				int nid = elem->getNodeOrMidNodeId(n);
				if (nid == nidElMax)
				{
					elem->sacrificial = 1;
					break;
				}
			}
		}
	}
	return true;
}

void SRelement::CalculateHStress(double r, double s, double t, double stress[])
{
	//calculate stress at r,s,t for an element using h-basis functions
	FillMappingNodesBreakoutOnly(true);//needShapeDerivs is true
	FillMapping(r, s, t);

	double dbdx, dbdy, dbdz, dbdr, dbds, dbdt, uel, vel, wel, ex, ey, ez, gamxy, gamxz, gamyz;
	ex = 0.0;
	ey = 0.0;
	ez = 0.0;
	gamxy = 0.0;
	gamxz = 0.0;
	gamyz = 0.0;
	SRvec3 disp;
	for (int n = 0; n < numNodesTotal; n++)
	{
		int nid = getNodeOrMidNodeId(n);
		dbdr = dNdr[n];
		dbds = dNds[n];
		dbdt = dNdt[n];
		XyzDerivatives(dbdr, dbds, dbdt, dbdx, dbdy, dbdz);
		disp = model.post.nodeDisps.Get(nid);
		uel = disp.d[0];
		vel = disp.d[1];
		wel = disp.d[2];
		ex += dbdx*uel;
		ey += dbdy*vel;
		ez += dbdz*wel;
		gamxy += (dbdy*uel + dbdx*vel);
		gamxz += (dbdz*uel + dbdx*wel);
		gamyz += (dbdz*vel + dbdy*wel);
	}
	double strain[6];
	strain[0] = ex;
	strain[1] = ey;
	strain[2] = ez;
	strain[3] = gamxy;
	strain[4] = gamxz;
	strain[5] = gamyz;
	StraintoStress(r, s, t, strain, stress);
}

double SRelement::CalculateSvmHCentroid()
{
	//calculate raw svm in an element at centroid using h-basis functions (shape functions)
	double stress[6];
	double r, s, t;
	model.map.ElementCentroid(this, r, s, t);
	CalculateHStress(r, s, t, stress);
	double svm = model.math.GetSvm(stress);
	return svm;
}













