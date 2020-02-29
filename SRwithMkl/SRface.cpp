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
// SRface.cpp: implementation of the SRface class.
//
//////////////////////////////////////////////////////////////////////

#include "SRmodel.h"

extern SRmodel model;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

//////////////////////////////////////////////////////////////////////////
// SRfaceUtil static functions are used during creation of global faces //
//////////////////////////////////////////////////////////////////////////

int SRfaceUtil::GetLocalEdge(int n1, int n2)
{
	//find the local edge number of the face corresponding to local node numbers n1, n2
	//return:
		// local edge number 

	if (n1 == 0)
	{
		if (n2 == 1) //1st edge of tri or quad
			return 0;
		else if (n2 == 2) //3rd edge of tri
			return 2;
		else if (n2 == 3) //4th edge of quad
			return 3;
	}
	else if (n1 == 1)
	{
		if (n2 == 0) //1st edge of tri or quad
			return 0;
		else if (n2 == 2) //2nd edge of tri or quad
			return 1;
	}
	else if (n1 == 2)
	{
		if (n2 == 1) //2nd edge of tri or quad
			return 1;
		else if (n2 == 0)// 3rd edge of tri
			return 2;
		else if (n2 == 3) // 3rd edge of quad
			return 2;
	}
	else if (n1 == 3)
	{
		if (n2 == 2) // 3rd edge of quad
			return 2;
		else if (n2 == 0) // 4th edge of quad
			return 3;
	}
	return -1;
}

bool SRfaceUtil::FaceMatch(int n1, int n2, int n3, int n4, int g1, int g2, int g3, int g4, int &gn1, int &gn2, int &gn3, int &gn4)
{
	//see if face with corner nodes n1,n2,n3,(n4) matches 
	//face with corner nodes g1,g2,g3,(g4)
	//input:
		//n1,n2,n3,(n4) = nodes at corners of face (n4 = -1 for triangle)
		//g1,g2,g3,(g4) = nodes at corners of existing global face (g4 = -1 for triangle)
	//output:
		//gn1,gn2,gn3,(gn4) = node order on face g1,g2,g3,(g4) (gn4 = -1 for triangle)
	//return:
		//true if match else false

	if (n4 == -1)
	{
		//tri face
		gn4 = -1;
		if (n1 == g1)
		{
			gn1 = 0;
			if (n2 == g2 && n3 == g3)
			{
				gn2 = 1;
				gn3 = 2;
				return true;
			}
			if (n2 == g3 && n3 == g2)
			{
				gn2 = 2;
				gn3 = 1;
				return true;
			}
		}
		else if (n1 == g2)
		{
			gn1 = 1;
			if (n2 == g1 && n3 == g3)
			{
				gn2 = 0;
				gn3 = 2;
				return true;
			}
			if (n2 == g3 && n3 == g1)
			{
				gn2 = 2;
				gn3 = 0;
				return true;
			}
		}
		else if (n1 == g3)
		{
			gn1 = 2;
			if (n2 == g2 && n3 == g1)
			{
				gn2 = 1;
				gn3 = 0;
				return true;
			}
			if (n2 == g1 && n3 == g2)
			{
				gn2 = 0;
				gn3 = 1;
				return true;
			}
		}
	}
	else
		return QuadFaceMatch(n1, n2, n3, n4, g1, g2, g3, g4, gn1, gn2, gn3, gn4);
	return false;
}

bool SRfaceUtil::FaceMatch(int nv[], int gv[], int gnv[])
{
	//overload with nodes in a vector
	//see if face with corner nodes in vector nv matches 
	//face with corner nodes in vector gv
	//input:
		//nv = vector of corner nodes
		//gv = vector of nodes at corners of existing global face (g4 = -1 for triangle)
	//output:
		//gnv = node order on face g1,g2,g3,(g4)
	//return:
		//true if match else false

	return FaceMatch(nv[0], nv[1], nv[2], nv[3], gv[0], gv[1], gv[2], gv[3], gnv[0], gnv[1], gnv[2], gnv[3]);
}

int SRfaceUtil::GlobalFaceMatch(int &gn1, int &gn2, int &gn3, int &gn4, int n1, int n2, int n3, int n4)
{
	//find which global face has nodes n1,n2,n3,(n4). n4 is -1 for tri face
	//input:
		//n1,n2,n3,(n4) = nodes at corners of face
	//output:
		//gn1,gn2,gn3,(gn4) = node order on corresponding global face
	//return:
		//face id

	gn4 = -1;
	int g1, g2, g3, g4 = -1;
	int nn1face = model.input.GetNumNodeFaces(n1);
	int f = -1;
	if (nn1face == -1)
	{
		//global node n1 overflowed its nodeFaces array, need to search all model faces:
		for (f = 0; f < model.GetNumFaces(); f++)
		{
			SRface* face = model.GetFace(f);
			g1 = face->nodeIds[0];
			g2 = face->nodeIds[1];
			g3 = face->nodeIds[2];
			if (n4 != -1)
			{
				if (face->GetNumNodes() == 3)
					continue;
				else
					g4 = face->nodeIds[3];
			}
			else
			{
				if (face->GetNumNodes() == 4)
					continue;
			}
			if (FaceMatch(n1, n2, n3, n4, g1, g2, g3, g4, gn1, gn2, gn3, gn4))
				return f;
		}
	}
	else
	{
		for (int f = 0; f < nn1face; f++)
		{
			int n1face = model.input.GetNodeFace(n1, f);
			SRface* face = model.GetFace(n1face);
			g1 = face->nodeIds[0];
			g2 = face->nodeIds[1];
			g3 = face->nodeIds[2];
			if (n4 != -1)
			{
				if (face->GetNumNodes() == 3)
					continue;
				else
					g4 = face->nodeIds[3];
			}
			else
			{
				if (face->GetNumNodes() == 4)
					continue;
			}
			if (FaceMatch(n1, n2, n3, n4, g1, g2, g3, g4, gn1, gn2, gn3, gn4))
				return n1face;
		}
	}

	return -1;
}


void SRfaceUtil::GetLocalEdgeNodes(int nej, int lej, int &n1, int &n2)
{
	//get local node numbers at ends of a local edge
	//input:
		//nej = number of edges on face
		//lej = local edge number
	//output:
		//n1,n2 = local node numbers on the face at ends 1 and 2 of the edge

	if (nej == 3)
	{
		if (lej == 0)
		{
			n1 = 0;
			n2 = 1;
		}
		else if (lej == 1)
		{
			n1 = 1;
			n2 = 2;
		}
		else if (lej == 2)
		{
			n1 = 0;
			n2 = 2;
		}
	}
	else
		QuadLocalEdgeNodes(lej, n1, n2);
}

///////////////////////////////////////////////////////////////
//               SRface member functions                     //
///////////////////////////////////////////////////////////////


SRface::SRface()
{
	constraintId = -1;
	flat = -1;
	saveForBreakout = false;
	multifaceForceGroupId = -1;
	unsupported = false;
	onPartBdry = false;
};


void SRface::Create(int idt, int n1, int n2, int n3, int n4)
{
	//create global face
	//input:
		//n1,n2,n3,n4 = node numbers of corners of face; n4=-1 for triangle
	id = idt;
	nodeIds[0] = n1;
	nodeIds[1] = n2;
	nodeIds[2] = n3;
	nodeIds[3] = n4;
	int ncorner = 4;
	if (n4 == -1)
		ncorner = 3;
	localEdges.Allocate(ncorner);
	elementOwners[0] = -1;
	elementOwners[1] = -1;
    tractionJump = 0.0;
	for (int i = 0; i < ncorner; i++)
	{
		int nid = nodeIds[i];
		int n = model.input.GetNumNodeFaces(nid);
		if (n != -1)
		{
			if (n == MAXNODEFACEOWNERS)
				model.input.PutNumNodeFaces(nid, -1);
			else
			{
				model.input.PutNodeFace(nid, n, id);
				n++;
				model.input.PutNumNodeFaces(nid, n);
			}
		}
	}
}

int SRface::GetNodeId(int i)
{
	if (i < 4)
		return nodeIds[i];
	else
		return -1;
}

SRedge * SRface::GetEdge(int localedgenum)
{
	//get global edge corresponding to local edge of this face
	//input:
		//localedgenum = local edge number
	//return:
		//pointer to global edge

	int gej = localEdges.GetPointer(localedgenum)->GetGlobalEdgeId();
	return model.GetEdge(gej);
}

SRnode * SRface::GetNode(int localnodenum)
{
	//get global node corresponding to local node of this face
	//input:
		//localnodenum = local node number
	//return:
		//pointer to global node

	int nodeid = GetNodeId(localnodenum);
	return model.GetNode(nodeid);
}

void SRface::GetForceValue(SRforce *pforce, double rf, double sf, double forceVal[3])
{
	//get force value  at a natural coordinate point on this face
	//input:
		//rf, sf = natural coordinates
	//output:
		//forceVal = 3 dof force vector

	for (int dof = 0; dof < 3; dof++)
		forceVal[dof] = 0.0;

	double N[8];
	int nn;
	model.map.FaceShapeFunctions(this, rf, sf, N);
	nn = 2 * localEdges.GetNum();
	for (int dof = 0; dof < 3; dof++)
		forceVal[dof] = 0.0;

	for (int j = 0; j < nn; j++)
	{
		if (pforce->isPressure())
		{
			double pressure = pforce->GetForceVal(j, 0);
			//pressure is in -outward normal direction:
			SRvec3 norm;
			OutwardNormal(rf, sf, norm);
			for (int dof = 0; dof < 3; dof++)
			{
				double ft = -pressure*norm.d[dof];
				forceVal[dof] += (ft*N[j]);
			}
		}
		else
		{
			for (int dof = 0; dof < 3; dof++)
			{
				double ft = pforce->GetForceVal(j, dof);
				forceVal[dof] += (ft*N[j]);
			}
		}
	}

	if (pforce->GetCoordId() != -1)
	{
		SRcoord* coord = model.GetCoord(pforce->GetCoordId());
		SRvec3 pos, e1l, e2l, e3l, t;
		Position(rf, sf, pos);
		coord->CalculateBasisVectors(pos, e1l, e2l, e3l);

		t.Zero();
		for (int dof = 0; dof < 3; dof++)
			t.d[dof] += (forceVal[0] * e1l.d[dof] + forceVal[1] * e2l.d[dof] + forceVal[2] * e3l.d[dof]);
		t.Copy(forceVal);
	}
}

void SRface::ProcessForce(SRforce* force, SRvec3& ResF)
{
	//fill up contribution of this face's force to global force vector
	//input:
		//mforceId = Id of force in model.forces
		//ResF = resultant force vector on model
	//output:
		//ResF = updated

	//note:
		//adds to global force vector currently stored in model.solution
	double rf, sf, w, detj, bw, forceVal[3];
	int i, dof, nint, gp, eq, gfun;
	double* globalForce = model.GetSolutionVector();
	SRdoubleVector bv;
	bv.Allocate(GetNumGlobalFunctions());
	double* basis = bv.GetVector();
	SRvec3 tmpRes, p;
	nint = model.math.FillGaussPoints(this);
	for (gp = 0; gp < nint; gp++)
	{
		model.math.GetGP2d(gp, rf, sf, w);

		model.basis.FaceBasisFuncs(rf, sf, this, basis);
		detj = Jacobian(rf, sf);
		w *= detj;
		GetForceValue(force, rf, sf, forceVal);

		tmpRes.Assign(forceVal);
		tmpRes.Scale(w);
		ResF.PlusAssign(tmpRes);
		for (i = 0; i < globalFunctionNumbers.GetNum(); i++)
		{
			gfun = globalFunctionNumbers.Get(i);
			bw = basis[i]*w;
			for( dof = 0; dof < 3; dof++)
			{
				eq = model.GetFunctionEquation(gfun,dof);
				if(eq >= 0)
					globalForce[eq] += (bw*forceVal[dof]);
			}
		}
	}
}

void SRface::Clear()
{
	//free memory of an element that is no longer needed
	localEdges.Free();
	globalFunctionNumbers.Free();
}

void SRface::GetLocalEdgeNodes(int lej, int &n1, int &n2)
{
	//get local node numbers at ends of a local edge
	//input:
		//lej = local edge number
	//output:
		//n1,n2 = local node numbers on the face at ends 1 and 2 of the edge

	SRfaceUtil::GetLocalEdgeNodes(localEdges.GetNum(), lej, n1, n2);
}

SRconstraint * SRface::GetConstraint()
{
	//get constraint associated with this face
	//return:
		//pointer to constraint, NULL if none.
	if(constraintId == -1)
		return NULL;
	else
		return model.GetConstraint(constraintId);
}

void SRface::NodeNaturalCoordinates(int lnode, double &rf, double &sf)
{
	//get the natural coordinates of a node of a face
	//input:
		//lnode = corner number
	//output:
		//rf, sf = natural coordinates at that corner
	if (GetNumLocalEdges() == 3)
	{
		rf = model.map.GetTriNoder(lnode);
		sf = model.map.GetTriNodes(lnode);
	}
	else
	{
		rf = model.map.GetQuadNoder(lnode);
		sf = model.map.GetQuadNodes(lnode);
	}
}


void SRface::FillMappingNodes()
{
	//fill in mapping nodes for a face
	//note:
		//fills class variables xnode, ynode, znode, and size

	int i, n, nn, ns;
	ns = GetNumLocalEdges();
	SRnode* node;
	SRnode* node1;
	for (i = 0; i < ns; i++)
	{
		n = GetNodeId(i);
		node = model.GetNode(n);
		xnode[i] = node->GetXyz(0);
		ynode[i] = node->GetXyz(1);
		znode[i] = node->GetXyz(2);
	}
	SRedge* edge;
	nn = 2 * ns;
	int lej;
	for (i = ns; i < nn; i++)
	{
		lej = i - ns;
		edge = GetEdge(lej);
		n = edge->GetMidNodeId();
		node = model.GetNode(n);
		xnode[i] = node->GetXyz(0);
		ynode[i] = node->GetXyz(1);
		znode[i] = node->GetXyz(2);
	}
	//fill map face size; approximate as average of edge lengths:
	size = 0.0;
	for (i = 0; i < ns; i++)
	{
		edge = GetEdge(i);
		node = edge->GetNode(0);
		node1 = edge->GetNode(1);
		size += model.math.Distance(node->Position(), node1->Position());
	}
	size /= ns;
}

double SRface::UnitTriad(double rf, double sf, SRvec3 &p, SRvec3 &e1, SRvec3 &e2, SRvec3 &e3, bool detjonly)
{
	//calculate unit triad at natural point rf,sf on face
	//input:
		//rf,sf = natural coordinates on face
		//detjonly = true if only determinant of jacobian is needed else false
	//output (if not detjonly):
		//p = coordinates at rf,sf
		//e1,e2 = tangent vectors to face at rf,sf
		//e3 = normal vector to face at rf,sf
	//return:
		//determinant of face jacobian

	double N[8], detj;
	int i, nn, ns;
	ns = GetNumNodes();
	nn = 2 * ns;
	if (ns == 3)
	{
		double l1, l2, l3;
		model.map.TriangleAreaCoords(rf, sf, l1, l2, l3);
		double l14 = 4.0*l1 - 1.0;
		double l24 = 4.0*l2 - 1.0;
		double l34 = 4.0*l3 - 1.0;
		if (!detjonly)
		{
			N[0] = l1*(2.0*l1 - 1.0);
			N[1] = l2*(2.0*l2 - 1.0);
			N[2] = l3*(2.0*l3 - 1.0);
			N[3] = 4.0*l1*l2;
			N[4] = 4.0*l2*l3;
			N[5] = 4.0*l1*l3;
		}
		dNdrf[0] = DL1DR*l14;
		dNdsf[0] = DL1DS*l14;
		dNdrf[1] = DL2DR*l24;
		dNdsf[1] = DL2DS*l24;
		dNdrf[2] = 0.0; //DL3DR*l34 but DL3DR=0
		dNdsf[2] = DL3DS*l34;
		dNdrf[3] = 4.0*(DL1DR*l2 + DL2DR*l1);
		dNdsf[3] = 4.0*(DL1DS*l2 + DL2DS*l1);
		dNdrf[4] = 4.0*DL2DR*l3; //4.0*(DL2DR*l3+DL3DR*l2) but DL3DR=0
		dNdsf[4] = 4.0*(DL2DS*l3 + DL3DS*l2);
		dNdrf[5] = 4.0*DL1DR*l3; // 4.0*(DL1DR*l3+DL3DR*l1) but DL3DR=0
		dNdsf[5] = 4.0*(DL1DS*l3 + DL3DS*l1);
	}
	else
	{
		QuadShapeDerivs(rf, sf);
		if (!detjonly)
			model.map.QuadShapeFunctions(rf, sf, N);
	}

	//e1 = dXdrf, normalized
	//e2 = dXdsf
	//e3 = e1 cross e2, normalized,
	//then e2 = e3 cross e1 to find other tangent that is normal to e1.
	e1.Zero();
	e2.Zero();
	if (!detjonly)
	{
		p.Zero();
		for (i = 0; i < nn; i++)
		{
			p.d[0] += xnode[i] * N[i];
			p.d[1] += ynode[i] * N[i];
			p.d[2] += znode[i] * N[i];
		}
	}
	for (i = 0; i < nn; i++)
	{
		e1.d[0] += xnode[i] * dNdrf[i];
		e1.d[1] += ynode[i] * dNdrf[i];
		e1.d[2] += znode[i] * dNdrf[i];
		e2.d[0] += xnode[i] * dNdsf[i];
		e2.d[1] += ynode[i] * dNdsf[i];
		e2.d[2] += znode[i] * dNdsf[i];
	}
	e1.Cross(e2, e3);
	detj = e3.Normalize();
	if (!detjonly)
	{
		e1.Normalize();
		e3.Cross(e1, e2);
	}
	return detj;
}

bool SRface::isFlat()
{
	//determine if a face is flat
	//return:
		//true if flat else false
	if (flat == -1)
	{
		flat = 1;
		//flat not calculated yet:
		if (GetNumNodes() == 3)
		{
			for (int e = 0; e < 3; e++)
			{
				SRedge* edge = GetEdge(e);
				if (!edge->isStraight())
				{
					flat = 0;
					return false;
				}
			}
			return true;
		}

		double rc, sc;
		model.map.FaceCentroid(this, rc, sc);
		SRvec3 norm;
		OutwardNormal(rc, sc, norm);
		//note: for quad, it's not enough to see if edges are straight. face could still be twisted
		// or edges curved in plane. Test if vector from 1st node to all others is normal
		// to normal at centroid:
		SRvec3 pos0;
		LocalNodePosition(0, pos0);
		int nn = 6;
		if (GetNumNodes() == 4)
			nn = 8;
		SRvec3 pos, v01;
		double dot;
		for (int i = 1; i < nn; i++)
		{
			LocalNodePosition(i, pos);
			pos.Subtract(pos0, v01);
			dot = v01.Dot(norm);
			if (fabs(dot) > 0.00001) //not normal to within 0.00001 degrees
			{
				flat = 0;
				break;
			}
		}
	}
	return (bool) flat;
}

void SRface::LocalNodePosition(int lnode, SRvec3& pos)
{
	//look up the position of a local node (corner or midside
	pos.d[0] = xnode[lnode];
	pos.d[1] = ynode[lnode];
	pos.d[2] = znode[lnode];
}

void SRface::NaturalCoordinatesNearMidedge(int lej, double& r, double& s)
{
	//calculate natural coordinates near the midedge of a local edge of this face
	//input:
		//lej = local edge number
	//output:
		//r, s = natural coordinates
	double re, se;
	naturalCoordinatesFromEdge(lej, 0.0, re, se);
	SRvec3 ejnorm;
	model.map.FaceEdgeNaturalNormal(this, lej, ejnorm);
	r = re - 0.01*ejnorm.d[0]; //this is a point near center of edge, slightly inward
	s = se - 0.01*ejnorm.d[1]; //this is a point near center of edge, slightly inward
}

void SRface::Position(double rf, double sf, SRvec3 &p)
{
	//calculate position at natural point rf,sf on face using quadratic mapping
	//input:
		//face = pointer to face
		//rf,sf = natural coordinates on face
	//output:
		//p = position at rf,sf
	double N[8];
	int i, nn, ns;
	ns = GetNumNodes();
	nn = 2 * ns;
	model.map.FaceShapeFunctions(this, rf, sf, N);
	p.Zero();
	for (i = 0; i < nn; i++)
	{
		p.d[0] += xnode[i] * N[i];
		p.d[1] += ynode[i] * N[i];
		p.d[2] += znode[i] * N[i];
	}
}

void SRface::OutwardNormal(double rf, double sf, SRvec3& norm, bool checkOutwardLocally)
{
	//determine outward normal at natural coordinate rf,sf on face
	//input:
		//rf,sf = natural coords on face
		//checkOutwardLocally = true to check that normal points outward at this position, false to use the "flipnormal" setting
	//output:
		//norm = outward normal at rf,sf

	SRvec3 e1, e2, pos;
	UnitTriad(rf, sf, pos, e1, e2, norm);
	bool flip;
	if (!checkOutwardLocally)
		flip = GetFlipNormal();
	else
	{
		//find a point slightly inside element from the point on face:
		int eid = GetElementOwner(0);
		SRelement* elem = model.GetElement(eid);
		int lface = GetElementLocalFace(0);
		SRvec3 onorm, pinside;
		model.map.ElementFaceNaturalNormal(elem, lface, onorm);
		double r, s, t, rinside, sinside, tinside;
		model.map.ElementNaturalCoordsFromFace(elem, lface, rf, sf, r, s, t);
		rinside = r - 0.01*onorm.d[0];
		sinside = s - 0.01*onorm.d[1];
		tinside = t - 0.01*onorm.d[2];
		elem->Position(rinside, sinside, tinside, pinside);
		pos.Subtract(pinside, onorm);
		onorm.Normalize();
		double dot = norm.Dot(onorm);
		if (dot < 0.0)
			flip = true;
		else
			flip = false;
	}
	if (flip)
	{
		//reverse sign so normal points outward from element:
		for (int j = 0; j < 3; j++)
			norm.d[j] = -norm.d[j];
	}
}

double SRface::Jacobian(double rf, double sf)
{
	//calculate the shape functions for a face at a natural coordinate point
	//input:
		//rf, sf = natural coordinates
	//return:
		//determinant of jacobian
	SRvec3 p, e1, e2, e3;
	return UnitTriad(rf, sf, p, e1, e2, e3, true);
}

bool SRface::natCoordsAtCorner(int cornerId, double& r, double& s)
{
	//calculate the natural coordinates of a face a one of it's corners
	//input:
		//node id at corner
	//output:
		//r,s = natural coordinates
	//return
		//true if pos touches one of the corners else false
	int n = GetNumNodes();
	for (int i = 0; i < n; i++)
	{
		if (nodeIds[i] == cornerId)
		{
			if (n == 3)
			{
				r = model.map.GetTriNoder(i);
				s = model.map.GetTriNodes(i);
			}
			else
			{
				r = model.map.GetQuadNoder(i);
				s = model.map.GetQuadNodes(i);
			}
			return true;
		}
	}
	return false;
}

void SRface::SetSaveForBreakout(int elOwner)
{
	if (elOwner == -1)
	{
		saveForBreakout = false;
		return;
	}
	saveForBreakout = true;
	if (elOwner == 1)
	{
		elementOwners[0] = elementOwners[1];
		elementLocalFace[0] = elementLocalFace[1];
		elementOwners[1] = -1;
		elementLocalFace[1] = -1;
	}
};

int SRface::midNodeMatch(int mid)
{
	for (int i = 0; i < GetNumLocalEdges(); i++)
	{
		if (GetEdge(i)->GetMidNodeId() == mid)
			return i + GetNumNodes();
	}
	return -1;
}

double SRface::normalRotation()
{
	if (!IsBoundaryFace())
		return 0.0;
	double rf, sf;
	model.map.FaceCentroid(this, rf, sf);
	SRvec3 norm0, norm;
	deformedNormal(rf, sf, norm0, norm);
	double dot;
	dot = fabs(norm.Dot(norm0));
	return (1.0 - dot);
}

void SRface::deformedNormal(double rf, double sf, SRvec3& norm0, SRvec3 &norm)
{
	//calculate deformed normal at a point rf, sf on a face
	//input:
		//rf,sf = natural coordinates on face
	//output:
		//norm0 = undeformed normal vector to face at rf,sf
		//norm = deformed normal vector to face at rf,sf

	int i, nn, ns;
	SRvec3 er0, es0, er, es;
	double unode[8], vnode[8], wnode[8];
	ns = GetNumNodes();
	nn = 2 * ns;
	if (ns == 3)
	{
		double l1, l2, l3;
		model.map.TriangleAreaCoords(rf, sf, l1, l2, l3);
		double l14 = 4.0*l1 - 1.0;
		double l24 = 4.0*l2 - 1.0;
		double l34 = 4.0*l3 - 1.0;
		dNdrf[0] = DL1DR*l14;
		dNdsf[0] = DL1DS*l14;
		dNdrf[1] = DL2DR*l24;
		dNdsf[1] = DL2DS*l24;
		dNdrf[2] = 0.0; //DL3DR*l34 but DL3DR=0
		dNdsf[2] = DL3DS*l34;
		dNdrf[3] = 4.0*(DL1DR*l2 + DL2DR*l1);
		dNdsf[3] = 4.0*(DL1DS*l2 + DL2DS*l1);
		dNdrf[4] = 4.0*DL2DR*l3; //4.0*(DL2DR*l3+DL3DR*l2) but DL3DR=0
		dNdsf[4] = 4.0*(DL2DS*l3 + DL3DS*l2);
		dNdrf[5] = 4.0*DL1DR*l3; // 4.0*(DL1DR*l3+DL3DR*l1) but DL3DR=0
		dNdsf[5] = 4.0*(DL1DS*l3 + DL3DS*l1);
	}
	else
		QuadShapeDerivs(rf, sf);

	//fill displacements at nodes:
	SRvec3 disp;
	int n = 0;
	for (int i = 0; i < ns; i++)
	{
		GetNode(i)->GetDisp(disp);
		unode[n] = disp.d[0];
		vnode[n] = disp.d[1];
		wnode[n] = disp.d[2];
		n++;
	}
	for (int i = 0; i < ns; i++)
	{
		GetEdge(i)->getDisp(0.0, disp);
		unode[n] = disp.d[0];
		vnode[n] = disp.d[1];
		wnode[n] = disp.d[2];
		n++;
	}

	for (i = 0; i < nn; i++)
	{
		er0.d[0] += xnode[i] * dNdrf[i];
		er0.d[1] += ynode[i] * dNdrf[i];
		er0.d[2] += znode[i] * dNdrf[i];
		es0.d[0] += xnode[i] * dNdsf[i];
		es0.d[1] += ynode[i] * dNdsf[i];
		es0.d[2] += znode[i] * dNdsf[i];
		er.d[0] += unode[i] * dNdrf[i];
		er.d[1] += vnode[i] * dNdrf[i];
		er.d[2] += wnode[i] * dNdrf[i];
		es.d[0] += unode[i] * dNdsf[i];
		es.d[1] += vnode[i] * dNdsf[i];
		es.d[2] += wnode[i] * dNdsf[i];
	}
	er0.Cross(es0, norm0);
	norm0.Normalize();
	er.PlusAssign(er0);
	es.PlusAssign(es0);
	er.Cross(es, norm);
	norm.Normalize();
}

void SRface::RotateVec(int coordid, double rf, double sf, double forceVal[])
{
	SRcoord* coord = model.GetCoord(coordid);
	SRvec3 pos, e1l, e2l, e3l, t;
	Position(rf, sf, pos);
	coord->CalculateBasisVectors(pos, e1l, e2l, e3l);

	t.Zero();
	for (int dof = 0; dof < 3; dof++)
		t.d[dof] += (forceVal[0] * e1l.d[dof] + forceVal[1] * e2l.d[dof] + forceVal[2] * e3l.d[dof]);
	t.Copy(forceVal);
}

void SRface::GetSummedForceValue(double rf, double sf, double forceVal[3])
{
	int dof;
	for (dof = 0; dof < 3; dof++)
		forceVal[dof] = 0.0;
	double ft[3];
	for (int i = 0; i < forceIds.GetNum(); i++)
	{
		SRforce* force = model.forces.GetPointer(forceIds.Get(i));
		GetForceValue(force, rf, sf, ft);
		for (dof = 0; dof < 3; dof++)
			forceVal[dof] += ft[dof];
	}
}

int SRface::GetNodeOrMidNodeId(int i)
{
	int nn = GetNumNodes();
	if (i < nn)
		return nodeIds[i];
	else
	{
		int lej = i - nn;
		return GetEdge(lej)->midnodeId;
	}
}

SRnode* SRface::GetNodeOrMidnode(int localnodenum)
{
	int nid = GetNodeOrMidNodeId(localnodenum);
	return model.GetNode(nid);
}

void SRface::naturalCoordinatesFromEdge(int lej, double rej, double &rf, double &sf, bool useDirection)
{
	//get natural coordinates on face from natural coordinate on local edge 
	//input:
	//lej = local edge number
	//rej = natural coordinate on edge
	//useDirection = true to use the global edge sign else fals
	//output:
	//rf,sf = natural coordinates on face


	double rejt = rej;
	if (useDirection)
		rejt *= localEdges.Get(lej).GetDirection();
	if (localEdges.GetNum() == 3)
	{
		if (lej == 0)
		{
			rf = rejt;
			sf = 0.0;
		}
		else if (lej == 1)
		{
			sf = SQRT3OVER2*(1.0 + rejt);
			rf = 1.0 - sf*SQRT3OVER3;
		}
		else if (lej == 2)
		{
			sf = SQRT3OVER2*(1.0 + rejt);
			rf = SQRT3OVER3*sf - 1.0;
		}
	}
	else
		QuadFaceNaturalCoordinatesFromEdge(lej, rejt, rf, sf);

}

void SRface::OutwardNormalLinear(double rf, double sf, SRvec3& norm)
{
	//determine outward normal at natural coordinate rf,sf on face
	//input:
	//rf,sf = natural coords on face
	//output:
	//norm = outward normal at rf,sf

	SRvec3 e1, e2, pos;
	UnitTriadLinear(rf, sf, pos, e1, e2, norm);
	if (flipNormal)
	{
		//reverse sign so normal points outward from element:
		for (int j = 0; j < 3; j++)
			norm.d[j] = -norm.d[j];
	}
}

void SRface::SetFlipNormal()
{
	//flip the normal to a face if necessary, so it points outward from the 1st element that owns the face
	//note:
	//sets class variable flipNormal
	flipNormal = false;

	//position at centroid of face:
	SRvec3 pos, elcentroid, facenorm, e1, e2, velface;
	double rc, sc;
	model.map.FaceCentroid(this, rc, sc);
	UnitTriadLinear(rc, sc, pos, e1, e2, facenorm);

	int eid = elementOwners[0];
	SRelement* elem = model.GetElement(eid);
	double r, s, t;
	model.map.ElementCentroid(elem, r, s, t);
	elem->PositionLinear(r, s, t, elcentroid);
	pos.Subtract(elcentroid, velface);
	velface.Normalize();
	double dot = facenorm.Dot(velface);
	if (dot < 0.0)
		flipNormal = true;
	else
		flipNormal = false;
}

double SRface::UnitTriadLinear(double rf, double sf, SRvec3 &p, SRvec3 &e1, SRvec3 &e2, SRvec3 &e3)
{
	//calculate unit triad at natural point rf,sf on face using linear mapping
	//input:
	//rf,sf = natural coordinates on face
	//output:
	//p = coordinates at rf,sf
	//e1,e2 = tangent vectors to face at rf,sf
	//e3 = normal vector to face at rf,sf
	//return:
	//determinant of face jacobian

	double N[4], detj;
	int i, nn;
	nn = GetNumNodes();
	if (nn == 3)
	{
		double l1, l2, l3;
		model.map.TriangleAreaCoords(rf, sf, l1, l2, l3);
		N[0] = l1;
		N[1] = l2;
		N[2] = l3;
		dNdrf[0] = DL1DR;
		dNdsf[0] = DL1DS;
		dNdrf[1] = DL2DR;
		dNdsf[1] = DL2DS;
		dNdrf[2] = 0.0; //DL3DR=0
		dNdsf[2] = DL3DS;
	}
	else
	{
		QuadLinearShapeDerivs(rf, sf);
		model.map.QuadLinearShapeFunctions(rf, sf, N);
	}

	//e1 = dXdrf, normalized
	//e2 = dXdsf
	//e3 = e1 cross e2, normalized,
	//then e2 = e3 cross e1 to find other tangent that is normal to e1.
	e1.Zero();
	e2.Zero();
	p.Zero();
	for (i = 0; i < nn; i++)
	{
		p.d[0] += xnode[i] * N[i];
		p.d[1] += ynode[i] * N[i];
		p.d[2] += znode[i] * N[i];
	}
	for (i = 0; i < nn; i++)
	{
		e1.d[0] += xnode[i] * dNdrf[i];
		e1.d[1] += ynode[i] * dNdrf[i];
		e1.d[2] += znode[i] * dNdrf[i];
		e2.d[0] += xnode[i] * dNdsf[i];
		e2.d[1] += ynode[i] * dNdsf[i];
		e2.d[2] += znode[i] * dNdsf[i];
	}
	e1.Cross(e2, e3);
	detj = e3.Normalize();
	e1.Normalize();
	e3.Cross(e1, e2);
	return detj;
}

void SRface::QuadLinearShapeDerivs(double rf, double sf)
{
	//support routine for FaceUnitTriad for quad faces.
	//calculate dNdrf, dNdsf at rf,sf
	//input:
	//rf, sf = natural coordinates
	//note:
	//fills up class variables dNdrf, dNdsf

	int i;
	double r0, s0, ri, si;
	//corners:
	for (i = 0; i < 4; i++)
	{
		ri = model.map.quadnoder[i];
		si = model.map.quadnodes[i];
		r0 = rf*ri;
		s0 = sf*si;
		dNdrf[i] = 0.25*ri*(1.0 + s0);
		dNdsf[i] = 0.25*si*(1.0 + r0);
	}
}

void SRface::NaturalCoordinatesNearNode(int lnode, double &rf, double &sf)
{
	double r, s, rc, sc;
	model.map.FaceCentroid(this, rc, sc);
	NodeNaturalCoordinates(lnode, r, s);
	rf = 0.01*rc + 0.99*r;
	sf = 0.01*sc + 0.99*s;
}

bool SRface::checkPartBdry(bool cornersOnly)
{
	onPartBdry = true;
	int nn = GetNumNodesTotal();
	if (cornersOnly)
		nn = GetNumNodes();

	for (int i = 0; i < nn; i++)
	{
		SRnode* node = GetNodeOrMidnode(i);
		if (node->sacrificialType != onBsurf)
		{
			onPartBdry = false;
			break;
		}
	}
	return onPartBdry;
}



