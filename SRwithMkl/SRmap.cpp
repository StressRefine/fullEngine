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
// SRmap.cpp: implementation of the SRmap class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SRmodel.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

extern SRmodel model;

SRmap::SRmap()
{
	setupWasCalled = false;
	//r,s at nodes of tri. faces:
	//corners:
	trinoder[0] = -1.0;
	trinoder[1] = 1.0;
	trinoder[2] = 0.0;
	trinodes[0] = 0.0;
	trinodes[1] = 0.0;
	trinodes[2] = SQRT3;
	//midsides:
	trinoder[3] = 0.0;
	trinoder[4] = 0.5;
	trinoder[5] = -0.5;
	trinodes[3] = 0.0;
	trinodes[4] = SQRT3OVER2;
	trinodes[5] = SQRT3OVER2;

	//r,s,t at tet corners:
	tetCorners[0].Assign(-1.0, 0.0, 0.0);
	tetCorners[1].Assign(1.0, 0.0, 0.0);
	tetCorners[2].Assign(0.0, SQRT3, 0.0);
	tetCorners[3].Assign(0.0, SQRT3OVER3, TWOSQRTTWOTHIRDS);

	BrickWedgeQuadInit();

	//tet L derivatives (constant with r,s,t):
	dLtdr[0] = DLT1DR;
	dLtdr[1] = DLT2DR;
	dLtdr[2] = DLT3DR;
	dLtdr[3] = DLT4DR;
	dLtds[0] = DLT1DS;
	dLtds[1] = DLT2DS;
	dLtds[2] = DLT3DS;
	dLtds[3] = DLT4DS;
	dLtdt[0] = DLT1DT;
	dLtdt[1] = DLT2DT;
	dLtdt[2] = DLT3DT;
	dLtdt[3] = DLT4DT;

	//mapping from local node numbers on faces to local node numbers
	//in tet element:
	int tetfacelnodetmp0[3] = { 1, 2, 3 };
	int tetfacelnodetmp1[3] = { 0, 2, 3 };
	int tetfacelnodetmp2[3] = { 0, 1, 3 };
	int tetfacelnodetmp3[3] = { 0, 1, 2 };
	for (int i = 0; i < 3; i++)
	{
		tetFaceLocalNodes[0][i] = tetfacelnodetmp0[i];
		tetFaceLocalNodes[1][i] = tetfacelnodetmp1[i];
		tetFaceLocalNodes[2][i] = tetfacelnodetmp2[i];
		tetFaceLocalNodes[3][i] = tetfacelnodetmp3[i];
	}

	//local node numbers at ends of quad, tri, and tet edges:
	quadEdgeLocalNodes[0][0] = 0;
	quadEdgeLocalNodes[0][1] = 1;
	quadEdgeLocalNodes[1][0] = 1;
	quadEdgeLocalNodes[1][1] = 2;
	quadEdgeLocalNodes[2][0] = 3;
	quadEdgeLocalNodes[2][1] = 2;
	quadEdgeLocalNodes[3][0] = 0;
	quadEdgeLocalNodes[3][1] = 3;

	triEdgeLocalNodes[0][0] = 0;
	triEdgeLocalNodes[0][1] = 1;
	triEdgeLocalNodes[1][0] = 1;
	triEdgeLocalNodes[1][1] = 2;
	triEdgeLocalNodes[2][0] = 0;
	triEdgeLocalNodes[2][1] = 2;

	tetEdgeLocalNodes[0][0] = 0;
	tetEdgeLocalNodes[0][1] = 1;
	tetEdgeLocalNodes[1][0] = 1;
	tetEdgeLocalNodes[1][1] = 2;
	tetEdgeLocalNodes[2][0] = 0;
	tetEdgeLocalNodes[2][1] = 2;
	tetEdgeLocalNodes[3][0] = 0;
	tetEdgeLocalNodes[3][1] = 3;
	tetEdgeLocalNodes[4][0] = 1;
	tetEdgeLocalNodes[4][1] = 3;
	tetEdgeLocalNodes[5][0] = 2;
	tetEdgeLocalNodes[5][1] = 3;

	int triFaceEdgeNode0[3] = { 0, 1, 3 };
	int triFaceEdgeNode1[3] = { 1, 2, 4 };
	int triFaceEdgeNode2[3] = { 0, 2, 5 };
	int quadFaceEdgeNode0[3] = { 0, 1, 4 };
	int quadFaceEdgeNode1[3] = { 1, 2, 5 };
	int quadFaceEdgeNode2[3] = { 3, 2, 6 };
	int quadFaceEdgeNode3[3] = { 0, 3, 7 };
	for (int lnode = 0; lnode < 3; lnode++)
	{
		triFaceEdgeNodes[0][lnode] = triFaceEdgeNode0[lnode];
		triFaceEdgeNodes[1][lnode] = triFaceEdgeNode1[lnode];
		triFaceEdgeNodes[2][lnode] = triFaceEdgeNode2[lnode];
		quadFaceEdgeNodes[0][lnode] = quadFaceEdgeNode0[lnode];
		quadFaceEdgeNodes[1][lnode] = quadFaceEdgeNode1[lnode];
		quadFaceEdgeNodes[2][lnode] = quadFaceEdgeNode2[lnode];
		quadFaceEdgeNodes[3][lnode] = quadFaceEdgeNode3[lnode];
	}

	//natural normals:
	SRvec3 e1, e2, norm, p1, p2;
	//brick face 0 = 0 1 2 3 (t = -1)
	brickNatNormals[0].Assign(0.0, 0.0, -1.0);
	//brick face 1 = 4 5 6 7 (t = 1)
	brickNatNormals[1].Assign(0.0, 0.0, 1.0);
	//brick face 2 = 0 3 7 4 ( r = -1)
	brickNatNormals[2].Assign(-1.0, 0.0, 0.0);
	//brick face 3 = 1 2 6 5 (r = 1)
	brickNatNormals[3].Assign(1.0, 0.0, 0.0);
	//brick face 4 = 0 1 5 4 (s = -1)
	brickNatNormals[4].Assign(0.0, -1.0, 0.0);
	//brick face 5 = 3 2 6 7 (s = 1)
	brickNatNormals[5].Assign(0.0, 1.0, 0.0);

	//wedge face 0 = 0 1 2 (t = -1)
	wedgeNatNormals[0].Assign(0.0, 0.0, -1.0);
	//wedge face 1 = 3 4 5 (t = 1)
	wedgeNatNormals[1].Assign(0.0, 0.0, 1.0);
	//wedge face 2 = 0 1 4 3 (s = 0)
	wedgeNatNormals[2].Assign(0.0, -1.0, 0.0);
	//wedge face 3 = 1 2 5 4 
	wedgeNatNormals[3].Assign(SQRT3OVER2, 0.5, 0.0);
	//wedge face 4 = 0 2 5 3
	wedgeNatNormals[4].Assign(-SQRT3OVER2, 0.5, 0.0);

	//tet face 0 = 1 2 3
	p1.Copy(tetCorners[2]);
	p2.Copy(tetCorners[1]);
	p1.Subtract(p2, e1);
	p1.Copy(tetCorners[3]);
	p1.Subtract(p2, e2);
	e1.Cross(e2, norm);
	norm.Normalize();
	tetNatNormals[0].Copy(norm);

	//tet face 1 = 0 2 3
	p1.Copy(tetCorners[3]);
	p2.Copy(tetCorners[0]);
	p1.Subtract(p2, e1);
	p1.Copy(tetCorners[2]);
	p1.Subtract(p2, e2);
	e1.Cross(e2, norm);
	norm.Normalize();
	tetNatNormals[1].Copy(norm);

	//tet face 2 = 0 1 3
	p1.Copy(tetCorners[1]);
	//p2.Copy(tetCorners[0]); -was already done above
	p1.Subtract(p2, e1);
	p1.Copy(tetCorners[3]);
	p1.Subtract(p2, e2);
	e1.Cross(e2, norm);
	norm.Normalize();
	tetNatNormals[2].Copy(norm);

	//tet face 3 = 0 1 2 ( t = 0)
	tetNatNormals[3].Assign(0.0, 0.0, -1.0);

	triNatNormals[0].Assign(0.0, -1.0, 0.0);
	triNatNormals[1].Assign(SQRT3OVER2, 0.5, 0.0);
	triNatNormals[2].Assign(-SQRT3OVER2, 0.5, 0.0);

	quadNatNormals[0].Assign(0.0, -1.0, 0.0);
	quadNatNormals[1].Assign(1.0, 0.0, 0.0);
	quadNatNormals[2].Assign(0.0, 1.0, 0.0);
	quadNatNormals[3].Assign(-1.0, 0.0, 0.0);
}

void SRmap::Setup()
{
	//miscellaneous mapping set up:
		//fill the mapping nodes of all edges, faces, and elements
		//set the flipnormal flag for all faces, so the face normal points
		//outward to 1st element that owns it

	for (int e = 0; e < model.GetNumEdges(); e++)
		model.GetEdge(e)->FillMappingNodes();

	for (int f = 0; f < model.GetNumFaces(); f++)
		model.GetFace(f)->FillMappingNodes();

	for (int e = 0; e < model.GetNumElements(); e++)
		model.GetElement(e)->FillMappingNodes();

	// flipnormal flag all faces and flat status for all faces:
	model.allFacesFlat = true;
	for (int f = 0; f < model.GetNumFaces(); f++)
	{
		SRface* face = model.GetFace(f);
		face->SetFlipNormal();
		if (!face->isFlat())
			model.allFacesFlat = false;
	}
}

double SRmap::BrickTetMap(double rb, double sb, double tb, double &r, double &s, double &t)
{
	//map natural coordinates in a degenerated brick to tet. see notes p 21
	//input:
		//rb,sb,tb = natural coordinates in brick
	//output:
		//r,s,t = natural coordinates in tet

	//fill up brick linear shape functions:
	double tbm = 1.0 - tb;
	double tbp = 1.0 + tb;
	double sbm = 1.0 - sb;
	double sbp = 1.0 + sb;
	r = 0.25*rb*sbm*tbm;
	s = SQRT3OVER2*(ONETHIRD*tbp + 0.5*sbp*tbm);
	t = 0.816496580927726*tbp; //sqrt(2/3)
	double jac = 0.08838834764831844*sbm*tbm*tbm; //0.088... = sqrt(2)/16
	return jac;
}

double SRmap::QuadTriMap(double rq, double sq, double &r, double &s)
{
	//map natural coordinates in a degenerated quad to triangle. see notes p 20
	//input:
		//rq,sq=natural coordinates in quad
	//output:
		//r,s=natural coordinates in triangle
	//return:
		//Jacobian of mapping

    double sqm = 1.0 - sq;
	r = 0.5*rq*sqm;
	s = SQRT3OVER2*(1.0 + sq);
	return 0.5*SQRT3OVER2*sqm;
}

double SRmap::TetMidEdgeMappingDeriv(int n0, int n1, double dldra[], double L[])
{
	//calculate derivative of quadratic mapping on an edge, see notes p 17
	//input:
		//n0, n1 = local node numbers at the ends of the edge for the tet
		//dldra = derivative of tet volume functions with respect to "ra" = appropriate
			//coordinate for the edge, see notes
		//L = tet volume coordinate function values at a natural coordinate in the tet
	//return:
		//derivative of quadratic mapping with respect to "ra"

	double dndr = 4.0*(dldra[n0] * L[n1] + L[n0] * dldra[n1]);
	return dndr;
}

void SRmap::TriangleAreaCoords(double r, double s, double &l1, double &l2, double &l3)
{
	//calculate area coordinates of a triangle at natural coordinate r,s
	//output:
		//l1, l2, l3 = area coordinates;

	l1 = 0.5*(1.0 - r - SQRT3OVER3*s);
	l2 = l1 + r;
	l3 = SQRT3OVER3*s;
}

void SRmap::TriangleAreaCoords(double r, double s, double L[])
{
	//calculate area coordinates of a triangle at natural coordinate r,s
	//output:
		//L[0], L[1], L[2] = area coordinates;

	double l1, l2, l3;
	TriangleAreaCoords(r, s, l1, l2, l3);
	L[0] = l1;
	L[1] = l2;
	L[2] = l3;
}

void SRmap::TetVolumeCoords(double r, double s, double t, double &l1, double &l2, double &l3, double &l4)
{
	//calculate volume coordinates of a tet at natural coordinate r,s,t
	//output:
		//l1, l2, l3, l4 = volume coordinates;
	l1 = 0.5*(1.0 - r - SQRT3OVER3*s - SQRT6OVER6*t);
	l2 = l1 + r;
	l3 = SQRT3OVER3*(s - SQRT8OVER8*t);
	l4 = SQRT3OVERSQRT8*t;
}


void SRmap::TetVolumeCoords(double r, double s, double t, double L[])
{
	//calculate volume coordinates of a tet at natural coordinate r,s,t
	//output:
		//L[0], L[1], L[2], L[3] = tet coordinates;

	double l1, l2, l3, l4;
	TetVolumeCoords(r, s, t, l1, l2, l3, l4);
	L[0] = l1;
	L[1] = l2;
	L[2] = l3;
	L[3] = l4;
}

double SRmap::EdgeArcLength(SRedge* edge,double r)
{
	//determine the arc length for an edge at a natural coordinate
	//input:
		//edge = pointer to edge
		//r = natural coordinate
	//return:
		//arc length

	SRvec3 er;
	return EdgeTangent(edge,r,er);
}

void SRmap::ElementNaturalCoordsAtMidedge(SRelement* elem, int localEdge, double& r, double& s, double& t)
{
	//determine element natural coordinate at mid-edge of a local edge
	//input:
		//elem = pointer to element
		//localEdge = local edge number
	//output:
		//r,s,t = natural coordinates in element
	SRvec3 tmpCorners[2];

	//storage of edge corners:
	int ln0, ln1;
	if (elem->GetType() == tet)
	{
		ln0 = tetEdgeLocalNodes[localEdge][0];
		ln1 = tetEdgeLocalNodes[localEdge][1];
		tmpCorners[0].Copy(tetCorners[ln0]);
		tmpCorners[1].Copy(tetCorners[ln1]);
	}
	else if (elem->GetType() == wedge)
	{
		ln0 = wedgeEdgeLocalNodes[localEdge][0];
		ln1 = wedgeEdgeLocalNodes[localEdge][1];
		tmpCorners[0].Copy(wedgeCorners[ln0]);
		tmpCorners[1].Copy(wedgeCorners[ln1]);
	}
	else
	{
		ln0 = brickEdgeLocalNodes[localEdge][0];
		ln1 = brickEdgeLocalNodes[localEdge][1];
		tmpCorners[0].Copy(brickNodes[ln0]);
		tmpCorners[1].Copy(brickNodes[ln1]);
	}

	//linearly interpolate between r,s,t of the two ends
	//of the local edge. 
	r = 0.5*(tmpCorners[0].d[0] + tmpCorners[1].d[0]);
	s = 0.5*(tmpCorners[0].d[1] + tmpCorners[1].d[1]);
	t = 0.5*(tmpCorners[0].d[2] + tmpCorners[1].d[2]);
}

void SRmap::ElementNaturalCoordsFromEdge(SRelement* elem, int localEdge, double re, double& r, double& s, double& t)
{
	//determine element natural coordinate at a natural coordinate on a local edge
	//input:
		//elem = pointer to element
		//localEdge = local edge number
		//re = natural coordinate on localEdge
	//output:
		//r,s,t = natural coordinates in element
	SRvec3 tmpCorners[2];

	//storage of edge corners:
	int ln0, ln1;
	if (elem->GetType() == tet)
	{
		ln0 = tetEdgeLocalNodes[localEdge][0];
		ln1 = tetEdgeLocalNodes[localEdge][1];
		tmpCorners[0].Copy(tetCorners[ln0]);
		tmpCorners[1].Copy(tetCorners[ln1]);
	}
	else if (elem->GetType() == wedge)
	{
		ln0 = wedgeEdgeLocalNodes[localEdge][0];
		ln1 = wedgeEdgeLocalNodes[localEdge][1];
		tmpCorners[0].Copy(wedgeCorners[ln0]);
		tmpCorners[1].Copy(wedgeCorners[ln1]);
	}
	else
	{
		ln0 = brickEdgeLocalNodes[localEdge][0];
		ln1 = brickEdgeLocalNodes[localEdge][1];
		tmpCorners[0].Copy(brickNodes[ln0]);
		tmpCorners[1].Copy(brickNodes[ln1]);
	}

	//linearly interpolate between r,s,t of the two ends
	//of the local edge. 
	double blend0 = (1.0 - re) / 2.0;
	double blend1 = (1.0 + re) / 2.0;
	r = blend0*tmpCorners[0].d[0] + blend1*tmpCorners[1].d[0];
	s = blend0*tmpCorners[0].d[1] + blend1*tmpCorners[1].d[1];
	t = blend0*tmpCorners[0].d[2] + blend1*tmpCorners[1].d[2];
}



void SRmap::ElementNaturalCoordsFromEdge(SRelement* elem, SRvec3& edgePos, double& r, double& s, double& t)
{
	//determine element natural coordinates from natural coordinate along global edge
	//input:
		//elem = pointer to element
		//edgePos = physical position on edge
		//r, s, t = starting guess for corresponding element natural coordinates. use element centroid if unknown
	//output:
		//r,s,t = natural coordinates in element
	ElementInverseMap(elem, edgePos, r, s, t, 1.e-10);

}

void SRmap::ElementNaturalCoordsFromFace(SRelement* elem, SRvec3& facePos, double& r, double& s, double& t)
{
	//find element natural coordinates corresponding to natural coordinates on a global face
	//input:
		//elem = pointer to element
		//facePos = physical position on face
	//r, s, t = starting guess for corresponding element natural coordinates. use element centroid if unknown
	//output:
		//r, s, t = natural coordinates in element corresponding to facePos
	ElementInverseMap(elem, facePos, r, s, t, 1.e-10);

}

void SRmap::ElementNaturalCoordsFromFace(SRelement* elem, int lface, double rf, double sf, double& r, double& s, double& t)
{
	//find element natural coordinates corresponding to coordinates on a face
	//input:
		//elem = pointer to element
		//lface = local face number on element
		//rf,sf = natural coordinates on global face
	//return:
		//r,s, t = natural coordinates in element

	SRvec3 faceCorners[4];
	double N[4];
	int i;
	int nn = elem->GetFace(lface)->GetNumLocalEdges();
	int gnv[4];
	for (i = 0; i < nn; i++)
	{
		int ln = elem->GetLocalFaceGlobalNodeOrder(lface, i);
		gnv[ln] = i;
	}

	switch (elem->GetType())
	{
	case tet:
		switch (lface)
		{
		case 0:
			faceCorners[gnv[0]] = tetCorners[1];
			faceCorners[gnv[1]] = tetCorners[2];
			faceCorners[gnv[2]] = tetCorners[3];
			break;
		case 1:
			faceCorners[gnv[0]] = tetCorners[0];
			faceCorners[gnv[1]] = tetCorners[2];
			faceCorners[gnv[2]] = tetCorners[3];
			break;
		case 2:
			faceCorners[gnv[0]] = tetCorners[0];
			faceCorners[gnv[1]] = tetCorners[1];
			faceCorners[gnv[2]] = tetCorners[3];
			break;
		case 3:
			faceCorners[gnv[0]] = tetCorners[0];
			faceCorners[gnv[1]] = tetCorners[1];
			faceCorners[gnv[2]] = tetCorners[2];
			break;
		default:
			ERROREXIT;
		}

		break;
	case wedge:
		switch (lface)
		{
		case 0:
			faceCorners[gnv[0]] = wedgeCorners[0];
			faceCorners[gnv[1]] = wedgeCorners[1];
			faceCorners[gnv[2]] = wedgeCorners[2];
			break;
		case 1:
			faceCorners[gnv[0]] = wedgeCorners[3];
			faceCorners[gnv[1]] = wedgeCorners[4];
			faceCorners[gnv[2]] = wedgeCorners[5];
			break;
		case 2:
			faceCorners[gnv[0]] = wedgeCorners[0];
			faceCorners[gnv[1]] = wedgeCorners[1];
			faceCorners[gnv[2]] = wedgeCorners[4];
			faceCorners[gnv[3]] = wedgeCorners[3];
			break;
		case 3:
			faceCorners[gnv[0]] = wedgeCorners[1];
			faceCorners[gnv[1]] = wedgeCorners[2];
			faceCorners[gnv[2]] = wedgeCorners[5];
			faceCorners[gnv[3]] = wedgeCorners[4];
			break;
		case 4:
			faceCorners[gnv[0]] = wedgeCorners[0];
			faceCorners[gnv[1]] = wedgeCorners[2];
			faceCorners[gnv[2]] = wedgeCorners[5];
			faceCorners[gnv[3]] = wedgeCorners[3];
			break;
		default:
			ERROREXIT;
		}
		break;
	case brick:
		switch (lface)
		{
		case 0:
			faceCorners[gnv[0]] = brickNodes[0];
			faceCorners[gnv[1]] = brickNodes[1];
			faceCorners[gnv[2]] = brickNodes[2];
			faceCorners[gnv[3]] = brickNodes[3];
			break;
		case 1:
			faceCorners[gnv[0]] = brickNodes[4];
			faceCorners[gnv[1]] = brickNodes[5];
			faceCorners[gnv[2]] = brickNodes[6];
			faceCorners[gnv[3]] = brickNodes[7];
			break;
		case 2:
			faceCorners[gnv[0]] = brickNodes[0];
			faceCorners[gnv[1]] = brickNodes[3];
			faceCorners[gnv[2]] = brickNodes[7];
			faceCorners[gnv[3]] = brickNodes[4];
			break;
		case 3:
			faceCorners[gnv[0]] = brickNodes[1];
			faceCorners[gnv[1]] = brickNodes[2];
			faceCorners[gnv[2]] = brickNodes[6];
			faceCorners[gnv[3]] = brickNodes[5];
			break;
		case 4:
			faceCorners[gnv[0]] = brickNodes[0];
			faceCorners[gnv[1]] = brickNodes[1];
			faceCorners[gnv[2]] = brickNodes[5];
			faceCorners[gnv[3]] = brickNodes[4];
			break;
		case 5:
			faceCorners[gnv[0]] = brickNodes[3];
			faceCorners[gnv[1]] = brickNodes[2];
			faceCorners[gnv[2]] = brickNodes[6];
			faceCorners[gnv[3]] = brickNodes[7];
			break;
		default:
			ERROREXIT;
		}
		break;
	}

	if (nn == 3)
		TriangleAreaCoords(rf, sf, N);
	else
		model.map.QuadLinearShapeFunctions(rf, sf, N);

	//linearly interpolate r,s,t from the corner values:
	r = 0.0;
	s = 0.0;
	t = 0.0;
	for (i = 0; i < nn; i++)
	{
		r += (N[i] * faceCorners[i].d[0]);
		s += (N[i] * faceCorners[i].d[1]);
		t += (N[i] * faceCorners[i].d[2]);
	}
}



void SRmap::FaceCentroid(SRface* face, double& rc, double& sc)
{
	//calculate the natural coordinates of the centroid of a face
	//input:
		//face = pointer to face
	//output:
		//rc,sc = natural coordinates of centroid

	if(face->GetNumNodes() == 3)
	{
		rc = TRI_CENTROID_R;
		sc = TRI_CENTROID_S;
	}
	else
	{
		rc = 0.0;
		sc = 0.0;
	}
}

void SRmap::FaceCentroid(SRface* face, SRvec3& pos)
{
	//calculate the position of the centroid of face
	//input:
		//face = pointer to face
	//output:
		//pos =  position of centroid

	double rc, sc;
	FaceCentroid(face, rc, sc);
	face->Position(rc, sc, pos);
}


void SRmap::ElementCentroid(SRelement* elem,SRvec3& pos)
{
	//calculate the position of the centroid of an element
	//input:
		//elem = pointer to element
	//output:
		//pos = centroid position

	double r, s, t;
	ElementCentroid(elem, r, s, t);
	elem->Position(r, s, t, pos);
}

void SRmap::ElementCentroid(SRelement* elem,double& r, double& s, double& t)
{
	//calculate the natural coordinates of the centroid of an element
	//input:
		//elem = pointer to element
	//output:
		//r,s,t = natural coordinates
	if (elem->GetType() == tet)
	{
		r = TET_CENTROID_R;
		s = TET_CENTROID_S;
		t = TET_CENTROID_T;
	}
	else if (elem->GetType() == wedge)
	{
		r = TRI_CENTROID_R;
		s = TRI_CENTROID_S;
		t = 0.0;
	}
	else
	{
		r = 0.0;
		s = 0.0;
		t = 0.0;
	}
}

double SRmap::EdgeTangent(SRedge* edge, double r, SRvec3& er,bool normalize)
{
	//calculate the tangent vector to an edge at natural coordinate r
	//input:
		//edge = pointer to edge
		//r = natural coordinates
        //normalize = true to normalize the tangent vector else false
	//output:
        // er = tangent vector
    //return:
        //length of er = ||dX/dr||

	double dNdr[3];
	EdgeShapeDerives(r, dNdr);

	er.Zero();
	for(int i = 0; i < 3; i++)
	{
		er.d[0] += edge->getXnode(i) * dNdr[i];
		er.d[1] += edge->getYnode(i) * dNdr[i];
		er.d[2] += edge->getZnode(i) * dNdr[i];
	}
	if(normalize)
		return er.Normalize();
	else
		return er.Length();
}

void SRmap::FaceShapeFunctions(SRface* face, double rf, double sf, double N[])
{
	//calculate the shape functions for a face at a natural coordinate point
	//input:
		//face = pointer to face
		//rf, sf = natural coordinates
	//output:
		//N = shape functions 

	if(face->GetNumLocalEdges() == 3)
	{
		double l1, l2, l3;
		TriangleAreaCoords(rf, sf, l1, l2, l3);
		N[0] = l1*(2.0*l1 - 1.0);
		N[1] = l2*(2.0*l2 - 1.0);
		N[2] = l3*(2.0*l3 - 1.0);
		N[3] = 4.0*l1*l2;
		N[4] = 4.0*l2*l3;
		N[5] = 4.0*l1*l3;
	}
	else
	{
		QuadShapeFunctions(rf, sf, N);
	}
}

int SRmap::EdgeShapeFunctions(double re, double N[])
{
	//calculate the shape functions for an edge at a natural coordinate
	//input:
		//edge = pointer to edge
		//re = natural coordinate
	//output:
		//N = shape functions

	N[0] = 0.5*re*(re - 1.0);
	N[1] = 0.5*re*(re + 1.0);
	N[2] = 1.0 - re*re;
	return 3;
}

double SRmap::ApproxFaceArea(SRface *face)
{
	//return:
		//approximate area of a face
	double rf, sf;
	FaceCentroid(face, rf, sf);
	double detj = face->Jacobian(rf, sf);
	double area;
	if (face->GetNumLocalEdges() == 3)
	{
		//tri face:
		area = detj * SQRT3;
	}
	else
	{
		//quad face:
		area = detj * 4.0;
	}
	return area;
}

double SRmap::FaceArea(SRface *face)
{
	//return:
		//accurate area of a face
	double area = 0.0;
	int nint = model.math.FillFaceGaussPoints(2, face->GetNumLocalEdges());
	double rf, sf, w;
	for (int gp = 0; gp < nint; gp++)
	{
		model.math.GetGP2d(gp, rf, sf, w);
		double detj = face->Jacobian(rf, sf);
		area += (w*detj);
	}
	return area;
}

bool SRmap::ElementInverseMap(SRelement* elem, SRvec3& x, double& r, double& s, double& t, double tolin)
{
	//invert element mapping to determine r,s,t corresponding to physical point x
	//input:
		//elem=pointer to element
		//x=coordinates of a point in or near element
		//startingGuess=true if starting guess has been assigned to r,s,t else false
	//output:
		//r,s,t=natural coordinates in elem corresponding to x
	//return
		//true if converged else false
	double err, atol = INVERSEMAP_ABSTOL, rtol = INVERSEMAP_RELTOL, tol, xl;
	int nits = 0, nitsmax = INVERSEMAP_MAX_ITS;
	SRvec3 xt, dx, dr;
	SRmat33 jactrans;
	xl = elem->GetSize();
	tol = rtol*xl;
	if (tol < atol)
		tol = atol;
	if (tolin > TINY)
		tol = tolin;

	//Newton's method loop
	while (nits <= nitsmax)
	{
		elem->Position(r, s, t, xt);
		x.Subtract(xt, dx);
		err = dx.Length();
		if (err < tol)
			return true;
		elem->FillJacobian(r, s, t);
		elem->GetElJac()->Transpose(jactrans);
		jactrans.Solve(dx.d, dr.d);
		r += dr.d[0];
		s += dr.d[1];
		t += dr.d[2];
		nits++;
	}
	ERROREXIT;
	return false;
}

int SRmap::GetFaceEdgeLocalNode(int Nedge, int lej, int node)
{
	//look up local node number of a face corresponding to node number at end of local edge
	//input:
		//Nedge = number of edges of the face
		//lej = local edge number
		//node = local node number of the end of the edge
	if (Nedge == 3)
		return triFaceEdgeNodes[lej][node];
	else
		return quadFaceEdgeNodes[lej][node];
}

void SRmap::EdgeShapeDerives(double r, double* dNdr)
{
	//calculate the shape derivatoves for an edge at natural coordinate r
	//input:
		//r = natural coordinates
	//output:
		//dNdr = shape derivatoves
	dNdr[0] = r - 0.5;
	dNdr[1] = r + 0.5;
	dNdr[2] = -2.0*r;
}

void SRmap::ElementFaceNaturalNormal(SRelement* elem, int lface, SRvec3& norm)
{
	//look up the normal to a local face of an element in natural space
	//input:
		//elem = pointer to element
		//lface = local face of element
	//output:
		//norm = normal vector in natural space
	if (elem->GetType() == brick)
		norm.Copy(brickNatNormals[lface]);
	else if (elem->GetType() == tet)
		norm.Copy(tetNatNormals[lface]);
	else if (elem->GetType() == wedge)
		norm.Copy(wedgeNatNormals[lface]);
}

void SRmap::FaceEdgeNaturalNormal(SRface* face, int lej, SRvec3& norm)
{
	//look up the normal to a local edge of an face in natural space
	//input:
		//face = pointer to face
		//lej = local edge of face
	//output:
		//norm = normal vector in natural space
	if (face->GetNumNodes() == 3)
		norm.Copy(triNatNormals[lej]);
	else
		norm.Copy(quadNatNormals[lej]);
}

void SRmap::FaceShapeFunctionsLinear(SRface* face, double rf, double sf, double N[])
{
	//calculate the shape functions for a face at a natural coordinate point
	//input:
	//face = pointer to face
	//rf, sf = natural coordinates
	//output:
	//N = shape functions 

	if (face->GetNumLocalEdges() == 3)
		TriangleAreaCoords(rf, sf, N[0], N[1], N[2]);
	else
		QuadLinearShapeFunctions(rf, sf, N);
}

