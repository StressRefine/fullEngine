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
// SRelementBrickWedge.cpp: implementation of the SRelement class
//	for bricks and wedges
//
//////////////////////////////////////////////////////////////////////


#include "SRmodel.h"

extern SRmodel model;


void SRelement::GetWedgeFaceNodes(int lface, int &n1, int &n2, int &n3, int &n4)
{
	//get the global numbers for a local face of a wedge
	//input:
		//lface = local face number
	//output:
		//n1, n2, n3, n4 = node numbers
        //n4 == -1 for triangular face

	if (lface == 0)
	{
		//local face 1 = 1-2-3
		n1 = nodeIds.Get(0);
		n2 = nodeIds.Get(1);
		n3 = nodeIds.Get(2);
		n4 = -1;
	}
	else if (lface == 1)
	{
		//local face 2 = 4-5-6
		n1 = nodeIds.Get(3);
		n2 = nodeIds.Get(4);
		n3 = nodeIds.Get(5);
		n4 = -1;
	}
	else if (lface == 2)
	{
		//local face 3 = 1-2-5-4
		n1 = nodeIds.Get(0);
		n2 = nodeIds.Get(1);
		n3 = nodeIds.Get(4);
		n4 = nodeIds.Get(3);
	}
	else if (lface == 3)
	{
		//local face 4 = 2-3-6-5
		n1 = nodeIds.Get(1);
		n2 = nodeIds.Get(2);
		n3 = nodeIds.Get(5);
		n4 = nodeIds.Get(4);
	}
	else if (lface == 4)
	{
		//local face 5 = 1-3-6-4
		n1 = nodeIds.Get(0);
		n2 = nodeIds.Get(2);
		n3 = nodeIds.Get(5);
		n4 = nodeIds.Get(3);
	}
}

void SRelement::GetBrickFaceNodes(int lface, int &n1, int &n2, int &n3, int &n4)
{
	//get the global numbers for a local face of a brick
	//input:
		//lface  = local face number
	//output:
		//n1, n2, n3, n4 = node numbers
	if (lface == 0)
	{
		//local face 1 = 1-2-3-4
		n1 = nodeIds.Get(0);
		n2 = nodeIds.Get(1);
		n3 = nodeIds.Get(2);
		n4 = nodeIds.Get(3);
	}
	else if (lface == 1)
	{
		//local face 2 = 5-6-7-8
		n1 = nodeIds.Get(4);
		n2 = nodeIds.Get(5);
		n3 = nodeIds.Get(6);
		n4 = nodeIds.Get(7);
	}
	else if (lface == 2)
	{
		//local face 3 = 1-4-8-5
		n1 = nodeIds.Get(0);
		n2 = nodeIds.Get(3);
		n3 = nodeIds.Get(7);
		n4 = nodeIds.Get(4);
	}
	else if (lface == 3)
	{
		//local face 4 = 2-3-7-6
		n1 = nodeIds.Get(1);
		n2 = nodeIds.Get(2);
		n3 = nodeIds.Get(6);
		n4 = nodeIds.Get(5);
	}
	else if (lface == 4)
	{
		//local face 5 = 1-2-6-5
		n1 = nodeIds.Get(0);
		n2 = nodeIds.Get(1);
		n3 = nodeIds.Get(5);
		n4 = nodeIds.Get(4);
	}
	else if (lface == 5)
	{
		//local face 6 = 4-3-7-8
		n1 = nodeIds.Get(3);
		n2 = nodeIds.Get(2);
		n3 = nodeIds.Get(6);
		n4 = nodeIds.Get(7);
	}
}

void SRelement::GetWedgeNodePos(int lnode, double &r, double &s, double &t)
{
	//get the element coordinates at a node of the element
	//input:
		//lnode = local node number
	//output:
        //r,s,t = element coordinates
	if(lnode < 3)
	{
		r = model.map.GetTriNoder(lnode);
		s = model.map.GetTriNodes(lnode);
		t = -1.0;
	}
	else
	{
		r = model.map.GetTriNoder(lnode - 3);
		s = model.map.GetTriNodes(lnode - 3);
		t = 1.0;
	}
}

void SRelement::GetBrickNodePos(int lnode, double &r, double &s, double &t)
{
	//get the element coordinates at a node of the element
	//input:
		//lnode = local node number
	//output:
        //r,s,t = element coordinates

	SRvec3 rst = model.map.GetBrickNodeCoord(lnode);
	r = rst.d[0];
	s = rst.d[1];
	t = rst.d[2];
}

void SRelement::GetQuadFaceDerivs(int lface, double rfl, double sfl, double &drfdrfl, double &drfdsfl, double &dsfdrfl, double &dsfdsfl)
{
	//get derivatives of global face coords w.r.t. local face coords for a quad face
	//input:
		//lface = local face number
		//rfl,sfl = local coordinates on face
	//output:
		//drfdrfl = derivative of global face coordinate rf w.r.t. local face coordinate rfl
		//drfdrsl = derivative of global face coordinate rf w.r.t. local face coordinate sfl
		//dsfdrfl = derivative of global face coordinate sf w.r.t. local face coordinate rfl
		//dsfdsfl = derivative of global face coordinate sf w.r.t. local face coordinate sfl

	int i, node1, node2, node3, node4;
	double dNdrfl[4], dNdsfl[4], ri, si, r0, s0;
	for (i = 0; i < 4; i++)
	{
		ri = model.map.GetQuadNoder(i);
		si = model.map.GetQuadNodes(i);
		r0 = rfl*ri;
		s0 = sfl*si;
		dNdrfl[i] = 0.25*ri*(1.0 + s0);
		dNdsfl[i] = 0.25*si*(1.0 + r0);
	}
	node1 = localFaces.GetPointer(lface)->globalNodeOrder[0];
	node2 = localFaces.GetPointer(lface)->globalNodeOrder[1];
	node3 = localFaces.GetPointer(lface)->globalNodeOrder[2];
	node4 = localFaces.GetPointer(lface)->globalNodeOrder[3];
	drfdrfl = dNdrfl[node2] + dNdrfl[node3] - dNdrfl[node1] - dNdrfl[node4];
	drfdsfl = dNdsfl[node2] + dNdsfl[node3] - dNdsfl[node1] - dNdsfl[node4];
	dsfdrfl = dNdrfl[node3] + dNdrfl[node4] - dNdrfl[node1] - dNdrfl[node2];
	dsfdsfl = dNdsfl[node3] + dNdsfl[node4] - dNdsfl[node1] - dNdsfl[node2];
}

void SRelement::GetQuadFaceGlobalCoords(int lface, double rfl, double sfl, double &rf, double &sf)
{
	//get global face coords for a quad face
	//input:
		//rfl,sfl = local coordinates on face
	//output:
		//rf,sf = global coordinates on face

	int i, node1, node2, node3, node4;
	double Nv[4], ri, si, r0, s0;
	for (i = 0; i < 4; i++)
	{
		ri = model.map.GetQuadNoder(i);
		si = model.map.GetQuadNodes(i);
		r0 = rfl*ri;
		s0 = sfl*si;
		Nv[i] = 0.25*(1.0 + r0)*(1.0 + s0);
	}

	node1 = localFaces.GetPointer(lface)->globalNodeOrder[0];
	node2 = localFaces.GetPointer(lface)->globalNodeOrder[1];
	node3 = localFaces.GetPointer(lface)->globalNodeOrder[2];
	node4 = localFaces.GetPointer(lface)->globalNodeOrder[3];
	rf = Nv[node2] + Nv[node3] - Nv[node1] - Nv[node4];
	sf = Nv[node3] + Nv[node4] - Nv[node1] - Nv[node2];
}

void SRelement::FillWedgeJacobian(double r, double s, double t, bool fillElJac)
{
	//calculate derivatives of mapping and jacobian at a natural coordinate point in a wedge
	//input:
		//r, s, t = natural coordinates
		//fillEljac = true to fill jacobian (eljac), false to just fill mapping derivatives
	//note:
		//fills up class variables dNdr, dNds, dNdt, and elJac

	double L[3], dLdr[3], dLds[3];
	dLdr[0] = DL1DR;
	dLdr[1] = DL2DR;
	dLdr[2] = DL3DR;
	dLds[0] = DL1DS;
	dLds[1] = DL2DS;
	dLds[2] = DL3DS;
	model.map.TriangleAreaCoords(r, s, L);
	double L2Lm1, L4m1, dLdr4Lm1, dLds4Lm1;
	double t2, tm1, tp1;
	t2 = 1.0 - t*t;
	tm1 = 1.0 - t;
	tp1 = 1.0 + t;

	//corner nodes + quad-edge midnodes:
	int itm1, itp1, iQuadMid;
	for (itm1 = 0; itm1 < 3; itm1++)
	{
		itp1 = itm1 + 3;
		iQuadMid = itm1 + 12;
		L2Lm1 = L[itm1] * (2.0*L[itm1] - 1.0);
		L4m1 = 4.0*L[itm1] - 1.0;
		dLdr4Lm1 = dLdr[itm1] * L4m1;//=dLIdr*(4LI-1)
		dLds4Lm1 = dLds[itm1] * L4m1;//=dLIds*(4LI-1)
		dNdr[iQuadMid] = dLdr[itm1] * t2;
		dNds[iQuadMid] = dLds[itm1] * t2;
		dNdt[iQuadMid] = -2.0*t*L[itm1];
		dNdr[itm1] = 0.5*(dLdr4Lm1*tm1 - dNdr[iQuadMid]);
		dNds[itm1] = 0.5*(dLds4Lm1*tm1 - dNds[iQuadMid]);
		dNdt[itm1] = -0.5*(L2Lm1 + dNdt[iQuadMid]);
		dNdr[itp1] = 0.5*(dLdr4Lm1*tp1 - dNdr[iQuadMid]);
		dNds[itp1] = 0.5*(dLds4Lm1*tp1 - dNds[iQuadMid]);
		dNdt[itp1] = 0.5*(L2Lm1 - dNdt[iQuadMid]);
	}
	//tri midedge nodes:
	static int liv[3] = { 0, 1, 0 };
	static int ljv[3] = { 1, 2, 2 };
	double lilj2, tmpr, tmps;
	int I, J;
	for (itm1 = 6; itm1 < 9; itm1++)
	{
		itp1 = itm1 + 3;
		I = liv[itm1 - 6];
		J = ljv[itm1 - 6];
		lilj2 = 2.0*L[I] * L[J];
		tmpr = 2.0*(dLdr[I] * L[J] + dLdr[J] * L[I]);
		tmps = 2.0*(dLds[I] * L[J] + dLds[J] * L[I]);
		dNdr[itm1]= tmpr*tm1;
		dNds[itm1]= tmps*tm1;
		dNdt[itm1]= -lilj2;
		dNdr[itp1]= tmpr*tp1;
		dNds[itp1]= tmps*tp1;
		dNdt[itp1]= lilj2;
	}
	if (fillElJac)
		FillJac();
}

void SRelement::FillBrickJacobian(double r, double s, double t, bool fillElJac)
{
	//calculate  jacobian at a natural coordinate point in a brick
	//input:
		//r, s, t = natural coordinates
		//fillEljac = true to fill jacobian (eljac), false to just fill mapping derivatives
	//note:
		//fills up class variable elJac

	BrickShapeDerivs(r, s, t);
	if (fillElJac)
		FillJac();
}

void SRelement::BrickShapeDerivs(double r, double s, double t)
{
	//calculate derivatives of mapping at a natural coordinate point in a brick
	//see Zienkiewicz p 185
	//input:
		//r, s, t = natural coordinates
	//note:
		//fills up class variable dNdr, dNds, dNdt

	double r0, s0, t0, ri, si, ti;
	double s01, t01, r01;
		int i;

	//derivatives of corner shape functions:
	for (i = 0; i < 8; i++)
	{
		ri = model.map.GetBrickNodeCoord(i).d[0];
		si = model.map.GetBrickNodeCoord(i).d[1];
		ti = model.map.GetBrickNodeCoord(i).d[2];
		r0 = ri*r;
		s0 = si*s;
		t0 = ti*t;
		r01 = 1.0 + r0;
		s01 = 1.0 + s0;
		t01 = 1.0 + t0;
		dNdr[i] = 0.125*ri*s01*t01*(2.0*r0 + s0 + t0 - 1.0);
		dNds[i] = 0.125*si*r01*t01*(2.0*s0 + r0 + t0 - 1.0);
		dNdt[i] = 0.125*ti*r01*s01*(2.0*t0 + r0 + s0 - 1.0);
	}
	double r2 = 1.0 - r*r;
	double s2 = 1.0 - s*s;
	double t2 = 1.0 - t*t;
	//derivatives of midedge shape functions:

	//edges parallel to r:
	for (i = 8; i <= 14; i += 2)
	{
		si = model.map.GetBrickNodeCoord(i).d[1];
		ti = model.map.GetBrickNodeCoord(i).d[2];
		s0 = si*s;
		t0 = ti*t;
		s01 = 1.0 + s0;
		t01 = 1.0 + t0;
		dNdr[i] = -0.5*r*s01*t01;
		dNds[i] = 0.25*si*r2*t01;
		dNdt[i] = 0.25*ti*r2*s01;
	}
	//edges parallel to s:
	for (i = 9; i <= 15; i += 2)
	{
		ri = model.map.GetBrickNodeCoord(i).d[0];
		ti = model.map.GetBrickNodeCoord(i).d[2];
		r0 = ri*r;
		t0 = ti*t;
		r01 = 1.0 + r0;
		t01 = 1.0 + t0;
		dNds[i] = -0.5*s*r01*t01;
		dNdr[i] = 0.25*ri*s2*t01;
		dNdt[i] = 0.25*ti*s2*r01;
	}
	//edges parallel to t:
	for (i = 16; i < 20; i++)
	{
		ri = model.map.GetBrickNodeCoord(i).d[0];
		si = model.map.GetBrickNodeCoord(i).d[1];
		r0 = ri*r;
		s0 = si*s;
		r01 = 1.0 + r0;
		s01 = 1.0 + s0;
		dNdt[i] = -0.5*t*r01*s01;
		dNdr[i] = 0.25*ri*t2*s01;
		dNds[i] = 0.25*si*t2*r01;
	}
}

int SRelement::WedgeShapeFns(double r, double s, double t, double N[])
{
	//calculate  shape functions at a natural coordinate point in a wedge
	//input:
		//r, s, t = natural coordinates
	//output:
		//N = shape functions
	//return:
		//number of functions

	double L[3];
	double tp1 = 1.0 + t;
	double tm1 = 1.0 - t;
	model.map.TriangleAreaCoords(r, s, L);
	double t2, L2Lm1;
	int itm1, itp1, imid;
	t2 = 1.0 - t*t;

	//corner nodes + quad-edge midnodes:
	for (itm1 = 0; itm1 < 3; itm1++)
	{
		itp1 = itm1 + 3;
		imid = itm1 + 12;
		L2Lm1 = L[itm1] * (2.0*L[itm1] - 1.0);
		N[imid] = L[itm1] * t2;
		N[itm1] = 0.5*(L2Lm1*tm1 - N[imid]);
		N[itp1] = 0.5*(L2Lm1*tp1 - N[imid]);
	}
	//tri midedge nodes:
	static int liv[3] = { 0, 1, 0 };
	static int ljv[3] = { 1, 2, 2 };
	double LI, LJ, lilj2;
	int I, J;
	for (itm1 = 6; itm1 < 9; itm1++)
	{
		itp1 = itm1 + 3;
		I = liv[itm1 - 6];
		J = ljv[itm1 - 6];
		LI = L[I];
		LJ = L[J];
		lilj2 = 2.0*LI*LJ;
		N[itm1] = lilj2*tm1;
		N[itp1] = lilj2*tp1;
	}
	return 15;
}


int SRelement::BrickShapeFns(double r, double s, double t, double N[])
{
	//calculate  shape functions at a natural coordinate point in a brick
	//see Zienkiewicz p 185
	//input:
		//r, s, t = natural coordinates
	//output:
		//N = shape functions
	//return:
		//number of functions

	int i;
	double r0, s0, t0;
	double r2m1 = 1.0 - r*r;
	double s2m1 = 1.0 - s*s;
	double t2m1 = 1.0 - t*t;
	//corner functions:
	for (i = 0; i < 8; i++)
	{
		r0 = r*model.map.GetBrickNodeCoord(i).d[0];
		s0 = s*model.map.GetBrickNodeCoord(i).d[1];
		t0 = t*model.map.GetBrickNodeCoord(i).d[2];
		N[i] = 0.125*(1.0 + r0)*(1.0 + s0)*(1.0 + t0)*(r0 + s0 + t0 - 2.0);
	}
	//mid-edge functions parallel to r:
	for (i = 8; i <= 14; i += 2)
	{
		s0 = s*model.map.GetBrickNodeCoord(i).d[1];
		t0 = t*model.map.GetBrickNodeCoord(i).d[2];
		N[i] = 0.25*r2m1*(1.0 + s0)*(1.0 + t0);
	}
	//mid-edge functions parallel to s:
	for (i = 9; i <= 15; i += 2)
	{
		r0 = r*model.map.GetBrickNodeCoord(i).d[0];
		t0 = t*model.map.GetBrickNodeCoord(i).d[2];
		N[i] = 0.25*s2m1*(1.0 + r0)*(1.0 + t0);
	}
	//mid-edge functions parallel to t:
	for (i = 16; i < 20; i++)
	{
		r0 = r*model.map.GetBrickNodeCoord(i).d[0];
		s0 = s*model.map.GetBrickNodeCoord(i).d[1];
		N[i] = 0.25*t2m1*(1.0 + r0)*(1.0 + s0);
	}
	return 20;
}
