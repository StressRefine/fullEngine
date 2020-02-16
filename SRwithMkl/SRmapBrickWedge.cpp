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
// SRmapBrickWedge.cpp: implementation of the SRmap class
//	for bricks and wedges
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SRmodel.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

extern SRmodel model;

void SRmap::BrickWedgeQuadInit()
{
	//miscellaneous initialization for mapping of quads, wedges, and bricks

	//r,s at mapping points of quad. faces:
	//corners:
	quadnoder[0] = -1.0;
	quadnoder[1] = 1.0;
	quadnoder[2] = 1.0;
	quadnoder[3] = -1.0;
	quadnodes[0] = -1.0;
	quadnodes[1] = -1.0;
	quadnodes[2] = 1.0;
	quadnodes[3] = 1.0;
	//midnodes:
	quadnoder[4] = 0.0;
	quadnoder[5] = 1.0;
	quadnoder[6] = 0.0;
	quadnoder[7] = -1.0;
	quadnodes[4] = -1.0;
	quadnodes[5] = 0.0;
	quadnodes[6] = 1.0;
	quadnodes[7] = 0.0;

	//r,s,t at corners and edge midnodes of brick elements
	//corners
	int i;
	for(i = 0; i < 4; i++)
	{
		brickNodes[i].Assign(quadnoder[i], quadnodes[i], -1.0);
		brickNodes[i + 4].Assign(quadnoder[i], quadnodes[i], 1.0);
	}
	//midnodes:
	AssignBrickMidNodes(8, 0, 1);//mid-edge 0-1
	AssignBrickMidNodes(9, 1, 2);//mid-edge 1-2
	AssignBrickMidNodes(10, 2, 3);//mid-edge 2-3
	AssignBrickMidNodes(11, 0, 3);//mid-edge 0-3
	AssignBrickMidNodes(12, 4, 5);//mid-edge 4-5
	AssignBrickMidNodes(13, 5, 6);//mid-edge 5-6
	AssignBrickMidNodes(14, 6, 7);//mid-edge 6-7
	AssignBrickMidNodes(15, 4, 7);//mid-edge 4-7
	AssignBrickMidNodes(16, 0, 4);//mid-edge 0-4
	AssignBrickMidNodes(17, 1, 5);//mid-edge 1-5
	AssignBrickMidNodes(18, 2, 6);//mid-edge 2-6
	AssignBrickMidNodes(19, 3, 7);//mid-edge 3-7

	//r,s,t of wedge corners:
	wedgeCorners[0].Assign(trinoder[0], trinodes[0], -1.0);
	wedgeCorners[1].Assign(trinoder[1], trinodes[1], -1.0);
	wedgeCorners[2].Assign(trinoder[2], trinodes[2], -1.0);
	wedgeCorners[3].Assign(trinoder[0], trinodes[0], 1.0);
	wedgeCorners[4].Assign(trinoder[1], trinodes[1], 1.0);
	wedgeCorners[5].Assign(trinoder[2], trinodes[2], 1.0);

	//local node numbers at ends of edges for wedges:
	wedgeEdgeLocalNodes[0][0] = 0;
	wedgeEdgeLocalNodes[0][1] = 1;
	wedgeEdgeLocalNodes[1][0] = 1;
	wedgeEdgeLocalNodes[1][1] = 2;
	wedgeEdgeLocalNodes[2][0] = 0;
	wedgeEdgeLocalNodes[2][1] = 2;
	wedgeEdgeLocalNodes[3][0] = 3;
	wedgeEdgeLocalNodes[3][1] = 4;
	wedgeEdgeLocalNodes[4][0] = 4;
	wedgeEdgeLocalNodes[4][1] = 5;
	wedgeEdgeLocalNodes[5][0] = 3;
	wedgeEdgeLocalNodes[5][1] = 5;
	wedgeEdgeLocalNodes[6][0] = 0;
	wedgeEdgeLocalNodes[6][1] = 3;
	wedgeEdgeLocalNodes[7][0] = 1;
	wedgeEdgeLocalNodes[7][1] = 4;
	wedgeEdgeLocalNodes[8][0] = 2;
	wedgeEdgeLocalNodes[8][1] = 5;

	//local node numbers at ends of edges for bricks:
	brickEdgeLocalNodes[0][0] = 0;
	brickEdgeLocalNodes[0][1] = 1;
	brickEdgeLocalNodes[1][0] = 1;
	brickEdgeLocalNodes[1][1] = 2;
	brickEdgeLocalNodes[2][0] = 3;
	brickEdgeLocalNodes[2][1] = 2;
	brickEdgeLocalNodes[3][0] = 0;
	brickEdgeLocalNodes[3][1] = 3;
	brickEdgeLocalNodes[4][0] = 4;
	brickEdgeLocalNodes[4][1] = 5;
	brickEdgeLocalNodes[5][0] = 5;
	brickEdgeLocalNodes[5][1] = 6;
	brickEdgeLocalNodes[6][0] = 7;
	brickEdgeLocalNodes[6][1] = 6;
	brickEdgeLocalNodes[7][0] = 4;
	brickEdgeLocalNodes[7][1] = 7;
	brickEdgeLocalNodes[8][0] = 0;
	brickEdgeLocalNodes[8][1] = 4;
	brickEdgeLocalNodes[9][0] = 1;
	brickEdgeLocalNodes[9][1] = 5;
	brickEdgeLocalNodes[10][0] = 2;
	brickEdgeLocalNodes[10][1] = 6;
	brickEdgeLocalNodes[11][0] = 3;
	brickEdgeLocalNodes[11][1] = 7;
}



void SRmap::AssignBrickMidNodes(int midnode, int node0, int node1)
{
	//assign brick midnode r,s,t values
	double r, s, t;
	r = 0.5*(brickNodes[node0].d[0] + brickNodes[node1].d[0]);
	s = 0.5*(brickNodes[node0].d[1] + brickNodes[node1].d[1]);
	t = 0.5*(brickNodes[node0].d[2] + brickNodes[node1].d[2]);
	brickNodes[midnode].Assign(r, s, t);
}


void SRmap::QuadShapeFunctions(double rf, double sf, double N[])
{
	//calculate shape functions for a quadrilateral face
	//input:
		//rf,sf=natural coordinates on face
	//output:
		//N=vector of shape functions

	double r0, s0;
	int i;


	//quadratic mapping (Zienkiewicz, p. 174):
	//corners:
	for (i = 0; i < 4; i++)
	{
		r0 = rf*quadnoder[i];
		s0 = sf*quadnodes[i];
		N[i] = 0.25*(1.0 + r0)*(1.0 + s0)*(r0 + s0 - 1.0);
	}
	double r2, s2;
	r2 = rf*rf;
	s2 = sf*sf;
	double r2m1 = 0.5*(1.0 - r2);
	double s2m1 = 0.5*(1.0 - s2);
	//nodes5,7 = midnodes of edges on which r varies
	for (i = 4; i <= 6; i += 2)
	{
		s0 = sf*quadnodes[i];
		N[i] = r2m1*(1.0 + s0);
	}
	//nodes6,8 =midnodes of edges on which s varies
	for (i = 5; i <= 7; i += 2)
	{
		r0 = rf*quadnoder[i];
		N[i] = s2m1*(1.0 + r0);
	}
}

void SRmap::QuadLinearShapeFunctions(double rf, double sf, double N[])
{
	//calculate linear shape functions for a quadrilateral face
	//input:
		//rf,sf=natural coordinates on face
	//output:
		//N=vector of shape functions

	double r0, s0;
	int i;

	for (int i = 0; i<4; i++)
	{
		r0 = rf*quadnoder[i];
		s0 = sf*quadnodes[i];
		N[i] = 0.25*(1.0 + r0)*(1.0 + s0);
	}
}

void SRmap::WedgeLinearShapeFunctions(double r, double s, double t, double N[])
{
	//calculate linear shape functions for a wedge
	//input:
	//rf,sf=natural coordinates on face
	//output:
	//N=vector of shape functions
	double r0, s0, t0, ri, si, ti;
	double L[3];
	int i;
	//derivatives of corner shape functions:
	double tp1 = 0.5*(1.0 + t);
	double tm1 = 0.5*(1.0 - t);
	TriangleAreaCoords(r, s, L);
	for (i = 0; i < 3; i++)
	{
		N[i] = L[i] * tm1;
		N[i + 3] = L[i] * tp1;
	}
}

void SRmap::BrickLinearShapeFunctions(double r, double s, double t, double N[])
{
	//calculate  shape functions at a natural coordinate point in a brick
	//input:
	//r, s, t = natural coordinates
	//output:
	//N = shape functions

	int i;
	double r0, s0, t0;
	for (i = 0; i < 8; i++)
	{
		r0 = r*brickNodes[i].d[0];
		s0 = s*brickNodes[i].d[1];
		t0 = t*brickNodes[i].d[2];
		N[i] = 0.125*(1.0 + r0)*(1.0 + s0)*(1.0 + t0);
	}
}
