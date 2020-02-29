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
// SRfaceQuad.cpp: implementation of the SRface class
//	for quad faces
//
//////////////////////////////////////////////////////////////////////

#include "SRmodel.h"

extern SRmodel model;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

///////////////////////////////////////////////////////////////////////
// SRfaceUil static classes are used during creation of global faces //
///////////////////////////////////////////////////////////////////////

bool SRfaceUtil::QuadFaceMatch(int n1, int n2, int n3, int n4, int g1, int g2, int g3, int g4, int &gn1, int &gn2, int &gn3, int &gn4)
{
	//see if quad face with corner nodes n1,n2,n3,n4 matches face with corner nodes g1,g2,g3,g4
	//input:
		//n1,n2,n3,n4 = nodes at corners of quad face
		//g1,g2,g3,g4 = nodes at corners of existing global face
	//output:
		//gn1,gn2,gn3,gn4 = node order on face g1,g2,g3,(g4)
	//return:
		//true if match else false

	if (n1 == g1)
	{
		gn1 = 0;
		if (n2 == g2 && n3 == g3 && n4 == g4)
		{
			gn2 = 1;
			gn3 = 2;
			gn4 = 3;
			return true;
		}
		if (n2 == g4 && n3 == g3 && n4 == g4)
		{
			gn2 = 3;
			gn3 = 2;
			gn4 = 1;
			return true;
		}
	}
	else if (n1 == g2)
	{
		gn1 = 1;
		if (n2 == g3 && n3 == g4 && n4 == g1)
		{
			gn2 = 2;
			gn3 = 3;
			gn4 = 0;
			return true;
		}
		if (n2 == g1 && n3 == g4 && n4 == g3)
		{
			gn2 = 0;
			gn3 = 3;
			gn4 = 2;
			return true;
		}
	}
	else if (n1 == g3)
	{
		gn1 = 2;
		if (n2 == g4 && n3 == g1 && n4 == g2)
		{
			gn2 = 3;
			gn3 = 0;
			gn4 = 1;
			return true;
		}
		if (n2 == g2 && n3 == g1 && n4 == g4)
		{
			gn2 = 1;
			gn3 = 0;
			gn4 = 3;
			return true;
		}
	}
	else if (n1 == g4)
	{
		gn1 = 3;
		if (n2 == g3 && n3 == g2 && n4 == g1)
		{
			gn2 = 2;
			gn3 = 1;
			gn4 = 0;
			return true;
		}
		if (n2 == g1 && n3 == g2 && n4 == g3)
		{
			gn2 = 0;
			gn3 = 1;
			gn4 = 2;
			return true;
		}
	}
	return false;
}

void SRfaceUtil::QuadLocalEdgeNodes(int lej, int &n1, int &n2)
{
	//look up local node numbers at ends of a local edge of a quad face
	//input:
		//lej = local edge number
	//output:
		//n1,n2 = local node numbers on the face at ends 1 and 2 of the edge
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
		n1 = 3;
		n2 = 2;
	}
	else if (lej == 3)
	{
		n1 = 0;
		n2 = 3;
	}
}


///////////////////////////////////////////////////////////////
//               SRface member functions                    //
///////////////////////////////////////////////////////////////


void SRface::QuadFaceNaturalCoordinatesFromEdge(int lej, double rej, double &rf, double &sf)
{
	//calculate natural coordinates on quad face given natural coordinate on local edge
	//input:
		//lej = local edge number
        //rej = natural coordinate on local edge
	//output:
		//rf, sf = natural coordinates on quad face
	if (lej == 0)
	{
		rf = rej;
		sf = -1.0;
	}
	else if (lej == 1)
	{
		sf = rej;
		rf = 1.0;
	}
	else if (lej == 2)
	{
		rf = rej;
		sf = 1.0;
	}
	else if (lej == 3)
	{
		sf = rej;
		rf = -1.0;
	}
}

void SRface::QuadShapeDerivs(double rf, double sf)
{
	//support routine for FaceUnitTriad for quad faces.
	//calculate dNdrf, dNdsf at rf,sf
	//input:
		//rf, sf = natural coordinates
	//note:
		//fills up class variables dNdrf, dNdsf

	int i;
	double r0, s0, r2, s2, ri, si;
	r2 = rf*rf;
	s2 = sf*sf;
	//corners:
	for (i = 0; i < 4; i++)
	{
		ri = model.map.GetQuadNoder(i);
		si = model.map.GetQuadNodes(i);
		r0 = rf*ri;
		s0 = sf*si;
		dNdrf[i] = 0.25*ri*(1.0 + s0)*(2.0*r0 + s0);
		dNdsf[i] = 0.25*si*(1.0 + r0)*(2.0*s0 + r0);
	}
	//nodes5,7 = midnodes of edges on which r varies
	for (i = 4; i <= 6; i += 2)
	{
		si = model.map.GetQuadNodes(i);
		s0 = sf*si;
		dNdrf[i] = -rf*(1.0 + s0);
		dNdsf[i] = 0.5*si*(1.0 - r2);
	}
	//nodes6,8 = midnodes of edges on which s varies
	for (i = 5; i <= 7; i += 2)
	{
		ri = model.map.GetQuadNoder(i);
		r0 = rf*ri;
		dNdsf[i] = -sf*(1.0 + r0);
		dNdrf[i] = 0.5*ri*(1.0 - s2);
	}
}

void SRface::dumpData()
{
	SRfile& f = model.dumpFile;
	f.Print("%d ", id);
	for (int i = 0; i < GetNumNodes(); i++)
		f.Print(" %d", nodeIds[i]);
	for (int i = 0; i < localEdges.GetNum(); i++)
		f.Print(" %d", GetEdge(i)->GetMidNodeId());
	for (int i = 0; i < forceIds.GetNum(); i++)
		f.Print(" %d forced id %d", i, forceIds.Get(i));
	f.Print(" %d", constraintId);
	f.PrintReturn();
}


