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
// SRconstraint.cpp: implementation of the SRconstraint class.
//
//////////////////////////////////////////////////////////////////////

#include "SRmodel.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

extern SRmodel model;

SRconstraint::SRconstraint()
{
	for (int i = 0; i < 3; i++)
		constrainedDof[i] = false;
	coordId = -1;
	breakout = false;
}

void SRconstraint::Clear()
{
	for (int i = 0; i < 3; i++)
	{
		constrainedDof[i] = false;
	}
	coordId = -1;
	entityId = -1;
	coordId = -1;
	enforcedDisplacementData.Free();
}


void SRconstraint::Copy(SRconstraint& that)
{
	//copy contents of constraint "that" to this one
	//input:
		//that = constraint
	for (int i = 0; i < 3; i++)
		constrainedDof[i] = that.constrainedDof[i];
	if (!that.enforcedDisplacementData.isEmpty())
	{
		int n, m;
		that.enforcedDisplacementData.getSize(n, m);
		enforcedDisplacementData.Allocate(n,m);
		enforcedDisplacementData.Copy(that.enforcedDisplacementData);
	}
	entityId = that.entityId;
	type = that.type;
	coordId = that.coordId;
}

void SRconstraint::AddNodalConstraint(SRconstraint& that, bool SummingSets)
{
	//add contents of constraint "that" to this one; valid for nodal constraints only
	//input:
		//that = constraint
	if (type != nodalCon)
		return;

	bool thisEnforced = !enforcedDisplacementData.isEmpty();
	bool thatEnforced = !that.enforcedDisplacementData.isEmpty();
	if (thatEnforced && !thisEnforced)
	{
		thisEnforced = true;
		enforcedDisplacementData.Allocate(1, 3);
	}
	for (int dof = 0; dof < 3; dof++)
	{
		double thatenfd = 0.0;
		if (thatEnforced)
			thatenfd = that.enforcedDisplacementData.Get(0, dof);
		if (!constrainedDof[dof])
		{
			if (!that.constrainedDof[dof])
				continue;
			else
			{
				constrainedDof[dof] = true;
				enforcedDisplacementData.Put(0, dof, thatenfd);
			}
		}
		else
		{
			if (!that.constrainedDof[dof])
				continue;
			else if (thisEnforced && thatEnforced)
			{
				//only sum the enfds if they are not the same or from different sets:
				double enfd = enforcedDisplacementData.Get(0, dof);
				if ( (fabs(enfd - thatenfd) > TINY) || SummingSets)
				{
					enfd += thatenfd;
					enforcedDisplacementData.Put(0, dof, enfd);
				}
			}
		}
	}
}

void SRconstraint::ProcessFaceConstraint()
{
	//process a gcs face constraint

	//note:
		//constrains dofs above p2 to 0

	SRedge* edge;
	SRface* face;
	SRmat33 R;
	face = model.GetFace(entityId);

	for (int dof = 0; dof < 3; dof++)
	{
		if (constrainedDof[dof] == 0)
			continue;

		for (int g = 0; g < face->GetNumGlobalFunctions(); g++)
		{
			int gfun = face->GetGlobalFunctionNumber(g);
			model.PutFunctionEquation(gfun, dof, -1);
		}

		if (!hasEnforcedDisp())
			continue;

		double dispvec[8];//vector of coefficients for the p-1 and p-2 functions of the face. 8 is worst case for p2-quad face
		//since the displacement variation is prescribed to be quadratic, the functions above p2 are all constrained to 0.

		FillFaceEnforcedDispCoeffs(dof, dispvec);
		int n = face->GetNumNodes();
		for (int node = 0; node < n; node++)
		{
			int gfun = face->GetGlobalFunctionNumber(node);
			double enfd = dispvec[node];
			model.PutEnforcedDisplacement(gfun, dof, enfd);
		}

		for (int ledge = 0; ledge < face->GetNumLocalEdges(); ledge++)
		{
			edge = face->GetEdge(ledge);
			if (edge->GetPorder() < 2)
				continue;
			//edge p2 function
			int lfun = ledge + n;
			int gfun = edge->GetGlobalFunctionNumber(2);
			double enfd = dispvec[lfun];
			model.PutEnforcedDisplacement(gfun, dof, enfd);
		}
	}
}


SRcoord* SRconstraint::GetCoord()
{
	//look up coordinate system associated with a constraint, if any
	//return:
		// pointer to coordinate system, NULL if none
	if (coordId == -1)
		return NULL;
	else
		return model.GetCoord(coordId);
}


double SRconstraint::GetFaceEnforcedDisp(double rf, double sf, int dof)
{
	//Determine enforced displacement on a face at this location
	//input:
	//rf,sf = natural coordinates on face; ignored if nodeNum != -1
	//dof = degree of freedom number (0-2)
	//return:
	//enforced displacement
	int i, nn;
	double dispt = 0.0, N[8];
	SRface *face;
	i = entityId;
	face = model.GetFace(i);
	nn = 2 * face->GetNumLocalEdges();
	double dispVec[8];
	for (i = 0; i < nn; i++)
		dispVec[i] = enforcedDisplacementData.Get(i, dof);

	dispt = 0.0;
	model.map.FaceShapeFunctions(face, rf, sf, N);
	for (i = 0; i < nn; i++)
		dispt += dispVec[i] * N[i];
	return dispt;
}

double SRconstraint::GetEdgeEnforcedDisp(double r, int dof)
{
	//Determine enforced displacement on an edge at this location
	//input:
		//r = natural coordinate on edge; ignored if nodeNum != -1
		//dof = degree of freedom number (0-2)
	//return:
		//enforced displacement
	int i;
	double dispt = 0.0, N[3];

	dispt = 0.0;
	model.map.EdgeShapeFunctions(r, N);
	for (i = 0; i < 3; i++)
	{
		dispt += enforcedDisplacementData.Get(i, dof)*N[i];
	}
	return dispt;
}

void SRconstraint::FillFaceEnforcedDispCoeffs(int dof, double *dispvec)
{
	//fill vector of enforced displacement coefficients for face basis functions
	//input:
	//dof=degree-of-freedom number
	//output:
	//dispvec=enforced displacement coefficients for all functions of face

	SRface *face = model.GetFace(entityId);

	int nn = face->GetNumNodes();
	//nodes:
	for (int i = 0; i < nn; i++)
		dispvec[i] = GetEnforcedDisp(i, dof);

	//Fit edges by collocating p2 function at midedge. the enforced disp value there is the same as the corresponding
	//midside node value on the face.
	//At midedge, nodal bfs corresponding to the local edge are 1/2, all other nodal bf's are 0.

	int nc = face->GetNumLocalEdges();
	for (int i = 0; i < nc; i++)
	{
		double enfdisp = GetEnforcedDisp(i + nc, dof);
		dispvec[nn + i] = enfdisp;
	}
}

double SRconstraint::GetEnforcedDisp(int nodeNum, int dof)
{
	//look up enforced displacement value at a node of an edge or face constraint
	//input:
		//nodeNum = node number of corner of midside
		//dof = degree of freedom
	//return:
		//displacement value;
	double disp = enforcedDisplacementData.Get(nodeNum, dof);

	return disp;
}

void SRconstraint::PutEnforcedDisplacementData(int n, int dof, double val)
{
	if (enforcedDisplacementData.isEmpty())
		enforcedDisplacementData.Allocate(n + 1, 3);
	enforcedDisplacementData.Put(n, dof, val);
};

bool SRconstraint::allDofsConstrained()
{
	for (int i = 0; i < 3; i++)
	{
		if (constrainedDof[i] == 0)
			return false;
	}
	return true;
}

void SRconstraint::dumpData()
{
	SRfile& f = model.dumpFile;
	f.Print("%d ", entityId);
	if (type == nodalCon)
		f.Print(" node constrainedDofs: ");
	else if (type == faceCon)
	{
		SRface* face = model.GetFace(entityId);
		f.Print(" face ");
		for (int i = 0; i < face->GetNumNodes(); i++)
			f.Print(" %d", face->GetNode(i)->GetUserid());
		f.Print(" constrainedDofs: ");
	}
	for (int dof = 0; dof < 3; dof++)
		f.Print(" %d", constrainedDof[dof]);
	f.PrintReturn();
	if (coordId == -1)
		f.PrintLine("gcs");
	else
		f.PrintLine("%s", GetCoord()->GetName());
	if (hasEnforcedDisp())
	{
		f.PrintLine(" enforced: ");
		for (int i = 0; i < enforcedDisplacementData.getNumCols(); i++)
		{
			for (int dof = 0; dof < 3; dof++)
					f.Print("   %lg", GetEnforcedDisp(i, dof));
			f.PrintReturn();
		}
	}
	f.PrintReturn();
}

void SRconstraint::getDisp(int n, SRvec3& enfd)
{
	if (enforcedDisplacementData.isEmpty())
		enfd.Zero();
	else
	{
		for (int dof = 0; dof < 3; dof++)
			enfd.d[dof] = enforcedDisplacementData.Get(n, dof);
	}
}

double SRconstraint::getDisp(int n, int d)
{
	if (!enforcedDisplacementData.isEmpty())
		return enforcedDisplacementData.Get(n, d);
	else
		return 0.0;
}
