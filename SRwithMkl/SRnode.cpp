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
// SRnode.cpp: implementation of the SRnode class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SRmodel.h"

extern SRmodel model;

SRnode::SRnode()
{
	globalFunctionNumber = -1;
	userId = -1;
	firstElementOwner = -1;
	midSideEdgeOwner = -1;
	constraintId = -1;
	sacrificialType = notSacrificial;
	dispCoordid = -1;
	hasDisp = false;
}



SRconstraint* SRnode::GetConstraint()
{
	//return pointer to constraint acting on this node, NULL if not constrained
	if (constraintId == -1)
		return NULL;
	else
		return model.GetConstraint(constraintId);
}

bool SRnode::GetDisp(SRvec3& disp)
{
	if (isOrphan())
	{
		disp.Zero();
		return false;
	}
	if (isMidSide())
	{
		SRedge* edge = model.GetEdge(midSideEdgeOwner);
		edge->getDisp(0.0, disp);
		return true;
	}
	int gfun = globalFunctionNumber;
	for (int dof = 0; dof < 3; dof++)
	{
		double d = model.GetDisplacementCoeff(gfun, dof);
		disp.d[dof] = d;
	}
	return true;
}

bool SRnode::InsideBoundingBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
	double x = pos.d[0];
	double y = pos.d[1];
	double z = pos.d[2];
	if (x < xmin || x > xmax || y < ymin || y > ymax || z < zmin || z > zmax)
		return false;
	else
		return true;
}

bool SRnode::checkUnsupOrBktCon()
{
	if (sacrificialType == unsupported || sacrificialType == breakoutCon || sacrificialType == onBsurf || sacrificialType == ownedByBktEl)
		return true;
	if (model.breakout && model.faceBktConNodes.Get(id) == 1)
		return true;
	return false;
}

bool SRnode::checkUnsup()
{
	return sacrificialType == unsupported || sacrificialType == onBsurf;
}

void SRnode::dumpData()
{
	SRfile& f = model.dumpFile;
	f.PrintLine("%d %d %lg %lg %lg", id, userId, pos.d[0], pos.d[1], pos.d[2]);
}
