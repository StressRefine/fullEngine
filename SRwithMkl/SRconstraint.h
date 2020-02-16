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
// SRconstraint.h: interface for the SRconstraint class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRCONSTRAINT_INCLUDED)
#define SRCONSTRAINT_INCLUDED

enum SRconstraintType{ nodalCon, faceCon, inactiveCon };

class SRcoord;
class SRintVector;

class SRconstraint
{

	friend class SRinput;

public:
	void ProcessFaceConstraint();
	void FillFaceEnforcedDispCoeffs(int dof, double *dispvec);
	bool IsConstrainedDof(int i){ return constrainedDof[i]; };
	void SetConstrainedDof(int i){ constrainedDof[i] = true; };
	int GetEntityId(){ return entityId; };
	SRconstraintType GetType(){ return type; };
	SRcoord* GetCoord();
	bool hasEnforcedDisp(){	return !enforcedDisplacementData.isEmpty();	};
	int GetCoordId(){ return coordId; };
	bool isGcs(){ return (coordId == -1); };
	void Copy(SRconstraint& that);
	void AddNodalConstraint(SRconstraint& that, bool SummingSets = false);
	void operator = (SRconstraint& c2){ Copy(c2); };
	void getDisp(int n, SRvec3& enfd);
	double getDisp(int n, int d);
	void allocateEnforcedDisplacementData(int n){ enforcedDisplacementData.Allocate(n, 3); };
	void PutEnforcedDisplacementData(int n, int dof, double val);
	double GetEdgeEnforcedDisp(double re, int dof);
	double GetFaceEnforcedDisp(double rf, double sf, int dof);
	double GetEnforcedDisp(int nodeNum, int dof);
	int GetNumEnforcedDisp(){ return enforcedDisplacementData.getNumCols(); };
	bool allDofsConstrained();
	void Clear();
	bool isBreakout(){ return breakout; };
	void SetBreakout(){ breakout = true; };
	void SetType(SRconstraintType typeIn){ type = typeIn; };
	void SetEntityId(int i){ entityId = i; };
	void dumpData();


	SRconstraint();

	bool constrainedDof[3];
	bool breakout;
	int id;
	int entityId; //node, edge, or face
	SRconstraintType type;
	int coordId;
	SRdoubleMatrix enforcedDisplacementData;
};


#endif // !defined(SRCONSTRAINT_INCLUDED)
