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
// SRedge.h: interface for the SRedge class.
//
//
//////////////////////////////////////////////////////////////////////

	

#if !defined(SREDGE_INCLUDED)
#define SREDGE_INCLUDED

class SRoutputData;
struct SRpostEdgeData;

class SRedge;
class SRforce;

class SRedgeUtil
{
public:
	//static functions used during creation of global edges
	static bool EdgeMatch(int n1, int n2, int gn1, int gn2, int& direction);
	static int GlobalEdgeMatchN2(int n1, int n2, int& direction);
	static int GlobalEdgeMatch(int n1, int n2, int& direction);
};


class SRlocalEdge
{
	friend class SRinput;

public:
	SRlocalEdge::SRlocalEdge();
	int GetPOrder();
	SRedge* GetEdge();
	int GetGlobalEdgeId(){ return globalEdgeId; };
	int GetMidNodeId();
	int GetDirection(){ return direction; };
	void Assign(int id, int directiont){ globalEdgeId = id; direction = directiont; };
	int globalEdgeId;
	int direction;
};

class SRboundaryFaceData
{
public:
	int faceId;
	int localEdgeId;
};

class SRedge
{
	friend class SRedgeUtil;
	friend class SRinput;
public:
	SRconstraint* GetConstraint();
	void GetDisplacement(double r, SRvec3 &disp);
	void Clear();
	void GetForceValue(double r, SRforce* force, double forceVal[]);
	void Create(int n1, int n2, int& midnode, int idt);
	SRnode* GetNode(int localnodenum);
	int GetNodeId(int i);
	int GetMidNodeId(){ return midnodeId; };
	int GetNodeOrMidNodeId(int n)
	{
		if (n < 2)
			return GetNodeId(n);
		else
			return midnodeId;
	};
	int GetMidNodeUserId();
	int GetPorder(){ return pOrder; };
	int GetPrevPorder(){ return prevpOrder; };
	int GetNumGlobalFunctions(){ return globalFunctionNumbers.GetNum(); };
	int GetGlobalFunctionNumber(int i){ return globalFunctionNumbers.Get(i); };
	void FillMappingNodes();
	void AllocateGlobalFunctionNumbers(int n){ globalFunctionNumbers.Allocate(n); };
	void putPorder(int p);
	void AssignGlobalFunctionNumbers(int& fun, int pmin, int pmax);
	bool isStraight(){ return straight; };
	void Position(double r, SRvec3& p);
	bool PChanged(){ return pOrder != prevpOrder; };
	int GetConstraintId(){ return constraintId; };
	void SetConstraintId(int c){ constraintId = c; };
	double GetSize(){ return size; };
	void SetThin(bool thinin){ thin = thinin; };
	void Straighten(double fraction = 1.0);
	void ScaleStraightenVsEdgeLength();
	void basisFunctions(double r, double* basis);
	bool isSacrificial(){ return (sacrificial != 0); };
	void setSacrificial(int pf = 1){ sacrificial = pf; }
	void getDisp(double r, SRvec3& disp);
	int GetId(){ return id; };
	double getXnode(int i){ return xnode[i]; };
	double getYnode(int i){ return ynode[i]; };
	double getZnode(int i){ return znode[i]; };
	bool isOrphan();
	void ProcessForce(SRforce* force, SRvec3& ResF);
	void PutBoundaryFaceId(int f, int lej);
	bool checkKinkOK(SRface* face0, int lej0, SRface* face1, int lej1);
	void dumpData();

	SRedge();
	int id;
	int nodeIds[2];
	int pOrder;
	int prevpOrder;
	SRintVector globalFunctionNumbers;
	int midnodeId;
	int constraintId;
	SRvector <SRboundaryFaceData> boundaryfaceData;
	double xnode[3], ynode[3], znode[3];
	double size;
	bool straight;
	bool straightened;
	bool thin;
	int sacrificial;
	int forceId;
	double initialBow;
	double straightenFraction;
};

#endif //!defined(SREDGE_INCLUDED)
