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
// SRface.h: interface for the SRface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRFACE_INCLUDED)
#define SRFACE_INCLUDED

#include "SRutil.h"

class SRlocalEdge;
class SRconstraint;
class SRelement;
class SRforce;

class SRfaceUtil
{
public:
	//static functions used during creation of global faces
	static int GetLocalEdge(int n1, int n2);
	static int GlobalFaceMatch(int& gn1, int& gn2, int& gn3, int& gn4, int n1, int n2, int n3, int n4);
	static bool QuadFaceMatch(int n1, int n2, int n3, int n4, int g1, int g2, int g3, int g4, int& gn1, int& gn2, int& gn3, int& gn4);
	static bool FaceMatch(int nv[], int gv[], int gnv[]);
	static bool FaceMatch(int n1, int n2, int n3, int n4, int g1, int g2, int g3, int g4, int& gn1, int& gn2, int& gn3, int& gn4);
	static void GetLocalEdgeNodes(int nej, int lej, int& n1, int& n2);
	static void QuadLocalEdgeNodes(int lej, int& n1, int& n2);
};

class SRface  
{

	friend class SRinput;
	friend class SRfaceUtil;

public:
	SRface();
	SRconstraint* GetConstraint();
	void GetLocalEdgeNodes(int lej, int &n1, int &n2);
	void QuadFaceNaturalCoordinatesFromEdge(int lej, double rej, double& rf, double& sf);
	void naturalCoordinatesFromEdge(int lej, double rej, double& rf, double& sf, bool useDirection = false);
	void Clear();
	void ProcessForce(SRforce *force, SRvec3& ResF);
	void GetForceValue(SRforce *force, double rf, double sf, double forceVal[3]);
	void GetSummedForceValue(double rf, double sf, double forceVal[3]);
	SRnode* GetNode(int localnodenum);

	SRedge* GetEdge(int localedgenum);
	void Create(int idt, int n1, int n2, int n3, int n4 = -1);
	int GetNumNodes(){ return (nodeIds[3] == -1 ? 3 : 4); };
	int GetNodeId(int i);
	int GetNumLocalEdges(){ return localEdges.GetNum(); };
	int GetLocalEdgePOrder(int lej){ return localEdges.Get(lej).GetPOrder(); }
	int GetLocalEdgeGlobalId(int lej){ return localEdges.Get(lej).GetGlobalEdgeId(); };
	int GetLocalEdgeDirection(int lej){	return localEdges.Get(lej).GetDirection(); };
	int GetNumGlobalFunctions(){ return globalFunctionNumbers.GetNum();	};
	int GetGlobalFunctionNumber(int i){ return globalFunctionNumbers.Get(i); };
	int GetConstraintId(){ return constraintId; };
	int GetElementOwner(int i){ return elementOwners[i]; };
	int GetElementLocalFace(int i){ return elementLocalFace[i]; };
	int GetMidNodeId(int lej){ return localEdges.GetPointer(lej)->GetMidNodeId(); };
	bool hasForce(){ return !forceIds.isEmpty(); };
	bool hasConstraint(){ return constraintId != -1; };
	double GetTractionJump(){ return tractionJump; };
	void SetTractionJump(double jump){ tractionJump = jump; };
	void NodeNaturalCoordinates(int lnode, double &rf, double &sf);
	void NaturalCoordinatesNearNode(int lnode, double &rf, double &sf);
	void LocalNodePosition(int lnode, SRvec3& pos);
	bool GetFlipNormal(){ return flipNormal; };
	void SetFlipNormal();
	void FillMappingNodes();
	double UnitTriad(double rf, double sf, SRvec3 &p, SRvec3 &e1, SRvec3 &e2, SRvec3 &e3, bool detjonly = false);
	void QuadShapeDerivs(double rf, double sf);
	void PutElementOwner(int n, int el){ elementOwners[n] = el; };
	void PutElementLocalFace(int n, int lf){ elementLocalFace[n] = lf; }
	void AssignLocalEdge(int lej, int gej, int direction){ localEdges.GetPointer(lej)->Assign(gej, direction); };
	bool IsBoundaryFace(){ return (elementOwners[1] == -1); };
	void AllocateGlobalFunctionNumbers(int n){ globalFunctionNumbers.Allocate(n); };
	void PutGlobalFunctionNumber(int lfun, int gfun){ globalFunctionNumbers.Put(lfun, gfun); }
	bool isFlat();
	void NaturalCoordinatesNearMidedge(int lej, double& r, double& s);
	void Position(double rf, double sf, SRvec3& p);
	void OutwardNormal(double rf, double sf, SRvec3& norm, bool checkOutwardLocally = true);
	void SetConstraintId(int c){ constraintId = c; };
	double Jacobian(double rf, double sf);
	double GetSize(){ return size; };
	int GetId(){ return id; };
	bool natCoordsAtCorner(int cornerId, double& r, double& s);
	bool isSaveForBreakout(){ return saveForBreakout; };
	void SetSaveForBreakout(int elOwner);
	int GetNumNodesTotal(){ return 2 * localEdges.GetNum(); };
	int midNodeMatch(int mid);
	double normalRotation();
	void deformedNormal(double rf, double sf, SRvec3& norm0, SRvec3 &norm);
	void RotateVec(int coordid, double rf, double sf, double forceVal[]);
	int GetNodeOrMidNodeId(int i);
	SRnode* GetNodeOrMidnode(int localnodenum);
	void OutwardNormalLinear(double rf, double sf, SRvec3& norm);
	double UnitTriadLinear(double rf, double sf, SRvec3 &p, SRvec3 &e1, SRvec3 &e2, SRvec3 &e3);
	void QuadLinearShapeDerivs(double rf, double sf);
	bool checkPartBdry(bool cornersOnly = false);

	void dumpData();

	int nodeIds[4];
	SRvector <SRlocalEdge> localEdges;
	SRintVector globalFunctionNumbers;
	int constraintId;
	int elementOwners[2];
	int elementLocalFace[2];
	SRintVector forceIds;
	double tractionJump;
	bool flipNormal;
	double xnode[8], ynode[8], znode[8];
	double dNdrf[8], dNdsf[8];
	double size;
	int flat;
	int id;
	int multifaceForceGroupId;
	bool saveForBreakout;
	bool unsupported;
	bool onPartBdry;
#ifdef _DEBUG
	int nodeUids[8]; //(for debugging)
#endif
};

#endif // !defined(SRFACE_INCLUDED)
