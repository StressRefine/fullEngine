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
// SRnode.h: interface for the SRnode class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRNODE_INCLUDED)
#define SRNODE_INCLUDED

class SRvec3;
class SRoutputData;

enum SRSacrificialNodeType{ notSacrificial, onBsurf, unsupported, shellOrBeamNode, breakoutCon, singular, ownedByBktEl };

class SRnode  
{
	friend class SRinput;

public:
	SRnode();
	void Create(int useridt, double xt, double yt, double zt){ userId = useridt; pos.Assign(xt, yt, zt); forceId = -1; };
	int GetId(){ return id; };
	int GetUserid(){ return userId; };
	SRvec3& Position(){ return pos; };
	void SetPosition(SRvec3 newPos){ pos.Copy(newPos); };
	double GetXyz(int i){ return pos.d[i]; };
	int GetGlobalFunctionNumber(){ return globalFunctionNumber; };
	void PutGlobalFunctionNumber(int g){globalFunctionNumber = g; };
	bool isMidSide() { return (midSideEdgeOwner != -1); };
	int GetMidSideEdgeOwner() { return midSideEdgeOwner; };
	bool isOrphan() { return (firstElementOwner == -1); };
	void SetAsMidside(int edgeId) { midSideEdgeOwner = edgeId; };
	void SetFirstElementOwner(int e) { firstElementOwner = e; };
	int GetFirstElementOwner() { return firstElementOwner; };
	SRconstraint* GetConstraint();
	void SetConstraintId(int i){ constraintId = i; };
	void SetUserId(int u){ userId = u; };
	bool GetDisp(SRvec3& disp);
	bool InsideBoundingBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
	void dumpData();
	bool checkUnsup();
	bool checkUnsupOrBktCon();
	bool checkSacr(){ return (sacrificialType != notSacrificial); };
	bool checkBktCon(){ return(sacrificialType == breakoutCon); };

	int id;
	int userId; //user original nodes numbers in case non-contiguous
	SRvec3 pos;
	int globalFunctionNumber;
	int midSideEdgeOwner;
	int firstElementOwner;
	int constraintId;
	SRSacrificialNodeType sacrificialType;
	int forceId;
	int dispCoordid;
	bool hasDisp;
};

#endif //!defined(SRNODE_INCLUDED)
