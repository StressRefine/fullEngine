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
// SRforce.h: interface for the SRforce class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRFORCE_INCLUDED)
#define SRFORCE_INCLUDED

class SRforce;
class SRvolumeForce;
class SRelement;
class SRvec3;
class SRdoubleVector;

enum SRforceType { nodalForce, edgeForce, faceForce, inactiveForce };
enum SRvolumeForceType { gravity, centrifugal };

class SRforce
{
	friend class SRinput;

public:
	SRforceType GetType(){ return type; };
	bool isPressure(){ return pressure; };
	bool isGcs(){ return coordId == -1; };
	int GetCoordId(){ return coordId; };
	int GetEntityId(){ return entityId; };
	void SetEntityId(int i){ entityId = i; };
	double GetForceVal(int i, int j){ return forceVals.Get(i, j); };
	int GetNumForces(){ return forceVals.getNumCols(); };
	void Copy(SRforce& that, bool copyForceVals = true);
	void dumpData();
	void Clear();
	void AddNodalForce(SRforce& that, bool SummingSets = false);

	SRforce();

	SRforceType type;
	int coordId;
	int entityId; //node, edge, or face.
	int elemId; //element id for pressure on faces
	int nv[4]; //nodes at corners of face for pressure on faces
	SRdoubleMatrix forceVals;
	bool pressure;
	bool faceFromNodal;
	SRdoubleVector StiffMat;
	SRintVector stiffDiag;
};


class SRthermalForce
{
	friend class SRinput;
public:
	double GetTemp(SRelement* elem, double r, double s, double t);
	bool CeTMult(SRmaterial* mat, double eTx, double eTy, double eTz, double ceT[]);
	void Process();
	bool constantTemp;
	double temp;
	SRdoubleVector nodalTemp;
};

class SRvolumeForce
{
	friend class SRinput;

public:
	SRvolumeForce(){ g1 = g2 = g3 = 0.0; omega2 = 0.0; }
	void GetForceValue(SRelement* elem, SRvec3& p, double val[]);
	SRvolumeForceType GetType(){ return type; };
	SRintVector elList;

	SRvolumeForceType type;
	//for gravity loads:
	double g1, g2, g3;
	//for centrifugal loads:
	double omega2;
	double alpha;
	SRvec3 axis;
	SRvec3 origin;
	int numSets;
};

#endif // !defined(SRFORCE_INCLUDED)
