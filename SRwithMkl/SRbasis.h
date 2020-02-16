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
// SRbasis.h: interface for the SRbasis class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRBASIS_INCLUDED)
#define SRBASIS_INCLUDED

enum SRBasisCallType { basisonly, derivonly, both };

class SRbasis  
{
public:
	int FaceGetpmax(SRface* face);
	int FaceGetpmax(int gface);
	int CountFaceTotalFunctions(SRelement* elem, int lface, int &nint);
	int CountFaceInternalFunctions(SRelement* elem, int lface);
	int CountFaceInternalFunctions(SRface* face);
	int CountFaceTotalFunctions(SRface* face, int &nint);
	int CountElementFunctions(SRelement* elem, int &nint);
	int CountWedgeFunctions(SRelement* elem, int &nint);
	int CountTetFunctions(SRelement* elem, int &nint);
	int CountBrickFunctions(SRelement* elem,int &nint);
	int CountFaceInternalFunctions(int numsides, int pmax);
	void Basisfuns1d(int pej, double rej, double* h);
	void Basisfuns1d(int pej, double rej, double* basisej, double* dbasisdrej, bool deriv = false);
	int FaceBasisFuncs(double r, double s, SRface* face, double* basisvec);
	int ElementBasisFuncs(double r, double s, double t, SRelement* elem, double* basisvec, double* dbasisdr, double* dbasisds, double* dbasisdt, SRBasisCallType calltype = both);
	int ElementBasisFuncs(double r, double s, double t, SRelement* elem, double* basisvec);
	int ElementBasisFuncs(double r, double s, double t, SRelement* elem, double* dbasisdr, double* dbasisds, double* dbasisdt);
	int TetBasisFuncs(double r, double s, double t, SRelement* elem, double* basisvec, double* dbasisdr, double* dbasisds, double* dbasisdt,	SRBasisCallType calltype);
	int WedgeBasisFuncs(double r, double s, double t, SRelement* elem, double* basisvec, double* dbasisdr, double* dbasisds, double* dbasisdt, SRBasisCallType calltype);
	int BrickBasisFuncs(double r, double s, double t, SRelement* elem, double* basisvec, double* dbasisdr, double* dbasisds, double* dbasisdt, SRBasisCallType calltype);
	int QuadFaceBasisFuncs(double r, double s, int* directionv, int* pejv, double* basisvec, SRface* face);
	int FaceBasisFuncs(double r, double s, int nej, int* directionv, int* pejv, double* basisvec, SRface* face);
	void TriEdgeBlend(int localedge, int direction, double r, double s, double l1, double l2, double l3, int& fun, int pej, double* basisvec);
	void QuadEdgeBlend(int localedge, int direction, double r, double s, int &fun, int pej, double* basisvec);
	void TriBubbleFuncs(double l1, double l2, double l3, int pmax, int& fun, double* basisvec);
	void QuadBubbleFuncs(int pmax, double r, double s, int& fun, double* basisvec);
	void LegendrePns(int pej, double r, double* legpns);
	void LegendrePnDerivs(int n, double* legpn, double* dpdr);
};

#endif // !defined(SRBASIS_INCLUDED)
