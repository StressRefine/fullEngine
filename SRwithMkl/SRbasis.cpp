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
// SRbasis.cpp: implementation of the SRbasis class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SRmodel.h"
#include "SRbasis.h"


extern SRmodel model;


void SRbasis::LegendrePns(int p, double r, double *Pn)
{
	//Legendre Polynomials up to order p
	//input:
		//p = polynomial order
		//r = coordinate from -1 to 1
	//output:
		//Pn = vector of Legendre Polynomials at r
	int i;
	Pn[0] = 1.0;
	Pn[1] = r;
	for(i = 2; i <= p; i++)
	{
		//see Szabo regression formula p 355.
		Pn[i] = ((2.0*i - 1.0)*r*Pn[i - 1] - (i - 1.0)*Pn[i - 2]) / i;
	}
}

void SRbasis::Basisfuns1d(int pej, double rej, double *basisej, double *dbasisdr, bool deriv)
{
	//1d basis functions on an edge at rej
	//input:
		//pej = polynomial order for edge
		//rej = natural coordinate on edge
		//deriv = true if derivatives of basis functions desired else false
	//output:
		//basisej = vector of basis functions
		//dbasisdr = vector of derivatives of basis functions if deriv=true else not used
	double P[MAX1DFUNCTIONS];
	int i;
	LegendrePns(pej,rej,P);
	basisej[0] = 0.5*(1.0 - rej);
	basisej[1] = 0.5*(1.0 + rej);
	for (i = 2; i <= pej; i++)
	{
		//see Szabo, p 38. multiplicative factor is 1/(sqrt(2*(2j-1);
		//I factored out 1/sqrt(2), leaving 1/sqrt(2j-1).
		//handle zero base: my basisej[i] = his N(i+1). but N(i+1) = phi(i)
		//(formula 3,12a) so basisej[i] = phi(i)
		basisej[i] = (P[i] - P[i - 2])*SQRT2OVER2 / (sqrt(2.0*i - 1.0));
	}
	if (deriv)
	{
		//Szabo 3.12b
		dbasisdr[0] = -0.5;
		dbasisdr[1] = 0.5;
		for (i = 2; i <= pej; i++)
			dbasisdr[i] = SQRT2OVER2*sqrt(2.0*i - 1.0)*P[i - 1];
	}

}

void SRbasis::Basisfuns1d(int pej, double rej, double *h)
{
	//1d basis functions on an edge at rej, no derivatives
	//input:
		//pej = polynomial order for edge
		//rej = natural coordinate on edge
	//output:
		//basisvec = vector of basis functions
	
	double tmp;
	Basisfuns1d(pej,rej,h,&tmp,false);
}

void SRbasis::TriEdgeBlend(int localedge, int direction, double r, double s, double l1, double l2, double l3, int &fun, int pej, double *basisvec)
{
	//blend edge basis functions to form element functions for triangle
	//input:
		//localedge = local edge number
		//direction = 1 if local edge direction matches global direction else -1
		//r,s = natural coordinates on element
		//l1,l2,l3 = area coordinates for triangle
		//fun = function number to use to store next element function
		//pej = polynomial order for this edge
	//output:
		//fun: updated
		//basisvec = face basis function vector

	double basisej[MAX1DFUNCTIONS],lv[3];
	double rej,LI,LJ;
	int i,nodeI,nodeJ;
	lv[0] = l1;
	lv[1] = l2;
	lv[2] = l3;
	nodeI = model.map.triEdgeLocalNodes[localedge][0];
	nodeJ = model.map.triEdgeLocalNodes[localedge][1];
	LJ = lv[nodeJ];
	LI = lv[nodeI];
	rej = LJ - LI;
	rej *= direction;

	//edge basis functions go to zero at corners. catch this case explicitly
	if(((1.0 - rej) < RELSMALL) || ((rej + 1.0) < RELSMALL) )
	{
		for(i = 3; i <= pej; i++)
		{
			basisvec[fun] = 0.0;
			fun++;
		}
	}
	else
	{
		Basisfuns1d(pej, rej, basisej);

		double tmp = 4.0*(LI / (1.0 - rej))*(LJ / (1.0 + rej));
		for (i = 3; i <= pej; i++)
		{
			basisvec[fun] = tmp*basisej[i];
			fun++;
		}
	}
}

void SRbasis::TriBubbleFuncs(double l1, double l2, double l3,int pmax, int &fun, double *basisvec)
{
	//bubble basis functions for triangular face
	//input:
		//l1,l2,l3 = area coordinates
		//pmax = max polynomial order on any edge of the face
		//fun = function number to use to store next element function
	//output:
		//fun: updated
		//basisvec = face basis function vector

	if(pmax < 3)
		return;
	double lll = l1*l2*l3;
	int i,j,n;
	double r = l2 - l1;
	double s11 = 2.0*l3 - 1.0; //s11 = coordinate in s direction scaled from -1 to 1
	double Pr[MAX1DFUNCTIONS],Ps[MAX1DFUNCTIONS];
	n = pmax - 3;
	LegendrePns(pmax,r,Pr);
	LegendrePns(pmax,s11,Ps);

	int fun0 = fun;
	for(i = 0; i <= n; i++)
	{
		for(j = 0; j <= n;j++)
		{
			if((i+j) > n)
				break;
			basisvec[fun] = lll*Pr[i] * Ps[j];
			fun++;
		}
	}
}

int SRbasis::FaceBasisFuncs(double r, double s, int nej, int *directionv, int *pejv, double *basisvec, SRface* face)
{
	//basis functions for a tri or quad face
	//input:
		//r,s = natural coordinates on face
		//nej = no. of edges of face
		//directionv = direction for each edge (-1 or +1)
		//pejv = polynomial orders for each edge
	//output:
		//basisvec = vector of basis functions
	//return:
		//total number of basis functions;
	int edge, fun;
	int pmax = 0;
	if (nej == 3)
	{
		double N[6];
		model.map.FaceShapeFunctions(face, r, s, N);
		double l1, l2, l3;
		model.map.TriangleAreaCoords(r, s, l1, l2, l3);
		//corner and midedge functions:
		for (int i = 0; i < 6; i++)
			basisvec[i] = N[i];
		fun = 6;
		//higher edge functions:
		for (edge = 0; edge < 3; edge++)
		{
			TriEdgeBlend(edge, directionv[edge], r, s, l1, l2, l3, fun, pejv[edge], basisvec);

			if (pejv[edge] > pmax)
				pmax = pejv[edge];
		}
		//bubble functions:
		TriBubbleFuncs(l1, l2, l3, pmax, fun, basisvec);
	}
	else
		fun = QuadFaceBasisFuncs(r, s, directionv, pejv, basisvec, face);

	return fun;
}

int SRbasis::FaceBasisFuncs(double r, double s, SRface *face, double *basisvec)
{
	//basis functions for a tri or quad face
	//input:
		//r,s = natural coordinates on face
		//face = pointer to face
	//output:
		//basisvec = vector of basis functions
	//return:
		//total number of basis functions;

	int i, nej, directionv[4];
	int pejv[4];
	SRedge* edge;
	nej = face->GetNumLocalEdges();
	for(i = 0; i < nej; i++)
	{
		directionv[i] = face->GetLocalEdgeDirection(i);
		edge = face->GetEdge(i);
		pejv[i] = edge->pOrder;
	}
	return FaceBasisFuncs(r,s,nej,directionv,pejv,basisvec, face);
}

int SRbasis::ElementBasisFuncs(double r, double s, double t, SRelement *elem, double *basisvec, double *dbasisdr, double *dbasisds, double *dbasisdt, SRBasisCallType calltype)
{
	//basis functions and their derivatives at r,s,t for 
	//tet, wedge, and brick elements
	//input:
		//r,s,t = natural coordinates in the elements
		//calltype = "basisonly" if don't need derivatives,
		//			 "derivonly" if don't need basis function, 
 		//        or "both"
	//output:
		//basisvec = vector of basis functions
		//dbasisdr,dbasisds,dbasisdt = vector of basis function derivatives
	//return:
		//total number of basis functions;
	int nfunc;

	if ((elem->type) == tet)
		nfunc = TetBasisFuncs(r, s, t, elem, basisvec, dbasisdr, dbasisds, dbasisdt, calltype);
	else if ((elem->type) == wedge)
		nfunc = WedgeBasisFuncs(r, s, t, elem, basisvec, dbasisdr, dbasisds, dbasisdt, calltype);
	else if ((elem->type) == brick)
		nfunc = BrickBasisFuncs(r, s, t, elem, basisvec, dbasisdr, dbasisds, dbasisdt, calltype);

	int n = elem->GetNumFunctions();
	SRASSERT(nfunc == n);

	return nfunc;
}

int SRbasis::ElementBasisFuncs(double r, double s, double t, SRelement *elem, double *dbasisdr, double *dbasisds, double *dbasisdt)
{
	//basis function derivatives at r,s,t for an element
	//input:
		//r,s,t = natural coordinates in the elements
	//output:
		//dbasisdr,dbasisds,dbasisdt = vector of basis function derivatives
	//return:
		//total number of basis functions;
	double ht;
	return ElementBasisFuncs(r, s, t, elem, &ht, dbasisdr, dbasisds, dbasisdt, derivonly);
}

int SRbasis::ElementBasisFuncs(double r, double s, double t, SRelement *elem, double *basisvec)
{
	//basis functions at r,s,t for an element
	//input:
		//r,s,t = natural coordinates in the elements
	//output:
		//basisvec = vector of basis functions
	//return:
		//total number of basis functions;
	double hrt, hst, htt;
	return ElementBasisFuncs(r, s, t, elem, basisvec, &hrt, &hst, &htt, basisonly);
}

int SRbasis::TetBasisFuncs(double r, double s, double t, SRelement *elem, double *basisvec, double *dbasisdr, double *dbasisds, double *dbasisdt, SRBasisCallType calltype)
{
	//basis functions and their derivatives at r,s,t for a tetrahedral element
	//input:
		//r,s,t = natural coordinates in the elements
		//calltype = basisonly if don't need derivatives, derivonly if don't need basis function, or both
	//output:
		//basisvec = vector of basis functions
		//dbasisdr,dbasisds,dbasisdt = vector of basis function derivatives
	//return:
		//total number of basis functions

	//nodes 1 and 2 of local edges of tet:
	double lv[4], LI, LJ, LK, rej, basisejt;
	double dLIdr, dLIds, dLIdt, dLJdr, dLJds, dLJdt, dLKdr, dLKds, dLKdt, pipj;
	double basisej[MAX1DFUNCTIONS], dbasisdrej[MAX1DFUNCTIONS];
	double Pr[MAX1DFUNCTIONS], Ps[MAX1DFUNCTIONS], Pt[MAX1DFUNCTIONS], dpdr[MAX1DFUNCTIONS], dpds[MAX1DFUNCTIONS], dpdt[MAX1DFUNCTIONS];
	double rm1, rp1;
	int i, j, k, fun, pej, nodeI, nodeJ, nodeK, n, elpmax = 0, fun0;
	SRedge* edge;
	double N[10];

	//corner and midedge functions:
	model.map.TetVolumeCoords(r, s, t, lv);
	if (calltype != derivonly)
	{
		elem->ShapeFunctions(r, s, t, N);
		for (i = 0; i < 10; i++)
			basisvec[i] = N[i];

	}
	if(calltype != basisonly)
	{
		elem->FillJacobian(r, s, t, false);
		for (i = 0; i < 10; i++)
		{
			dbasisdr[i] = elem->dNdr[i];
			dbasisds[i] = elem->dNds[i];
			dbasisdt[i] = elem->dNdt[i];
		}
	}


	fun = 10;
	//higher edge functions:
	int lej, direction, n1, n2;

	for (lej = 0; lej < 6; lej++)
	{
		edge = elem->GetEdge(lej);
		pej = edge->pOrder;
		if(pej > elpmax)
			elpmax = pej;

		direction = elem->GetLocalEdgeDirection(lej);

		n1 = model.map.tetEdgeLocalNodes[lej][0];
		n2 = model.map.tetEdgeLocalNodes[lej][1];
		if(direction > 0)
		{
			nodeI = n1;
			nodeJ = n2;
		}
		else
		{
			nodeI = n2;
			nodeJ = n1;
		}

		LI = lv[nodeI];
		LJ = lv[nodeJ];
		rej = LJ - LI;

		//edge basis functions go to zero at corners. catch this case explicitly
		if (((1.0 - rej) < RELSMALL) || ((rej + 1.0) < RELSMALL))
		{
			fun0 = fun;
			if(calltype != derivonly)
			{
				for (i = 3; i <= pej; i++)
				{
					basisvec[fun] = 0.0;
					fun++;
				}
			}
			if(calltype != basisonly)
			{
				//asymptotic branch not supported for derivs
				fun = fun0;
				for (i = 3; i <= pej; i++)
				{
					dbasisdr[fun] = 0.0;
					dbasisds[fun] = 0.0;
					dbasisdt[fun] = 0.0;
					fun++;
				}

			}
			continue;
		}

		if (calltype == basisonly)
			Basisfuns1d(pej, rej, basisej);
		else
		{
			dLJdr = model.map.dLtdr[nodeJ];
			dLJds = model.map.dLtds[nodeJ];
			dLJdt = model.map.dLtdt[nodeJ];
			dLIdr = model.map.dLtdr[nodeI];
			dLIds = model.map.dLtds[nodeI];
			dLIdt = model.map.dLtdt[nodeI];
			Basisfuns1d(pej, rej, basisej, dbasisdrej, true);
		}

		double r2f = 4.0 / (1.0 - rej*rej);
		double llr24 = LI*LJ*r2f;
		fun0 = fun;
		if (calltype != derivonly)
		{
			for (i = 3; i <= pej; i++)
			{
				basisvec[fun] = llr24*basisej[i];
				fun++;
			}
		}
		if (calltype != basisonly)
		{
			fun = fun0;
			double drejdr, drejds, drejdt, tmp, tmp2;
			drejdr = dLJdr - dLIdr;
			drejds = dLJds - dLIds;
			drejdt = dLJdt - dLIdt;
			double halfrer2f = 0.5*rej*r2f;
			for (i = 3; i <= pej; i++)
			{
				basisejt = basisej[i];
				tmp = llr24*(dbasisdrej[i] + halfrer2f*basisejt);
				tmp2 = r2f*basisejt;
				dbasisdr[fun] = (dLIdr*LJ + dLJdr*LI)*tmp2 + drejdr*tmp;
				dbasisds[fun] = (dLIds*LJ + dLJds*LI)*tmp2 + drejds*tmp;
				dbasisdt[fun] = (dLIdt*LJ + dLJdt*LI)*tmp2 + drejdt*tmp;
				fun++;
			}
		}
	}

	//face functions:
	int lface, gface, pmax, nflI, nflJ, nflK;
	double rf11, sf11;
	for (lface = 0; lface < 4; lface++)
	{
		gface = elem->GetLocalFaceGlobalId(lface);
		pmax = FaceGetpmax(gface);
		if(pmax < 3)
			continue;
		nflI = elem->GetLocalFaceGlobalNodeOrder(lface, 0);
		nflJ = elem->GetLocalFaceGlobalNodeOrder(lface, 1);
		nflK = elem->GetLocalFaceGlobalNodeOrder(lface, 2);
		nodeI = model.map.tetFaceLocalNodes[lface][nflI];
		nodeJ = model.map.tetFaceLocalNodes[lface][nflJ];
		nodeK = model.map.tetFaceLocalNodes[lface][nflK];
		LI = lv[nodeI];
		LJ = lv[nodeJ];
		LK = lv[nodeK];
		//local coordinates on the face the run from -1 to 1:
		rf11 = LJ - LI;
		sf11 = 2.0*LK - 1.0;
		double lll = LI*LJ*LK;
		n = pmax - 3;
		LegendrePns(n, rf11, Pr);
		LegendrePns(n, sf11, Ps);
		double ht;
		fun0 = fun;
		if(calltype != derivonly)
		{
			for(i = 0; i <= n; i++)
			{
				ht = lll*Pr[i];
				for (j = 0; j <= n; j++)
				{
					if ((i + j) > n)
						break;
					basisvec[fun] = ht*Ps[j];
					fun++;
				}
			}
		}
		if(calltype != basisonly)
		{
			double t1r, t1s, t1t, t2r, t2s, t2t, t3r, t3s, t3t, dpipj, dpjpi;
			dLIdr = model.map.dLtdr[nodeI];
			dLIds = model.map.dLtds[nodeI];
			dLIdt = model.map.dLtdt[nodeI];
			dLJdr = model.map.dLtdr[nodeJ];
			dLJds = model.map.dLtds[nodeJ];
			dLJdt = model.map.dLtdt[nodeJ];
			dLKdr = model.map.dLtdr[nodeK];
			dLKds = model.map.dLtds[nodeK];
			dLKdt = model.map.dLtdt[nodeK];
			double LIJ = LI*LJ;
			double LIK = LI*LK;
			double LJK = LJ*LK;
			double lll2 = 2.0*lll;
			t1r = dLIdr*LJK + dLJdr*LIK + dLKdr*LIJ;
			t1s = dLIds*LJK + dLJds*LIK + dLKds*LIJ;
			t1t = dLIdt*LJK + dLJdt*LIK + dLKdt*LIJ;
			t2r = lll*(dLJdr - dLIdr);
			t2s = lll*(dLJds - dLIds);
			t2t = lll*(dLJdt - dLIdt);
			t3r = lll2*dLKdr;
			t3s = lll2*dLKds;
			t3t = lll2*dLKdt;
			fun = fun0;
			LegendrePnDerivs(n, Pr, dpdr);
			LegendrePnDerivs(n, Ps, dpds);
			for (i = 0; i <= n; i++)
			{
                double pri = Pr[i];
                double dpdri = dpdr[i];
				for(j = 0; j <= n; j++)
				{
					if((i+j) > n)
						break;
					pipj = pri*Ps[j];
					dpipj = dpdri*Ps[j];
					dpjpi = dpds[j]*pri;
					dbasisdr[fun] = t1r*pipj + t2r*dpipj + t3r*dpjpi;
					dbasisds[fun] = t1s*pipj + t2s*dpipj + t3s*dpjpi;
					dbasisdt[fun] = t1t*pipj + t2t*dpipj + t3t*dpjpi;
					fun++;
				}
			}
		}
	}

	//volume functions
	if(elpmax < 4)
		return fun;
	n = elpmax - 4;
	double llll, L1, L2, L3, L4, pi, pj;
	L1 = lv[0];
	L2 = lv[1];
	L3 = lv[2];
	L4 = lv[3];
	llll = L1*L2*L3*L4;
	LegendrePns(n, r, Pr);
	LegendrePns(n, s, Ps);
	LegendrePns(n, t, Pt);
	fun0 = fun;
	if(calltype != derivonly)
	{
		for (i = 0; i <= n; i++)
		{
			pi = llll*Pr[i];
			for(j = 0; j <= n; j++)
			{
				pipj = pi*Ps[j];
				for(k = 0; k <= n; k++)
				{
					if ((i + j + k) > n)
						break;
					basisvec[fun] = pipj*Pt[k];
					fun++;
				}
			}
		}
	}
	if(calltype != basisonly)
	{
		fun = fun0;
		LegendrePnDerivs(n, Pr, dpdr);
		LegendrePnDerivs(n, Ps, dpds);
		LegendrePnDerivs(n, Pt, dpdt);
		double tmpr, tmps, tmpt, ppp;
		double lllldri, lllldsj;
		double L234 = L2*L3*L4;
		double L134 = L1*L3*L4;
		double L124 = L1*L2*L4;
		double L123 = L1*L2*L3;
		tmpr = (DLT1DR*L234 + DLT2DR*L134 + DLT3DR*L124 + DLT4DR*L123);
		tmps = (DLT1DS*L234 + DLT2DS*L134 + DLT3DS*L124 + DLT4DS*L123);
		tmpt = (DLT1DT*L234 + DLT2DT*L134 + DLT3DT*L124 + DLT4DT*L123);
		for (i = 0; i <= n; i++)
		{
			pi = Pr[i];
			lllldri = llll*dpdr[i];
			for(j = 0; j <= n; j++)
			{
				pj = Ps[j];
				pipj = pi*pj;
				lllldsj = llll*dpds[j];
				for (k = 0; k <= n; k++)
				{
					if ((i + j + k) > n)
						break;
					ppp = pipj*Pt[k];
					dbasisdr[fun] = tmpr*ppp + lllldri * pj*Pt[k];
					dbasisds[fun] = tmps*ppp + lllldsj * pi*Pt[k];
					dbasisdt[fun] = tmpt*ppp + llll*dpdt[k] * pipj;
					fun++;
				}
			}
		}
	}

	return fun;
}

int SRbasis::FaceGetpmax(int gface)
{
	//get max polynomial order for global face
	//input:
		//gface =global face id
	//return:
		//max polynomial order
	SRface* face = model.GetFace(gface);
	return FaceGetpmax(face);
}

void SRbasis::LegendrePnDerivs(int n,double *P, double *dpdr)
{
	//derivatives of Legendre polynomials up to order n
	//input:
		//n = highest polynomial order
		//P = vector of Legendre polynomials
	//output:
		//dpdr = derivatives of Legendre polynomials at r

	//derivatives of legendre pns, see Szabo p 355
	dpdr[0] = 0.0;
	dpdr[1] = 1.0;
	for (int i = 2; i <= n; i++)
		dpdr[i] = (2.0*i - 1)*P[i - 1] + dpdr[i - 2];
}

int SRbasis::CountTetFunctions(SRelement *elem, int &nint)
{
	//count basis functions for a tet element
	//input:
		//elem = pointer to element
	//output:
		//nint = number of internal basis functions in element
	//return:
		//total number of basis functions in element

	int nfun, pej, elpmax = 0;
	//corner functions:
	nfun = 4;
	//edge functions:
	int lej;
	SRedge *edge;
	for (lej = 0; lej < 6; lej++)
	{
		edge = elem->GetEdge(lej);
		pej = edge->pOrder;
		if(pej > elpmax)
			elpmax = pej;
		nfun += (pej-1);
	}

	//face functions:
	SRface *face;
	for(int lface = 0; lface < 4; lface++)
	{
		face = elem->GetFace(lface);
		nfun += CountFaceInternalFunctions(face);
	}

	//volume functions
	//szabo, p 244
	if(elpmax < 4)
		nint = 0;
	else
		nint = (elpmax - 1)*(elpmax - 2)*(elpmax - 3) / 6;
	nfun += nint;
	return nfun;
}

int SRbasis::CountElementFunctions(SRelement *elem, int &nint)
{
	//count basis functions for an element
	//input:
		//elem = pointer to element
	//output:
		//nint = number of internal basis functions in element
	//return:
		//total number of basis functions in element

	if (elem->type == tet)
		return CountTetFunctions(elem, nint);
	else if (elem->type == wedge)
		return CountWedgeFunctions(elem, nint);
	else 
		return CountBrickFunctions(elem, nint);
}


int SRbasis::CountFaceTotalFunctions(SRface *face, int &nint)
{
	//count basis functions for a global face
	//input:
		//face = pointer to global face
	//output:
		//nint = number of internal basis functions on face
	//return:
		//total number of basis functions on face

	int i, nej, nfun, pej, pmax = 0;
	nej = face->GetNumLocalEdges();
	//number of nodes=number of edges.
	//nodal functions:
	nfun = nej;
	//edge functions:
	SRedge* edge;
	for(i = 0; i < nej; i++)
	{
		edge = face->GetEdge(i);
		pej = edge->pOrder;
		if(pej > pmax)
			pmax = pej;
		nfun += (pej - 1);
	}
	//internal face functions:
	nint = CountFaceInternalFunctions(nej, pmax);
	nfun += nint;
	return nfun;
}

int SRbasis::CountFaceInternalFunctions(int numsides, int pmax)
{
	//count number of internal face functions on a face given
	//the max p order of any of its edges
	//input:
		//numsides = number of sides of the face- 3 or 4
		//pmax = max p order of any of the face's edges
	//return:
		//number of internal face functions
	int nint = 0;
	if(numsides == 3)
	{
		//szabo, p103:
		if(pmax > 2)
			nint = (pmax - 2)*(pmax - 1) / 2;
		else
			nint = 0;
	}
	else
	{
		//szabo, p99:
		if(pmax > 3)
			nint = (pmax - 2)*(pmax - 3) / 2;
		else
			nint = 0;
	}
	return nint;
}

int SRbasis::CountFaceInternalFunctions(SRface *face)
{
	//count number of internal face functions on a face
	//input:
		//face = pointer to face
	//return:
		//number of internal face functions


	int nej, pmax;
	nej = face->GetNumLocalEdges();
	pmax = FaceGetpmax(face);
	return CountFaceInternalFunctions(nej, pmax);
}

int SRbasis::CountFaceInternalFunctions(SRelement *elem, int lface)
{
	//count number of internal face functions on a local face of an element
	//input:
		// elem = pointer to element
		// lface = local face number
	//return:
		// number of internal face functions

	SRface* gface = elem->GetFace(lface);
	return CountFaceInternalFunctions(gface);
}

int SRbasis::FaceGetpmax(SRface *face)
{
	//determine max p order on a face
	//input:
		//face = pointer to face
	//return:
		// max p order
	SRedge* edge;
	int p, pmax = 0;
	int i;
	for (i = 0; i < face->GetNumLocalEdges() ; i++)
	{
		edge = face->GetEdge(i);
		p = edge->pOrder;
		if(p > pmax)
			pmax = p;
	}
	return pmax;
}


int SRbasis::CountFaceTotalFunctions(SRelement *elem, int lface, int &nint)
{
	//count total number of face functions on a local face of an element
	//input:
		// elem = pointer to element
		// lface = local face number
	//output:
		// nint = number of internal face functions
	//return:
		// number of face functions
	SRface* gface = elem->GetFace(lface);
	return CountFaceTotalFunctions(gface, nint);
}
