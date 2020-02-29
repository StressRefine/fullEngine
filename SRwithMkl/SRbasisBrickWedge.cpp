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
// SRbasisBrickWedge.cpp: implementation of the SRbasis class- bricks and wedges
//
//////////////////////////////////////////////////////////////////////

#include "SRmodel.h"
#include "SRbasis.h"

extern SRmodel model;

int SRbasis::QuadFaceBasisFuncs(double r, double s, int *directionv, int *pejv, double *basisvec, SRface* face)
{
	//basis functions for a quad face
	//input:
		//r,s = natural coordinates on face
		//directionv = direction for each edge (-1 or +1)
		//pejv = polynomial orders for each edge
	//output:
		//basisvec = vector of basis functions
	//return:
		//total number of basis functions
	int fun;
	int pmax = 0;

	//corner and midedge functions:
	double N[8];
	model.map.FaceShapeFunctions(face, r, s, N);
	for (int n = 0; n < 8; n++)
		basisvec[n] = N[n];

	fun = 8;
	//higher edge functions:
	for(int l = 0; l < 4; l++)
	{
		QuadEdgeBlend(l, directionv[l], r, s, fun, pejv[l], basisvec);
		if (pejv[l] > pmax)
			pmax = pejv[l];
	}

	//bubble funcs:
	QuadBubbleFuncs(pmax, r, s, fun, basisvec);
	return fun;
}

void SRbasis::QuadEdgeBlend(int localedge, int direction, double r, double s, int &fun, int pej, double *basisvec)
{
	//blend edge basis functions to form element functions for quadrilateral
	//input:
	//localedge = local edge number
	//direction = 1 if local edge direction matches global direction else -1
		//r,s = natural coordinates on face
		//fun = function number to use to store next element function
		//pej = polynomial order for this edge
	//output:
		//fun: updated
		//basisvec = element shape function vector

	double rej, blend;
	double basisej[MAX1DFUNCTIONS];
	int i;
	double directionfactor = (double)direction;

	if (localedge == 0)
	{
		rej = r;
		blend = 0.5*(1.0 - s);
	}
	else if (localedge == 1)
	{
		rej = s;
		blend = 0.5*(1.0 + r);
	}
	else if (localedge == 2)
	{
		rej = r;
		blend = 0.5*(1.0 + s);
	}
	else if (localedge == 3)
	{
		rej = s;
		blend = 0.5*(1.0 - r);
	}
	else
		ERROREXIT;

	rej *= directionfactor;
	Basisfuns1d(pej, rej, basisej);
	for (i = 3; i <= pej; i++)
	{
		basisvec[fun] = blend*basisej[i];
		fun++;
	}
}

void SRbasis::QuadBubbleFuncs(int pmax, double r, double s, int &fun, double *basisvec)
{
	//bubble basis functions for quadrilateral
	//input:
		//pmax = max polynomial order on any edges
		//r,s = natural coordinates on element face
		//fun = function number to use to store next element function
	//output:
		//fun: updated
		//basisvec = element shape function vector

	double basis1dr[MAX1DFUNCTIONS], basis1ds[MAX1DFUNCTIONS];
	int n;
	if(pmax < 4)
		return;
	n = pmax - 2;

	Basisfuns1d(n, r, basis1dr);
	Basisfuns1d(n, s, basis1ds);

	int i,j;
	for (i = 2; i <= n; i++)
	{
		for(j = 2; j <= n; j++)
		{
			if((i+j) > pmax)
				break;
			basisvec[fun] = basis1dr[i] * basis1ds[j];
			fun++;
		}
	}
}

int SRbasis::BrickBasisFuncs(double r, double s, double t, SRelement *elem, double *basisvec, double *dbasisdr, double *dbasisds, double *dbasisdt, SRBasisCallType calltype)
{
	//basis functions and their derivatives at r,s,t for a brick element
	//input:
		//r,s,t = natural coordinates in the elements
		//calltype = basisonly if don't need derivatives, derivonly if don't need basis function, or both
	//output:
		//basisvec = vector of basis functions
		//dbasisdr,dbasisds,dbasisdt = vector of basis function derivatives
	//return:
		//total number of basis functions

	double rej, blend;
	double phiI[MAX1DFUNCTIONS], phiJ[MAX1DFUNCTIONS], phiK[MAX1DFUNCTIONS];
	double phiprimeI[MAX1DFUNCTIONS], phiprimeJ[MAX1DFUNCTIONS], phiprimeK[MAX1DFUNCTIONS];
	double r0, s0, t0, ri, si, ti;
	int i, j, k, fun, pej, n, fun0;
	double N[20];

	//corner and midedge functions:
	if (calltype != derivonly)
	{
		elem->ShapeFunctions(r, s, t, N);
		for (i = 0; i < 20; i++)
			basisvec[i] = N[i];
	}
	if (calltype != basisonly)
	{
		elem->FillJacobian(r, s, t, false);
		for (i = 0; i < 20; i++)
		{
			dbasisdr[i] = elem->dNdr[i];
			dbasisds[i] = elem->dNds[i];
			dbasisdt[i] = elem->dNdt[i];
		}
	}

	fun = 20;
	//higher edge functions:
	int lej;
	SRedge* edge;
	double directionfactor;
	int elpmax = 0;
	int corner;
	for(lej = 0; lej < 12; lej++)
	{
		edge = elem->GetEdge(lej);
		pej = edge->pOrder;
		if (pej > elpmax)
			elpmax = pej;
		if (pej < 2)
			continue;

		directionfactor = (double) elem->GetLocalEdgeDirection(lej);

		if (lej == 0 || lej == 3)
			corner = 0;
		else if (lej == 1)
			corner = 1;
		else if (lej == 2)
			corner = 3;
		else if (lej == 4 || lej == 7)
			corner = 4;
		else if (lej == 5)
			corner = 5;
		else if (lej == 6)
			corner = 7;
		else if (lej == 8 || lej == 9 || lej == 10 || lej == 11)
			corner = lej - 8;

		SRvec3& rst = model.map.brickNodes[corner];
		if (lej == 0 || lej == 2 || lej == 4 || lej == 6)
		{
			//parallel to r:
			rej = r*directionfactor;
			si = rst.d[1];
			ti = rst.d[2];
			s0 = si*s;
			t0 = ti*t;
		}
		else if (lej == 1 || lej == 3 || lej == 5 || lej == 7)
		{
			//parallel to s:
			rej = s*directionfactor;
			ri = rst.d[0];
			ti = rst.d[2];
			r0 = ri*r;
			t0 = ti*t;
		}
		else if (lej == 8 || lej == 9 || lej == 10 || lej == 11)
		{
			//parallel to t:
			rej = t*directionfactor;
			ri = rst.d[0];
			si = rst.d[1];
			r0 = ri*r;
			s0 = si*s;
		}

		if (calltype == basisonly)
			Basisfuns1d(pej, rej, phiI);
		else
			Basisfuns1d(pej, rej, phiI, phiprimeI, true);

		fun0 = fun;
		if(calltype != derivonly)
		{
			if (lej == 0 || lej == 2 || lej == 4 || lej == 6)
			{
				for (i = 3; i <= pej; i++)
				{
					basisvec[fun] = 0.25*phiI[i] * (1.0 + s0)*(1.0 + t0);
					fun++;
				}
			}
			else if (lej == 1 || lej == 3 || lej == 5 || lej == 7)
			{
				for (i = 3; i <= pej; i++)
				{
					basisvec[fun] = 0.25*phiI[i] * (1.0 + r0)*(1.0 + t0);
					fun++;
				}
			}
			else if (lej == 8 || lej == 9 || lej == 10 || lej == 11)
			{
				for(i = 3; i <= pej; i++)
				{
					basisvec[fun] = 0.25*phiI[i] * (1.0 + r0)*(1.0 + s0);
					fun++;
				}
			}
		}

		if(calltype != basisonly)
		{
			fun = fun0;
			if (lej == 0 || lej == 2 || lej == 4 || lej == 6)
			{
				//parallel to r:
				for (i = 3; i <= pej; i++)
				{
					dbasisdr[fun] = 0.25*directionfactor*phiprimeI[i] * (1.0 + s0)*(1.0 + t0);
					dbasisds[fun] = si*0.25*phiI[i] * (1.0 + t0);
					dbasisdt[fun] = ti*0.25*phiI[i] * (1.0 + s0);
					fun++;
				}
			}
			else if (lej == 1 || lej == 3 || lej == 5 || lej == 7)
			{
				//parallel to s:
				for (i = 3; i <= pej; i++)
				{
					dbasisds[fun] = 0.25*directionfactor*phiprimeI[i] * (1.0 + r0)*(1.0 + t0);
					dbasisdr[fun] = ri*0.25*phiI[i] * (1.0 + t0);
					dbasisdt[fun] = ti*0.25*phiI[i] * (1.0 + r0);
					fun++;
				}
			}
			else if (lej == 8 || lej == 9 || lej == 10 || lej == 11)
			{
				//parallel to t:
				for(i = 3; i <= pej; i++)
				{
					dbasisdt[fun] = 0.25*directionfactor*phiprimeI[i] * (1.0 + r0)*(1.0 + s0);
					dbasisdr[fun] = ri*0.25*phiI[i] * (1.0 + s0);
					dbasisds[fun] = si*0.25*phiI[i] * (1.0 + r0);
					fun++;
				}
			}
		}
	}

	//face functions
	int lface, gface, pmax;
	double rfl, sfl, rf, sf, drfdrfl, drfdsfl, dsfdrfl, dsfdsfl, drfldr, drflds, drfldt;
	double dsfldr, dsflds, dsfldt, drfdr, drfds, drfdt, dsfdr, dsfds, dsfdt, dblend;
	for (lface = 0; lface < 6; lface++)
	{
		gface = elem->GetLocalFaceGlobalId(lface);
		pmax = FaceGetpmax(gface);
		if(pmax < 4)
			continue;
		pej = pmax - 2;
		drfldr = 0.0;
		drflds = 0.0;
		drfldt = 0.0;
		dsfldr = 0.0;
		dsflds = 0.0;
		dsfldt = 0.0;

		if (lface == 0 || lface == 1)
		{
			rfl = r;
			sfl = s;
			drfldr = 1.0;
			dsflds = 1.0;
		}
		else if (lface == 2 || lface == 3)
		{
			rfl = s;
			sfl = t;
			drflds = 1.0;
			dsfldt = 1.0;
		}
		else if (lface == 4 || lface == 5)
		{
			rfl = r;
			sfl = t;
			drfldr = 1.0;
			dsfldt = 1.0;
		}

		if (lface == 0)
		{
			blend = 0.5*(1.0 - t);
			dblend = -0.5;
		}
		else if (lface == 1)
		{
			blend = 0.5*(1.0 + t);
			dblend = 0.5;
		}
		else if (lface == 2)
		{
			blend = 0.5*(1.0 - r);
			dblend = -0.5;
		}
		else if (lface == 3)
		{
			blend = 0.5*(1.0 + r);
			dblend = 0.5;
		}
		else if (lface == 4)
		{
			blend = 0.5*(1.0 - s);
			dblend = -0.5;
		}
		else if (lface == 5)
		{
			blend = 0.5*(1.0 + s);
			dblend = 0.5;
		}

		elem->GetQuadFaceGlobalCoords(lface, rfl, sfl, rf, sf);
		if (calltype == basisonly)
		{
			Basisfuns1d(pej, rf, phiI);
			Basisfuns1d(pej, sf, phiJ);
		}
		else
		{
			Basisfuns1d(pej, rf, phiI, phiprimeI, true);
			Basisfuns1d(pej, sf, phiJ, phiprimeJ, true);
		}
		double pipj;
		fun0 = fun;
		if(calltype != derivonly)
		{
			for(i = 2; i <= pej; i++)
			{
				double pi = phiI[i];
				for(j = 2; j <= pej; j++)
				{
					if ((i + j) <= pmax)
					{
						pipj = pi * phiJ[j];
						basisvec[fun] = pipj*blend;
						fun++;
					}
				}
			}
		}
		if(calltype != basisonly)
		{
			fun = fun0;
			elem->GetQuadFaceDerivs(lface, rfl, sfl, drfdrfl, drfdsfl, dsfdrfl, dsfdsfl);
			drfdr = drfdrfl*drfldr + drfdsfl*dsfldr;
			drfds = drfdrfl*drflds + drfdsfl*dsflds;
			drfdt = drfdrfl*drfldt + drfdsfl*dsfldt;
			dsfdr = dsfdrfl*drfldr + dsfdsfl*dsfldr;
			dsfds = dsfdrfl*drflds + dsfdsfl*dsflds;
			dsfdt = dsfdrfl*drfldt + dsfdsfl*dsfldt;
			if (lface == 0 || lface == 1)
			{
				for (i = 2; i <= pej; i++)
				{
					double pi = phiI[i];
					double dpi = phiprimeI[i];
					double drfdrdpi = drfdr*dpi;
					double drfdsdpi = drfds*dpi;
					double dsfdrpi = dsfdr*pi;
					double dsfdspi = dsfds*pi;
					for (j = 2; j <= pej; j++)
					{
						if ((i + j) <= pmax)
						{
							pipj = pi * phiJ[j];
							dbasisdr[fun] = blend*(drfdrdpi * phiJ[j] + dsfdrpi * phiprimeJ[j]);
							dbasisds[fun] = blend*(drfdsdpi * phiJ[j] + dsfdspi * phiprimeJ[j]);
							dbasisdt[fun] = dblend*pipj;
							fun++;
						}
					}
				}
			}
			else if (lface == 2 || lface == 3)
			{
				for (i = 2; i <= pej; i++)
				{
					double pi = phiI[i];
					double dpi = phiprimeI[i];
					double drfdsdpi = drfds*dpi;
					double drfdtdpi = drfdt*dpi;
					double dsfdspi = dsfds*pi;
					double dsfdtpi = dsfdt*pi;
					for (j = 2; j <= pej; j++)
					{
						if ((i + j) <= pmax)
						{
							pipj = phiI[i] * phiJ[j];
							dbasisds[fun] = blend*(drfdsdpi * phiJ[j] + dsfdspi * phiprimeJ[j]);
							dbasisdt[fun] = blend*(drfdtdpi * phiJ[j] + dsfdtpi * phiprimeJ[j]);
							dbasisdr[fun] = dblend*pipj;
							fun++;
						}
					}
				}
			}
			else if (lface == 4 || lface == 5)
			{
				for(i = 2; i <= pej; i++)
				{
					double pi = phiI[i];
					double dpi = phiprimeI[i];
					double drfdrdpi = drfdr*dpi;
					double drfdtdpi = drfdt*dpi;
					double dsfdrpi = dsfdr*pi;
					double dsfdtpi = dsfdt*pi;
					for(j = 2; j <= pej; j++)
					{
						if((i+j) <= pmax)
						{
							pipj = phiI[i] * phiJ[j];
							dbasisdr[fun] = blend*(drfdrdpi * phiJ[j] + dsfdrpi * phiprimeJ[j]);
							dbasisdt[fun] = blend*(drfdtdpi * phiJ[j] + dsfdtpi * phiprimeJ[j]);
							dbasisds[fun] = dblend*pipj;
							fun++;
						}
					}
				}
			}
		}
	}

//volume functions:
	if(elpmax < 6)
		return fun;
	n = elpmax - 4;
	double pipj;
	if(calltype == basisonly)
	{
		Basisfuns1d(n, r, phiI);
		Basisfuns1d(n, s, phiJ);
		Basisfuns1d(n, t, phiK);
	}
	else
	{
		Basisfuns1d(n, r, phiI, phiprimeI, true);
		Basisfuns1d(n, s, phiJ, phiprimeJ, true);
		Basisfuns1d(n, t, phiK, phiprimeK, true);
	}
	fun0 = fun;
	if(calltype != derivonly)
	{
		for(i = 2; i <= n; i++)
		{
			double pi = phiI[i];
			for (j = 2; j <= n; j++)
			{
				pipj = pi * phiJ[j];
				for (k = 2; k <= n; k++)
				{
					if ((i + j + k) > elpmax)
						break;
					basisvec[fun] = pipj*phiK[k];
					fun++;
				}
			}
		}
	}
	if(calltype != basisonly)
	{
		fun = fun0;
		double pidpj, pjdpi;
		for (i = 2; i <= n; i++)
		{
			double pi = phiI[i];
			double dpi = phiprimeI[i];
			for (j = 2; j <= n; j++)
			{
				pipj = pi * phiJ[j];
				pjdpi = dpi * phiJ[j];
				pidpj = phiprimeJ[j] * pi;
				for (k = 2; k <= n; k++)
				{
					if ((i + j + k) > elpmax)
						break;
					dbasisdr[fun] = pjdpi*phiK[k];
					dbasisds[fun] = pidpj*phiK[k];
					dbasisdt[fun] = pipj*phiprimeK[k];
					fun++;
				}
			}
		}
	}

	return fun;
}

int SRbasis::CountBrickFunctions(SRelement *elem, int &nint)
{
	//count basis functions for a brick element
	//input:
		//elem = pointer to element
	//output:
		//nint = number of internal basis functions in element
	//return:
		//total number of basis functions in element
	int nfun,pej,n;

	//corner functions:
	nfun = 8;
	//edge functions:
	int lej;
	SRedge* ej;
	int elpmax = 0;
	for(lej = 0; lej < 12; lej++)
	{
		ej = elem->GetEdge(lej);
		pej = ej->pOrder;
		if(pej > elpmax)
			elpmax = pej;
		n = pej - 1;
		nfun += n;
	}

	//face functions
	int lface;
	SRface* face;
	for(lface = 0; lface < 6; lface++)
	{
		face = elem->GetFace(lface);
		n = CountFaceInternalFunctions(face);
		nfun += n;
	}

//volume functions:
	//from Szabo, p241:
	if(elpmax < 6)
		nint = 0;
	else
		nint = (elpmax - 3)*(elpmax - 4)*(elpmax - 5) / 6;
	nfun += nint;
	return nfun;
}

int SRbasis::WedgeBasisFuncs(double r, double s, double t, SRelement *elem, double *basisvec, double *dbasisdr, double *dbasisds, double *dbasisdt, SRBasisCallType calltype)

{
	//basis functions and their derivatives at r,s,t for a wedge element
	//input:
		//r,s,t = natural coordinates in the elements
		//calltype = basisonly if don't need derivatives, derivonly if don't need basis function, or both
	//output:
		//basisvec = vector of basis functions
		//dbasisdr,dbasisds,dbasisdt = vector of basis function derivatives
	//return:
		//total number of basis functions

	double tm1 = 0.5*(1.0 - t);
	double tp1 = 0.5*(1.0 + t);
	int fun, fun0;
    int elpmax = 0;
	double L1, L2, L3, tmp, tmpt;
	double llr2f, re2, drejdr, drejds;
	double N[15];
	model.map.TriangleAreaCoords(r, s, L1, L2, L3);
	//corner and midedge functions:
	if (calltype != derivonly)
	{
		elem->ShapeFunctions(r, s, t, N);
		for (int i = 0; i < 15; i++)
			basisvec[i] = N[i];

	}
	if (calltype != basisonly)
	{
		elem->FillJacobian(r, s, t, false);
		for (int i = 0; i < 15; i++)
		{
			dbasisdr[i] = elem->dNdr[i];
			dbasisds[i] = elem->dNds[i];
			dbasisdt[i] = elem->dNdt[i];
		}
	}

	fun = 15;
	//higher edge functions for tri faces:
	double basisej[MAX1DFUNCTIONS], dbasisdrej[MAX1DFUNCTIONS], lv[3], dldr[3], dlds[3];
	double rej, tblend, r2f, basisejt, LI, LJ, dLIdr, dLJdr, dLIds, dLJds, dtblenddt;
	double tmpr, tmps;
	int i, nodeI, nodeJ, nodeK, lej, direction, pej;
	SRedge* edge;
	static int triLedge[6] = { 0, 1, 2, 0, 1, 2 };
	lv[0] = L1;
	lv[1] = L2;
	lv[2] = L3;
	dldr[0] = DL1DR;
	dldr[1] = DL2DR;
	dldr[2] = DL3DR;
	dlds[0] = DL1DS;
	dlds[1] = DL2DS;
	dlds[2] = DL3DS;

	for(lej = 0; lej < 6; lej++)
	{
		edge = elem->GetEdge(lej);
		pej = edge->pOrder;
		if (pej > elpmax)
			elpmax = pej;

		direction = elem->GetLocalEdgeDirection(lej);

		int triLej = triLedge[lej];
		if (direction > 0)
		{
			nodeI = model.map.triEdgeLocalNodes[triLej][0];
			nodeJ = model.map.triEdgeLocalNodes[triLej][1];
		}					
		else				
		{					
			nodeI = model.map.triEdgeLocalNodes[triLej][1];
			nodeJ = model.map.triEdgeLocalNodes[triLej][0];
		}
		if(lej < 3)
		{
			tblend = tm1;
			dtblenddt = -0.5;
		}
		else
		{
			tblend = tp1;
			dtblenddt = 0.5;
		}
		
		LJ = lv[nodeJ];
		LI = lv[nodeI];
		dLJdr = dldr[nodeJ];
		dLJds = dlds[nodeJ];
		dLIdr = dldr[nodeI];
		dLIds = dlds[nodeI];
		rej = LJ - LI;

		//edge basis functions go to zero at corners. catch this case explicitly
		if (((1.0 - rej) < RELSMALL) || ((rej + 1.0) < RELSMALL))
		{
			fun0 = fun;
			if (calltype != derivonly)
			{
				for(i = 3; i <= pej; i++)
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
			Basisfuns1d(pej, rej, basisej, dbasisdrej, true);

		double rm1 = 1.0 / (1.0 - rej);
		double rp1 = 1.0 / (1.0 + rej);
		r2f = 4.0*rm1*rp1;
		fun0 = fun;
		tmp = LI*LJ*r2f*tblend;
		if (calltype != derivonly)
		{
			for (i = 3; i <= pej; i++)
			{
				basisvec[fun] = tmp*basisej[i];
				fun++;
			}
		}
		if (calltype != basisonly)
		{
			fun = fun0;
			re2 = 0.5*rej*r2f;
			double tmprs;
			llr2f = LI*LJ*r2f;
			drejdr = dLJdr - dLIdr;
			drejds = dLJds - dLIds;
			tmpr = (dLIdr*LJ + dLJdr*LI)*r2f;
			tmps = (dLIds*LJ + dLJds*LI)*r2f;
			tmpt = dtblenddt*llr2f;
			for (i = 3; i <= pej; i++)
			{
				basisejt = basisej[i];
				tmprs = llr2f*(re2*basisejt + dbasisdrej[i]);
				dbasisdr[fun] = (tmpr*basisejt + drejdr*tmprs)*tblend;
				dbasisds[fun] = (tmps*basisejt + drejds*tmprs)*tblend;
				dbasisdt[fun] = tmpt*basisejt;
				fun++;
			}
		}
	}

	//higher edge functions for quad faces:
	for (lej = 6; lej < 9; lej++)
	{
		edge = elem->GetEdge(lej);
		pej = edge->pOrder;
		if (pej > elpmax)
			elpmax = pej;

		double directionfactor;
		nodeI = lej - 6;
		LI = lv[nodeI];
		dLIdr = dldr[nodeI];
		dLIds = dlds[nodeI];
		direction = elem->GetLocalEdgeDirection(lej);
		directionfactor = (double)direction;
		rej = t*directionfactor;

		if (calltype == basisonly)
			Basisfuns1d(pej, rej, basisej);
		else
			Basisfuns1d(pej, rej, basisej, dbasisdrej, true);

		fun0 = fun;
		if(calltype != derivonly)
		{
			for(i = 3; i <= pej; i++)
			{
				basisvec[fun] = basisej[i] * LI;
				fun++;
			}
		}
		if(calltype != basisonly)
		{
			fun = fun0;
			tmp = directionfactor*LI;
			for (i = 3; i <= pej; i++)
			{
				dbasisdr[fun] = dLIdr*basisej[i];
				dbasisds[fun] = dLIds*basisej[i];
				dbasisdt[fun] = tmp*dbasisdrej[i];
				fun++;
			}
		}
	}

	//face functions:

	int lface, gface, pmax, n, j;
	double Pr[MAX1DFUNCTIONS], Ps[MAX1DFUNCTIONS], dpdr[MAX1DFUNCTIONS], dpds[MAX1DFUNCTIONS];
	double phiI[MAX1DFUNCTIONS], phiprimeI[MAX1DFUNCTIONS], phiJ[MAX1DFUNCTIONS], phiprimeJ[MAX1DFUNCTIONS];
	double LK, dLKdr, dLKds, lll, tmpr1, tmpr2, tmpr3, tmps1, tmps2, tmps3, pipj, rf, sf;

	//tri faces:
	for(lface = 0; lface < 2; lface++)
	{
		if(lface == 0)
		{
			tblend = 0.5*(1.0 - t);
			dtblenddt = -0.5;
		}
		else
		{
			tblend = 0.5*(1.0 + t);
			dtblenddt = 0.5;
		}
		gface = elem->GetLocalFaceGlobalId(lface);
		pmax = FaceGetpmax(gface);
		n = pmax - 3;
		if (pmax < 3)
			continue;
		nodeI = elem->GetLocalFaceGlobalNodeOrder(lface, 0);
		nodeJ = elem->GetLocalFaceGlobalNodeOrder(lface, 1);
		nodeK = elem->GetLocalFaceGlobalNodeOrder(lface, 2);
		LI = lv[nodeI];
		dLIdr = dldr[nodeI];
		dLIds = dlds[nodeI];
		LJ = lv[nodeJ];
		dLJdr = dldr[nodeJ];
		dLJds = dlds[nodeJ];
		LK = lv[nodeK];
		dLKdr = dldr[nodeK];
		dLKds = dlds[nodeK];
		rf = LJ - LI;
		sf = 2.0*LK - 1.0;
		LegendrePns(pmax, rf, Pr);
		LegendrePns(pmax, sf, Ps);
		fun0 = fun;
		lll = LI*LJ*LK;
		tmp = lll*tblend;
		tmpt = dtblenddt*lll;
		if (calltype != derivonly)
		{
			for (i = 0; i <= n; i++)
			{
				for(j = 0; j <= n; j++)
				{
					if ((i + j) > n)
						break;
					pipj = Pr[i] * Ps[j];
					basisvec[fun] = tmp*pipj;
					fun++;
				}
			}
		}
		if(calltype != basisonly)
		{
			fun = fun0;
			LegendrePnDerivs(n, Pr, dpdr);
			LegendrePnDerivs(n, Ps, dpds);
			tmpr1 = (dLIdr*LJ*LK + dLJdr*LI*LK + dLKdr*LI*LJ)*tblend;
			tmps1 = (dLIds*LJ*LK + dLJds*LI*LK + dLKds*LI*LJ)*tblend;
			tmpr2 = (dLJdr - dLIdr)*tmp;
			tmps2 = (dLJds - dLIds)*tmp;
			tmpr3 = 2.0*dLKdr*tmp;
			tmps3 = 2.0*dLKds*tmp;
			for (i = 0; i <= n; i++)
			{
				double pi = Pr[i];
				double tmpr2dpi = tmpr2*dpdr[i];
				double tmps2dpi = tmps2*dpdr[i];
				double tmpr3pi = tmpr3*pi;
				double tmps3pi = tmps3*pi;
				for (j = 0; j <= n; j++)
				{
					if ((i + j) > n)
						break;
					pipj = pi * Ps[j];
					dbasisdr[fun] = tmpr1*pipj + tmpr2dpi * Ps[j] + tmpr3pi*dpds[j];
					dbasisds[fun] = tmps1*pipj + tmps2dpi * Ps[j] + tmps3pi*dpds[j];
					dbasisdt[fun] = tmpt*pipj;
					fun++;
				}
			}
		}

	}
	//quad faces:
	int k;
	for (lface = 2; lface < 5; lface++)
	{
		gface = elem->GetLocalFaceGlobalId(lface);
		pmax = FaceGetpmax(gface);
		n = pmax - 2;
		double sfl = t;
		//assign node numbers on tri face to use for determing rfl, see notes, p 5:
		if (lface == 2)
		{
			nodeI = 0;
			nodeJ = 1;
		}
		else if (lface == 3)
		{
			nodeI = 1;
			nodeJ = 2;
		}
		else if (lface == 4)
		{
			nodeI = 0;
			nodeJ = 2;
		}

		LI = lv[nodeI];
		dLIdr = dldr[nodeI];
		dLIds = dlds[nodeI];
		LJ = lv[nodeJ];
		dLJdr = dldr[nodeJ];
		dLJds = dlds[nodeJ];
		double rfl = LJ - LI;
		rej = rfl;
		if (((1.0 - rej)<RELSMALL) || ((rej + 1.0)<RELSMALL))
		{
			fun0 = fun;
			if(calltype != derivonly)
			{
				for(i = 2; i <= n; i++)
				{
					for(j = 2; j <= n; j++)
					{
						if ((i + j) <= pmax)
						{
							basisvec[fun] = 0.0;
							fun++;
						}
					}
				}
			}
			if(calltype != basisonly)
			{
				//asymptotic branch not supported for derivs
				fun = fun0;
				for (i = 2; i <= n; i++)
				{
					for (j = 2; j <= n; j++)
					{
						if ((i + j)>pmax)
							break;
						dbasisdr[fun] = 0.0;
						dbasisds[fun] = 0.0;
						dbasisdt[fun] = 0.0;
						fun++;
					}
				}

			}
			continue;
		}
		double drfldr = dLJdr - dLIdr;
		double drflds = dLJds - dLIds;
		elem->GetQuadFaceGlobalCoords(lface,rfl,sfl,rf,sf);
		if (calltype != derivonly)
		{
			Basisfuns1d(pmax, rf, phiI);
			Basisfuns1d(pmax, sf, phiJ);
		}
		if (calltype != basisonly)
		{
			Basisfuns1d(pmax, rf, phiI, phiprimeI, true);
			Basisfuns1d(pmax, sf, phiJ, phiprimeJ, true);
		}
		double rm1 = 1.0 / (1.0 - rej);
		double rp1 = 1.0 / (1.0 + rej);
		r2f = 4.0*rm1*rp1;
		double llr2f = LI*LJ*r2f;
		double dllr2fdr, dllr2fds;
		drejdr = drfldr;
		drejds = drflds;
		fun0 = fun;
		if(calltype != derivonly)
		{
			for (i = 2; i <= n; i++)
			{
                double pi = phiI[i];
				for(j = 2; j <= n; j++)
				{
					if ((i + j) <= pmax)
					{
						basisvec[fun] = llr2f*pi*phiJ[j];
						fun++;
					}
				}
			}
		}
		if(calltype != basisonly)
		{
			fun = fun0;
			dllr2fdr = r2f*(0.5*rej*llr2f*drejdr + (dLIdr*LJ + dLJdr*LI));
			dllr2fds = r2f*(0.5*rej*llr2f*drejds + (dLIds*LJ + dLJds*LI));
			double drfdrfl, drfdsfl, dsfdrfl, dsfdsfl;
			elem->GetQuadFaceDerivs(lface, rfl, sfl, drfdrfl, drfdsfl, dsfdrfl, dsfdsfl);
			double drfdr = drfdrfl*drfldr; //omitted term drfdsfl*dsfldr because dsfldr=0
			double drfds = drfdrfl*drflds; //omitted term drfdsfl*dsflds because dsflds=0
			double drfdt = drfdsfl; //omitted term drfdrfl*drfldt because drfldt=0, note dsfldt = 1
			double dsfdr = dsfdrfl*drfldr; //omitted term dsfdsfl*dsfldr because dsfldr=0
			double dsfds = dsfdrfl*drflds; //omitted term dsfdsfl*dsflds because dsflds=0
			double dsfdt = dsfdsfl; //omitted term dsfdrfl*drfldt because drfldt=0, note dsfldt = 1
			for (i = 2; i <= n; i++)
			{
                double pi = phiI[i];
                double dpi = phiprimeI[i];
				for(j = 2; j <= n; j++)
				{
					if ((i + j) > pmax)
						break;
					pipj = pi*phiJ[j];
					dbasisdr[fun] = dllr2fdr*pipj +	llr2f*(dpi*phiJ[j] * drfdr + pi*phiprimeJ[j] * dsfdr);
					dbasisds[fun] = dllr2fds*pipj +	llr2f*(dpi*phiJ[j] * drfds + pi*phiprimeJ[j] * dsfds);
					dbasisdt[fun] = llr2f*(dpi*phiJ[j] * drfdt + pi*phiprimeJ[j] * dsfdt); // (dll2rfdt = 0)
					fun++;
				}
			}
		}
	}


	//volume functions
	if(elpmax < 5)
		return fun;
	lll = L1*L2*L3;
	if (calltype == basisonly)
		Basisfuns1d(elpmax, t, phiI);
	else
		Basisfuns1d(elpmax, t, phiI, phiprimeI, true);
	n = elpmax - 5;
	int n3 = elpmax - 3;
	LegendrePns(n, r, Pr);
	LegendrePns(n, s, Ps);
	fun0 = fun;
	if(calltype != derivonly)
	{
		for (i = 0; i <= n; i++)
		{
            double pi = Pr[i];
			for(j = 0; j <= n; j++)
			{
				pipj = pi*Ps[j];
				for (k = 2; k <= n3; k++)
				{
					if ((i + j + k) > n3)
						break;
					basisvec[fun] = lll*pipj*phiI[k];
					fun++;
				}
			}
		}
	}
	if(calltype != basisonly)
	{
		fun = fun0;
		tmpr1 = DL1DR*L2*L3 + DL2DR*L1*L3 + DL3DR*L1*L2;
		tmps1 = DL1DS*L2*L3 + DL2DS*L1*L3 + DL3DS*L1*L2;
		LegendrePnDerivs(n, Pr, dpdr);
		LegendrePnDerivs(n, Ps, dpds);
		for (i = 0; i <= n; i++)
		{
            double pi = Pr[i];
            double dpi = dpdr[i];
			for (j = 0; j <= n; j++)
			{
				pipj = pi*Ps[j];
                double dpiPj = dpi*Ps[j];
                double dpjPi = pi*dpds[j];
				for (k = 2; k <= n3; k++)
				{
					if ((i + j + k)>n3)
						break;
					dbasisdr[fun] = ((tmpr1*pipj) + lll*dpiPj)*phiI[k];
					dbasisds[fun] = ((tmps1*pipj) + lll*dpjPi)*phiI[k];
					dbasisdt[fun] = lll*pipj*phiprimeI[k];
					fun++;
				}
			}
		}
	}

	return fun;
}

//number of wedge volume functions for p-orders 5 through 10.
//determined numerically in project "wejvolfuncs":
static int SRnumWedgeVolFuncs[6] = { 1, 4, 10, 20, 35, 56 };

int SRbasis::CountWedgeFunctions(SRelement *elem, int &nint)
{
	//count basis functions for a wedge element
	//input:
		//elem = pointer to element
	//output:
		//nint = number of internal basis functions in element
	//return:
		//total number of basis functions in element

	int nfun, lej, pej, lface, n;
	int elpmax = 0;
	SRedge* edge;
	//corner functions:
	nfun = 6;
	//edge functions:
	for(lej = 0; lej < 9; lej++)
	{
		edge = elem->GetEdge(lej);
		pej = edge->pOrder;
		if(pej > elpmax)
			elpmax = pej;
		nfun += (pej - 1);
	}

	//face functions:
	SRface* face;
	for(lface = 0; lface < 5; lface++)
	{
		face = elem->GetFace(lface);
		n = CountFaceInternalFunctions(face);
		nfun += n;
	}

	//volume functions
	if(elpmax < 5)
		nint=0;
	else
		nint = SRnumWedgeVolFuncs[elpmax - 5];
	nfun += nint;

	return nfun;
}
