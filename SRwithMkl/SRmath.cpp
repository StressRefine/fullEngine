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
// SRmath.cpp: implementation of the SRmath class.
//
//////////////////////////////////////////////////////////////////////

#include "SRmodel.h"
#ifndef NOSOLVER
#include "mkl.h"
#endif

extern SRmodel model;

int SRdoubleCompareFunc(const void* vd1, const void* vd2)
{
	double v1 = *((double*)vd1);
	double v2 = *((double*)vd2);
	if (v1 < v2)
		return -1;
	else if (v1 == v2)
		return 0;
	else
		return 1;
}

void SRvec3::Sort()
{
	//binary sort this vector
	qsort((void*)d, 3, sizeof(double), SRdoubleCompareFunc);
}

void SRvec3::localBasisFromFixedTang(SRvec3 tang, SRvec3 norm, SRvec3& e1, SRvec3& e2)
{
	tang.Cross(norm, e1);
	e1.Normalize();
	norm.Cross(e1, e2);
}

double SRvec3::Dot(double* x, double* y)
{
	//dot product between 3d vectors x and y
	//input:
		//x, y = vectors
	//return:
		//dot product

	double dot = 0.0;
	for (int i = 0; i < 3; i++)
		dot += x[i] * y[i];
	return dot;
}

double SRvec3::TwoNorm(double* v)
{
	//calculate the two-norm of a 3d vector
	//input:
		//v = vector
	//return:
		//two-norm of v

	double vn = 0.0;
	for(int i = 0; i < 3; i++)
		vn += v[i] * v[i];
	vn = sqrt(vn);
	return vn;
}

void SRvec3::Copy(double v[])
{
	//copy this into a vector of doubles
	for (int dof = 0; dof < 3; dof++)
		v[dof] = d[dof];
}

void SRvec3::Zero(double* v)
{
	//zero a vec3

	for (int i = 0; i < 3; i++)
		v[i] = 0.0;
}

void SRvec3::Cross(double* a, double* b, double* c)
{
	//cross 2 vectors
	//input:
		//a,b = vectors
	//output:
		//c = a cross b

	c[0] = a[1] * b[2] - b[1] * a[2];
	c[1] = a[2] * b[0] - b[2] * a[0];
	c[2] = a[0] * b[1] - b[0] * a[1];
}

double SRvec3::Normalize(double* v)
{
	//normalize a 3d vector
	//input:
		//v = vector
	//return:
		//norm of v

	double norm = TwoNorm(v);
	int i;
	double d = 1.0/norm;
	for(i = 0; i < 3; i++)
		v[i] *= d;
	return norm;
}

void SRvec3::Rotate(SRmat33& R, bool transpose)
{
	//v = this vector, perform v = R*v, overriding contents of this vector
	//input:
		//R = 3x3 rotation matrix
		//transpose = true to use R-transpose instead of R, else false

	SRvec3 v2;
	v2.Copy(*(this));
	Rotate(R,v2,transpose);
	Copy(v2);
}

void SRvec3::Rotate(SRmat33& R, SRvec3& v2, bool transpose)
{
	//v = this vector, perform v2 = R*v
	//input:
		//R = 3x3 rotation matrix
		//transpose = true to use R-transpose instead of R, else false
	//output:
		//v2 = rotated vector stored as SRvec3

	Rotate(R, v2.d, transpose);
}

void SRvec3::Rotate(SRmat33& R, double v2[], bool transpose)
{
	//v = this vector, perform v2 = R*v
	//input:
		//R = 3x3 rotation matrix
		//transpose = true to use R-transpose instead of R, else false
	//output:
		//v2 = rotated vector 

	int i, j;
	for (i = 0; i < 3; i++)
	{
		v2[i] = 0.0;
		if (transpose)
		{
			for (j = 0; j < 3; j++)
				v2[i] += (R.rows[j].d[i] * d[j]);
		}
		else
		{
			for (j = 0; j < 3; j++)
				v2[i] += (R.rows[i].d[j] * d[j]);
		}
	}
}

void SRvec3::Rotate(SRvec3 e1, SRvec3 e2, SRvec3 e3, SRvec3& v2)
{
	//rotate this vector to v2 using basis vectors e1,e2,e3
	v2.d[0] += d[0] * e1.d[0] + d[1] * e2.d[0] + d[2] * e3.d[0];
	v2.d[1] += d[0] * e2.d[1] + d[1] * e2.d[1] + d[2] * e3.d[1];
	v2.d[2] += d[0] * e3.d[2] + d[1] * e2.d[2] + d[2] * e3.d[2];
}

void SRvec3::Subtract(SRvec3& v2, SRvec3& v3)
{
	//v = this, perform v3 = v - v2
	//input:
		//v2 = vector to subtract
	//output:
		//v3 = vector resulting from subtraction

	v3.d[0] = d[0] - v2.d[0];
	v3.d[1] = d[1] - v2.d[1];
	v3.d[2] = d[2] - v2.d[2];
}

void SRvec3::Add(SRvec3& v2, SRvec3 &v3)
{
	//v = this, perform v3 = v + v2
	//input:
		//v2 = vector to add
	//output:
		//v3 = vector resulting from addition

	v3.d[0] = d[0] + v2.d[0];
	v3.d[1] = d[1] + v2.d[1];
	v3.d[2] = d[2] + v2.d[2];
}

double SRvec3::Distance(SRvec3& v2)
{
	//calculate distance from this vector to v2

	double dx, dy, dz;
	dx = v2.d[0] - d[0];
	dy = v2.d[1] - d[1];
	dz = v2.d[2] - d[2];
	return sqrt(dx*dx + dy*dy + dz*dz);
}

double SRvec3::Distance(double v1, double v2, double v3)
{
	//calculate distance from this vector to v2

	double dx, dy, dz;
	dx = v1 - d[0];
	dy = v2 - d[1];
	dz = v3 - d[2];
	return sqrt(dx*dx + dy*dy + dz*dz);
}

bool SRvec3::Equals(SRvec3& v2)
{
	//see if this vector equals v2

	SRvec3 dv;
	//dv = "this" - v2:
	Subtract(v2, dv);
	if (dv.Length() < TINY)
		return true;
	else
		return false;
}

void SRmat33::Assign(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33)
{
	//assign a mat33 using components m11...
	rows[0].d[0] = m11;
	rows[0].d[1] = m12;
	rows[0].d[2] = m13;
	rows[1].d[0] = m21;
	rows[1].d[1] = m22;
	rows[1].d[2] = m23;
	rows[2].d[0] = m31;
	rows[2].d[1] = m32;
	rows[2].d[2] = m33;
}

void SRmat33::setIdentity()
{
	Zero();
	rows[0].d[0] = 1.0;
	rows[1].d[1] = 1.0;
	rows[2].d[2] = 1.0;
}


void SRmat33::Mult(SRmat33& b, SRmat33& c, bool transpose)
{
	// a ="this" matrix
	//if(transpose) then c = a*b-transpose else c = a times b

	int i, j, k;
	double ct;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			ct = 0.0;
			if (transpose)
			{
				for (k = 0; k < 3; k++)
					ct += rows[i].d[k] * b.rows[j].d[k];
			}
			else
			{
				for (k = 0; k < 3; k++)
					ct += rows[i].d[k] * b.rows[k].d[j];
			}
			c.rows[i].d[j] = ct;
		}
	}
}

void SRmat33::Mult(SRvec3& vin, SRvec3& vout)
{
	//matrix-vector multiply
	// a ="this" matrix
	//input:
		//vin = vector to multiply by a
	//output:
		//vout = a*vin

	vout.Zero();
	int i;
	for (i = 0; i < 3; i++)
		vout.d[i] = vin.Dot(rows[i].d);
}

void SRmat33::Subtract(SRmat33& m2, SRmat33& m3)
{
	//m = this, perform m3 = m - m2
	//input:
		//m2 = matrix to subtract
	//output:
		//m3 = matrix resulting from subtraction

	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
			m3.rows[i].d[j] = rows[i].d[j] - m2.rows[i].d[j];
	}
}

void SRmat33::Scale(double s)
{
	//scale this matrix by scalar s

	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
			rows[i].d[j] *= s;
	}
}

void SRmat33::Copy(SRvec3& e1, SRvec3& e2, SRvec3& e3, bool columns)
{
	//copy vectors e1,e2,e3 into this matrix.
	//input:
		//e1,e2,e3 = 3-vectors
		//columns = true if e1,e2,e3 should be columns of matrix else false

	int i;
	if(columns)
	{
		for (i = 0; i < 3; i++)
		{
			rows[i].d[0] = e1.d[i];
			rows[i].d[1] = e2.d[i];
			rows[i].d[2] = e3.d[i];
		}
	}
	else
	{
		rows[0].Assign(e1.d);
		rows[1].Assign(e2.d);
		rows[2].Assign(e3.d);
	}
}

void SRmat33::vecCopy(int rc, SRvec3&v, bool copyRow)
{
	//copy a row or column of this marix into a vector
	//input:
		//rc = row or column number
		//copyRow = true to copy row, false for column
	if (copyRow)
	{
		for (int i = 0; i < 3; i++)
			v.d[i] = rows[rc].d[i];
	}
	else
	{
		for (int i = 0; i < 3; i++)
			v.d[i] = rows[i].d[rc];
	}
}


SRmath::SRmath()
{
	identityMatrix.setIdentity();
}

void SRmath::Setup()
{
	int np = model.GetMaxNumThreads();
	gpData.Allocate(np);
	for (int i = 0; i < np; i++)
	{
		SRGaussData *gpd = gpData.GetPointer(i);
		gpd->gp3d.Allocate(MAXGAUSSPOINTS);
		gpd->gp2d.Allocate(MAXGAUSSPOINTS2D);
	}
}


double SRmath::Distance(SRvec3& p1, SRvec3& p2)
{
	//distance between two points

	return p1.Distance(p2);
}


double SRmat33::Determinant()
{
	//Determinant of this 3x3 matrix

	return rows[0].d[0] * (rows[1].d[1] * rows[2].d[2] - rows[2].d[1] * rows[1].d[2])
		- rows[1].d[0] * (rows[0].d[1] * rows[2].d[2] - rows[2].d[1] * rows[0].d[2])
		+ rows[2].d[0] * (rows[0].d[1] * rows[1].d[2] - rows[1].d[1] * rows[0].d[2]);
}

bool SRmat33::Solve(double b[], double x[])
{
	//solve Ax = b for x where A is this mat33, using Cramer's rule

	double d, dx;
	double colsave[3];
	int i, j;
	d = Determinant();
	if (fabs(d) < TINY)
		return false;
	d = 1.0 / d;
	for (j = 0; j < 3; j++)
	{
		for (i = 0; i < 3; i++)
		{
			colsave[i] = rows[i].d[j];
			rows[i].d[j] = b[i];
		}
		dx = Determinant();
		x[j] = dx*d;
		for (i = 0; i < 3; i++)
			rows[i].d[j] = colsave[i];
	}
	return true;
	//tuning: could save some ops (< 1/3) by sharing co-factor mults between determinant calculation
	//and solution for x[0]. but code is then much less readable (see old stress refine source)
	//only do it this routine scores high in profile
}

int SRmath::Round(double d)
{
	//round off a double to an int

	double r;
	int i= (int) d;
	r = d-i;
	if(r < 0.5)
		return i;
	else
		return i + 1;
}



void SRmat33::Transpose(SRmat33& mout)
{
	//transpose this matrix
	//output:
		//mout = transpose
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
			mout.rows[i].d[j] = rows[j].d[i];
	}
}

double SRmat33::Invert()
{
	//invert this matrix
	//note:
		//this matrix is overwritten with its inverse

	double det, inv11, inv12, inv13, inv21, inv22, inv23, inv31, inv32, inv33;
	inv11 = rows[1].d[1] * rows[2].d[2] - rows[2].d[1] * rows[1].d[2];//co-factor of 11
	inv21 = rows[2].d[0] * rows[1].d[2] - rows[1].d[0] * rows[2].d[2];//co-factor of 12
	inv31 = rows[1].d[0] * rows[2].d[1] - rows[2].d[0] * rows[1].d[1];//co-factor of 13
	det = rows[0].d[0] * inv11 + rows[0].d[1] * inv21 + rows[0].d[2] * inv31;
	if(fabs(det) < TINY)
		ERROREXIT;
	double detm1 = 1.0 / det;
	inv12 = rows[2].d[1] * rows[0].d[2] - rows[0].d[1] * rows[2].d[2];//co-factor of 21
	inv13 = rows[0].d[1] * rows[1].d[2] - rows[1].d[1] * rows[0].d[2];//co-factor of 31
	inv22 = rows[0].d[0] * rows[2].d[2] - rows[2].d[0] * rows[0].d[2];//co-factor of 22
	inv23 = rows[1].d[0] * rows[0].d[2] - rows[0].d[0] * rows[1].d[2];//co-factor of 32
	inv32 = rows[2].d[0] * rows[0].d[1] - rows[0].d[0] * rows[2].d[1];//co-factor of 23
	inv33 = rows[0].d[0] * rows[1].d[1] - rows[1].d[0] * rows[0].d[1];//co-factor of 33
	rows[0].d[0] = inv11*detm1;
	rows[1].d[0] = inv21*detm1;
	rows[2].d[0] = inv31*detm1;
	rows[0].d[1] = inv12*detm1;
	rows[1].d[1] = inv22*detm1;
	rows[2].d[1] = inv32*detm1;
	rows[0].d[2] = inv13*detm1;
	rows[1].d[2] = inv23*detm1;
	rows[2].d[2] = inv33*detm1;
	return det;
}

int SRmath::FillGaussPoints(SRelement* elem)
{
	//Calculate and store gauss points for an element
	//input:
		//elem=pointer to element
	//return:
		//number of integration points

	if (elem->GetType() == tet)
		return FillTetGaussPoints(elem);
	else if (elem->GetType() == wedge)
		return FillWedgeGaussPoints(elem);
	else
		return FillBrickGaussPoints(elem);
}

int SRmath::FillBrickGaussPoints(SRelement* elem)
{
	//Calculate and store gauss points for a brick element
	//input:
		//elem = pointer to element
	//return:
		//number of integration points	

	int numgauss1d = elem->GetMaxPorder()+1;
	SRgaussPoint* g1d = gp1d.GetPoints(numgauss1d);
	int i, j, k, nint = 0;
	double wr, ws, wt, r, s, t, w;
	for (i = 0; i < numgauss1d; i++)
	{
		r = g1d[i].x;
		wr = g1d[i].w;
		for(j = 0; j < numgauss1d; j++)
		{
			s = g1d[j].x;
			ws = g1d[j].w;
			for(k = 0; k < numgauss1d; k++)
			{
				t = g1d[k].x;
				wt = g1d[k].w;
				w = wr*ws*wt;
				PutGP3d(nint, r, s, t, w, elem->GetThread());
				nint++;
			}
		}
	}
	return nint;
}

int SRmath::FillWedgeGaussPoints(SRelement* elem)
{
	//Calculate and store gauss points for a wedge element
	//for using degenerated gauss quadrature
	//input:
		//elem = pointer to element
	//return:
		//number of integration points

	int pmax = elem->GetMaxPorder();
	int numgauss1d = pmax+1;
	//fill up 2d Gauss points for the triangular cross section:
	int nint2d = FillFaceGaussPoints(pmax,3, elem->GetThread());
	int i, j, nint = 0;
	//blend 2d gauss points with one D in "t" direction:
	double r, s, t, wrs, wt, w;
	SRgaussPoint* g1d = gp1d.GetPoints(numgauss1d);
	for (i = 0; i < nint2d; i++)
	{
		GetGP2d(i, r, s, wrs, elem->GetThread());
		for (j = 0; j < numgauss1d; j++)
		{
			t = g1d[j].x;
			wt = g1d[j].w;
			w = wrs*wt;
			PutGP3d(nint, r, s, t, w, elem->GetThread());
			nint++;
		}
	}
	return nint;
}

int SRmath::FillTetGaussPoints(SRelement* elem)
{
	//Calculate and store gauss points for a tetrahedral element (by degenerating
	//brick gauss points)
	//input:
		//elem=pointer to element
	//return:
		//number of integration points

	int numgauss1d = elem->GetMaxPorder()+1;
	SRgaussPoint* g1d = gp1d.GetPoints(numgauss1d);
	int i, j, k, nint = 0;
	double wr, ws, wt, r, s, t, w, rb, sb, tb;
	for(i = 0; i < numgauss1d; i++)
	{
		rb = g1d[i].x;
		wr = g1d[i].w;
		for(j = 0; j < numgauss1d; j++)
		{
			sb = g1d[j].x;
			ws = g1d[j].w;
			for(k = 0; k < numgauss1d; k++)
			{
				tb = g1d[k].x;
				wt = g1d[k].w;
				w = wr*ws*wt;
				//map rb,sb,tb in degenerated brick to corresponding position in
				//tet:
				w *= model.map.BrickTetMap(rb, sb, tb, r, s, t);
				PutGP3d(nint, r, s, t, w, elem->GetThread());
				nint++;
			}
		}
	}
	return nint;
}

int SRmath::FillGaussPoints(SRface* face, int pIn, int thread)
{
	//Calculate and store gauss points for a face
	//input:
		//face = pointer to face
	//return:
		//number of integration points

	int maxp = model.basis.FaceGetpmax(face);
	int p = maxp;
	if (pIn != -1)
		p = pIn;
	return FillFaceGaussPoints(p, face->GetNumLocalEdges(), thread);
}

int SRmath::FillFaceGaussPoints(int maxPorder, int numSides, int thread)
{
	//Calculate and store gauss points for a face
	//input:
		//face = maxPorder = max p on face
		//numSides = 3 or 4
	//return:
		//number of integration points

	int numgauss1d = maxPorder + 1;
	SRgaussPoint* g1d = gp1d.GetPoints(numgauss1d);

	int i, j, nint = 0;
	double wr, ws, r, s, w, rq, sq;
	if (numSides == 3)
	{
		//triangular face
		for(i = 0; i < numgauss1d; i++)
		{
			rq = g1d[i].x;
			wr = g1d[i].w;
			for(j = 0; j < numgauss1d; j++)
			{
				sq = g1d[j].x;
				ws = g1d[j].w;
				w = wr*ws;
				//map rq,sq in degenerated quad to corresponding position in
				//triangle:
				w *= model.map.QuadTriMap(rq, sq, r, s);
				PutGP2d(nint, r, s, w, thread);
				nint++;
			}
		}
	}
	else
	{
		for(i = 0; i < numgauss1d; i++)
		{
			r = g1d[i].x;
			wr = g1d[i].w;
			for(j = 0; j < numgauss1d; j++)
			{
				s = g1d[j].x;
				ws = g1d[j].w;
				w = wr*ws;
				PutGP2d(nint, r, s, w, thread);
				nint++;
			}
		}
	}
	return nint;
}

int SRmath::CountGaussPoints(SRface* face)
{
	//Calculate number of gauss points for a face
	//input:
	 	//face = pointer to face
	//return:
		//number of integration points

	int maxp = model.basis.FaceGetpmax(face);
	int numgauss1d = maxp + 1;
		return numgauss1d*numgauss1d;
}

int SRmath::CountGaussPoints(SRelement* elem)
{
	//Calculate number of gauss points for an element
	//input:
		//elem = pointer to element
	//return:
		//number of integration points

	int maxp = 0;
	for (int i = 0; i < elem->GetNumLocalEdges(); i++)
	{
		int p = elem->GetLocalEdgePOrder(i);
		if (p > maxp)
			maxp = p;
	}
	int ngp1d = maxp + 1;
	return ngp1d*ngp1d*ngp1d;
}

void SRmath::PutGP2d(int i, double r, double s, double w, int thread)
{
	//store a 2d gauss point in the class variable gp2d
	//input:
		//r,s = coordinates
		//w = weight

	SRgaussPoint2D* pgp = gpData.GetPointer(thread)->gp2d.GetPointer(i);
	pgp->Assign(r, s, w);
}

void SRmath::GetGP2d(int i, double& r, double& s, double& w, int thread)
{
	//look up a 2d gauss point in the class variable gp2d
	//output:
		//r,s = coordinates
		//w = weight

	SRgaussPoint2D* pgp = gpData.GetPointer(thread)->gp2d.GetPointer(i);
	pgp->GetRsw(r, s, w);
}

void SRmath::PutGP3d(int i, double r, double s, double t, double w, int thread)
{
	//store a 3d gauss point in the class variable gp3d
	//input:
		//r,s,t = coordinates
		//w = weight
	SRGaussData* gpd = gpData.GetPointer(thread);
	SRgaussPoint3D* pgp = gpd->gp3d.GetPointer(i);
	pgp->Assign(r, s, t, w);
}

void SRmath::GetGP3d(int i, double& r, double& s, double& t, double& w, int thread)
{
	//look up a 3d gauss point in the class variable gp3d
	//output:
		//r,s,t = coordinates
		//w = weight

	SRGaussData* gpd = gpData.GetPointer(thread);
	SRgaussPoint3D* pgp = gpd->gp3d.GetPointer(i);
	pgp->GetRstw(r, s, t, w);
}

void SRmath::VectorCopy(int n, double a[], double b[])
{
	//copy double vector a to vector b
	//input:
		//a  = vector
		//n = length
	//output:
		//b = copy of a

	for(int i = 0; i < n; i++)
		b[i] = a[i];
}

double SRmath::VectorDot(int n, double a[], double b[])
{
	//dot product of double vector a and vector b
	//input:
		//a,b  = vectors
		//n = length
	//return:
		//a dot b

	double d = 0.0;
	for(int i = 0; i < n; i++)
		d += (b[i] * a[i]);
	return d;
}

double SRmath::Vector2Norm(int n, double a[])
{
	//two-norm of double vector a
	//input:
		//a  = vector
		//n = length
	//return:
		//two-norm of a

	double d = 0.0;
	for(int i = 0; i < n; i++)
		d += (a[i] * a[i]);
	d = sqrt(d);
	return d;
}

void SRmath::GetTraction(double norm[], double stress[], double tract[])
{
	//determine traction vector from stress and normal to a surface
	//input:
		//norm = surface normal (3d vector)
		//stress = stress tensor stored as vector
	//output:
		//tract = traction vector

	double sij[3][3];
	sij[0][0] = stress[0];
	sij[1][1] = stress[1];
	sij[2][2] = stress[2];
	sij[0][1] = sij[1][0] = stress[3];
	sij[0][2] = sij[2][0] = stress[4];
	sij[1][2] = sij[2][1] = stress[5];
	int i, j;
	for (i = 0; i < 3; i++)
	{
		tract[i] = 0.0;
		for (j = 0; j < 3; j++)
			tract[i] += sij[i][j] * norm[j];
	}
}


double SRmath::GetSvm(double stress[])
{
	//return svm of a stress tensor
	//input:
		//stress = stress tensor stored as vector
	//return:
		//svm = von Mise Stress


	double p, s11, s12, s13, s22, s23, s33, svm;
	p = ONETHIRD*(stress[0] + stress[1] + stress[2]);
	s11 = stress[0] - p;
	s22 = stress[1] - p;
	s33 = stress[2] - p;
	s12 = stress[3];
	s13 = stress[4];
	s23 = stress[5];
	svm = s11*s11 + s22*s22 + s33*s33 + 2.0*(s12*s12 + s13*s13 + s23*s23);
	svm = sqrt(1.5*svm);
	return svm;
}

void SRmath::GetPrinStress(double stress[], double& sp1, double& sp2)
{
	double sp3;
	//calculate sp1 = max prin, sp2 = min prin.
	GetPrinStress(stress, sp1, sp2, sp3);
	sp2 = sp3;
}

void SRmath::GetPrinStress(double stress[], double& sp1, double& sp2, double& sp3)
{
	//calculate sp1 = max prin, sp3 = min prin, sp2 = intermediate.
	//convention stressComp[6] = { "xx", "yy", "zz", "xy", "xz", "yz" };
	double I1, I2, I3, phi;
	double s1 = stress[0];
	double s2 = stress[1];
	double s3 = stress[2];
	double s12 = stress[3];
	double s13 = stress[4];
	double s23 = stress[5];
	double s1s2 = s1*s2;
	double s23sq = s23*s23;
	double s12sq = s12*s12;
	I1 = s1 + s2 + s3;
	I2 = s1s2 + s2*s3 + s3*s1 - s12sq - s23sq - s13*s13;
	I3  = s1s2*s3 - s1*s23sq - s2*s13*s13 - s3*s12sq + 2.0*s12*s23*s13;
	double den = I1*I1 - 3 * I2;
	den = den*den*den;
	den = 2.0*sqrt(den);
	phi = 0.3333333333333333333333*acos((2.0*I1*I1*I1 - 9.0*I1*I2 + 27.0*I3) / den);
	double I1Over3 = 0.3333333333333333333333*I1;
	double I1I2Factor = 0.6666666666666666666*sqrt(I1*I1 - 3.0*I2);
	s1 = I1Over3 + I1I2Factor*cos(phi);
	s2 = I1Over3 + I1I2Factor*cos(phi - 0.6666666666666666666*PI);
	s3 = I1Over3 + I1I2Factor*cos(phi - 1.3333333333333333333333*PI);

	//sort:
	SRvec3 svec;
	svec.Assign(s1, s2, s3);
	svec.Sort();
	sp1 = svec.d[2];
	sp2 = svec.d[1];
	sp3 = svec.d[0];
}


