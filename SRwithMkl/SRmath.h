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
// SRmath.h: interface for the SRmath class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRMATH_INCLUDED)
#define SRMATH_INCLUDED

#include <math.h>

#define MATHMIN(x,y) ((x<y)?x:y)
#define MATHMAX(x,y) ((x>y)?x:y)
    
//this is model.maxPorder+1:
#define MAX1DFUNCTIONS 9
//The following are all for maxPorder 8. if it is increased have to update
//these:
//max number of basis functions for an element.see Szabo, p 241
#define MAXELEMBASISFUNCTIONS 192
//This is for a quad face, See Szabo p 98.
#define MAXFACEBASISFUNCTIONS 51
#define MAXFACEEQUATIONS 153
#define MAXGAUSSPOINTS 729
#define MAXGAUSSPOINTS2D 81

#define BIG 1.e20
#define RELSMALL 1.e-7
#define SMALL 1.e-12
#define TINY 1.e-20
#define INVERSEMAP_ABSTOL 1.E-8
#define INVERSEMAP_RELTOL 1.E-6
#define INVERSEMAP_MAX_ITS 1000
#define SQRT3 1.732050807568877
#define SQRT3OVER2 0.866025403784438
#define SQRT3OVER4 0.433012701892219
#define SQRT3OVER3 0.5773502691896257
#define SQRT3OVER6 0.2886751345948128
#define SQRT2 1.414213562373095
#define SQRT2OVER2 0.707106781186547
#define SQRT6OVER6 0.408248290463862
#define SQRT6OVER12 0.204124145231931
#define SQRT3OVERSQRT8  0.6123724356957945
#define SQRT8OVER8 0.353553390593273
#define PI 3.141592653589793
#define TWOPI 6.283185307
#define PIOVER2 1.570796326794897
#define ONETHIRD 0.3333333333333333
#define TWOSQRTTWOTHIRDS 1.632993161855452
#define LOWSTRESSTOL 0.6666667
#define LOWSTRESSTOLFINAL 0.75

#define TRI_CENTROID_R 0.0
#define TRI_CENTROID_S SQRT3OVER3
#define TET_CENTROID_R 0.0
#define TET_CENTROID_S SQRT3OVER3
#define TET_CENTROID_T 0.5443310539518174 // = 2/3*sqrt(2/3)
#define TRI_FACE_S_FACTOR SQRT3
//derivatives of triangle area coordinates
#define DL1DR -0.5
#define DL1DS -SQRT3OVER6
#define DL2DR 0.5
#define DL2DS -SQRT3OVER6
#define DL3DR 0.0
#define DL3DS SQRT3OVER3
//derivatives of tetrahedral volume coordinates
#define DLT1DR -0.5
#define DLT1DS -SQRT3OVER6
#define DLT1DT -SQRT6OVER12
#define DLT2DR 0.5
#define DLT2DS -SQRT3OVER6
#define DLT2DT -SQRT6OVER12
#define DLT3DR 0.0
#define DLT3DS SQRT3OVER3
#define DLT3DT -SQRT6OVER12
#define DLT4DR 0.0
#define DLT4DS 0.0
#define DLT4DT SQRT3OVERSQRT8

class SRvec6
{
public:
	SRvec6(){ Zero(); };
	void Zero()
	{
		for (int i = 0; i < 6; i++)
			d[i] = 0.0;
	}
	double d[6];
};

class SRmat33;

class SRvec3
{
public:
	bool Equals(SRvec3& v2);
	void Assign(double v1, double v2, double v3){ d[0] = v1; d[1] = v2; d[2] = v3; };
	void Assign(double v[]){ d[0] = v[0]; d[1] = v[1]; d[2] = v[2]; };
	void PlusAssign(double v[]){ d[0] += v[0]; d[1] += v[1]; d[2] += v[2]; };
	void PlusAssign(SRvec3& v){ PlusAssign(v.d); };
	void PlusAssign(double v[], double scale){ d[0] += scale*v[0]; d[1] += scale*v[1]; d[2] += scale*v[2]; };
	void MinusAssign(double v[]){ d[0] -= v[0]; d[1] -= v[1]; d[2] -= v[2]; };
	void MinusAssign(SRvec3& v){ MinusAssign(v.d); };
	void Copy(SRvec3& v2){Assign(v2.d);};
	void Copy(double v[]);
	bool operator ==(SRvec3& v2){ return Equals(v2); };
	void operator =(SRvec3& v2){ Copy(v2); };
	void operator = (double v[]){ Assign(v); };
	void operator += (double v[]){ PlusAssign(v); };
	void operator += (SRvec3& v){ PlusAssign(v); };
	double Distance(SRvec3& v2);
	double Distance(double v1, double v2, double v3);
	void Add(SRvec3& v2, SRvec3& v3);
	void Subtract(SRvec3& v2,SRvec3& v3);
	void Cross(SRvec3& v2, SRvec3& v3){ Cross(d, v2.d, v3.d); };
	double Normalize(){ return Normalize(d); };
	void Zero(){ d[0] = 0.0; d[1] = 0.0; d[2] = 0.0; };
	double Length(){ return TwoNorm(d); };
	double Magnitude(){ return TwoNorm(d); };
	void Rotate(SRmat33& R, bool transpose = false);
	void Rotate(SRmat33& R, SRvec3& v2, bool transpose = false);
	void Rotate(SRmat33& R, double v2[], bool transpose = false);
	void Rotate(SRvec3 e1, SRvec3 e2, SRvec3 e3, SRvec3& v2);
	static double Normalize(double* v);
	static void Cross(double* a, double* b, double* c);
	static void Zero(double* v);
	static double TwoNorm(double* v);
	static double Dot(double* x, double* y);
	static void localBasisFromFixedTang(SRvec3 tang, SRvec3 norm, SRvec3& e1, SRvec3& e2);
	double Dot(SRvec3& v2){ return Dot( d, v2.d); };
	double Dot(double* x){ return Dot( d, x); };
	void Scale(double s)
	{
		for (int i = 0; i < 3; i++)
			d[i] *= s;
	};
	void Sort();
	SRvec3(){ Zero(); };
	SRvec3(double v[]){ Assign(v); };
	SRvec3(double v1, double v2, double v3){ Assign(v1, v2, v3); };
	~SRvec3(){};
	double d[3];
};

class SRmat33
{
public:
	double Invert();
	void Transpose(SRmat33& mout);
	void Assign(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33);
	bool Solve(double b[], double x[]);
	double Determinant();
	void Copy(SRvec3& e1, SRvec3& e2, SRvec3& e3, bool columns = true);
	void vecCopy(int rc, SRvec3&v, bool copyRow);
	void Scale(double s);
	void Subtract(SRmat33& m2, SRmat33& m3);
	void Mult(SRvec3& vin, SRvec3& vout);
	void Mult(SRmat33& b, SRmat33& c, bool transpose = false);
	void setIdentity();
	void Zero()
	{
		for (int i = 0; i < 3; i++)
			rows[i].Zero();
	};
	void operator +=(SRmat33& m2)
	{
		for (int i = 0; i < 3; i++)
			rows[i] += m2.rows[i];
	};
	void Copy(SRmat33& m2)
	{
		for (int i = 0; i < 3; i++)
			rows[i] = m2.rows[i];
	};
	void operator =(SRmat33& m2){ Copy(m2); };

	SRmat33(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33)
	{
		Assign(m11, m12, m13, m21, m22, m23, m31, m32, m33);
	};
	SRmat33(){};
	~SRmat33(){};

	SRvec3 rows[3];
};

class SRgaussPoint2D
{
public:
	void Assign(double rt, double st, double wt){ r = rt; s = st; w = wt; };
	void GetRsw(double& rt, double& st, double& wt){ rt = r; st = s; wt = w; };
private:
	double r;
	double s;
	double w;
};

class SRgaussPoint3D
{
public:
	void Assign(double rt, double st, double tt, double wt){ r = rt; s = st; t = tt; w = wt; };
	void GetRstw(double& rt, double& st, double& tt, double& wt){ rt = r; st = s; tt = t; wt = w; };
private:
	double r;
	double s;
	double t;
	double w;
};

struct SRgaussPoint
{
	double x;
	double w;
};

static SRgaussPoint gauss1d2[2]=
{
	{-0.577350269189626,1.0               },
	{0.577350269189626 ,1.0               },
};
static SRgaussPoint gauss1d3[3]=
{
	{-0.774596669241483,0.555555555555553 },
	{0.0               ,0.888888888888889 },
	{0.774596669241483 ,0.555555555555553 },
};
static SRgaussPoint gauss1d4[4]=
{
	{-0.861136311594053,0.347854845137448 },
	{-0.339981043584856,0.652145154862546 },
	{0.339981043584856 ,0.652145154862546 },
	{0.861136311594053 ,0.347854845137448 },
};
static SRgaussPoint gauss1d5[5]=
{
	{-0.906179845938664,0.236926885056182 },
	{-0.538469310105683,0.478628670499366 },
	{0.0               ,0.568888888888889 },
	{0.538469310105683 ,0.478628670499366 },
	{0.906179845938664 ,0.236926885056182 },
};
static SRgaussPoint gauss1d6[6]=
{
	{-0.932469514203152,0.171324492379162 },
	{-0.661209386466265,0.360761573048139 },
	{-0.238619186083197,0.467913934572689 },
	{0.238619186083197 ,0.467913934572689 },
	{0.661209386466265 ,0.360761573048139 },
	{0.932469514203152 ,0.171324492379162 },
};
static SRgaussPoint gauss1d7[7]=
{
	{-0.949107912342758,0.129484966168862 },
	{-0.741531185599394,0.279705391489277 },
	{-0.405845151377397,0.381830050505119 },
	{0.0               ,0.417959183673469 },
	{0.405845151377397 ,0.381830050505119 },
	{0.741531185599394 ,0.279705391489277 },
	{0.949107912342758 ,0.129484966168862 },
};
static SRgaussPoint gauss1d8[8]=
{
	{-0.960289856497536,0.10122853629037  },
	{-0.796666477413627,0.222381034453374 },
	{-0.525532409916329,0.313706645877887 },
	{-0.18343464249565 ,0.362683783378362 },
	{0.18343464249565  ,0.362683783378362 },
	{0.525532409916329 ,0.313706645877887 },
	{0.796666477413627 ,0.222381034453374 },
	{0.960289856497536 ,0.10122853629037  },
};
static SRgaussPoint gauss1d9[9]=
{
	{-0.968160239507626,0.0812743883615687},
	{-0.836031107326636,0.180648160694858 },
	{-0.61337143270059 ,0.260610696402936 },
	{-0.324253423403809,0.312347077040002 },
	{0.0               ,0.33023935500126  },
	{0.324253423403809 ,0.312347077040002 },
	{0.61337143270059  ,0.260610696402936 },
	{0.836031107326636 ,0.180648160694858 },
	{0.968160239507626 ,0.0812743883615687},
};
static SRgaussPoint gauss1d10[10]=
{
	{-0.973906528517172,0.0666713443086829},
	{-0.865063366688985,0.149451349150581 },
	{-0.679409568299024,0.219086362515982 },
	{-0.433395394129247,0.269266719309992 },
	{-0.148874338981631,0.295524224714753 },
	{0.148874338981631 ,0.295524224714753 },
	{0.433395394129247 ,0.269266719309992 },
	{0.679409568299024 ,0.219086362515982 },
	{0.865063366688985 ,0.149451349150581 },
	{0.973906528517172 ,0.0666713443086829},
};
static SRgaussPoint gauss1d11[11]=
{
	{-0.978228658146057,0.0556685671161696},
	{-0.887062599768095,0.125580369464905 },
	{-0.730152005574049,0.186290210927734 },
	{-0.519096129206812,0.23319376459199  },
	{-0.269543155952345,0.262804544510247 },
	{0.0               ,0.272925086777901 },
	{0.269543155952345 ,0.262804544510247 },
	{0.519096129206812 ,0.23319376459199  },
	{0.730152005574049 ,0.186290210927734 },
	{0.887062599768095 ,0.125580369464905 },
	{0.978228658146057 ,0.0556685671161696},
};

class SRgauss1d
{
public:
	//data:
	SRgaussPoint* gauss1d[10];
	SRgauss1d::SRgauss1d()
	{
		gauss1d[0] = gauss1d2; gauss1d[1] = gauss1d3; gauss1d[2] = gauss1d4;
		gauss1d[3] = gauss1d5; gauss1d[4] = gauss1d6; gauss1d[5] = gauss1d7;
		gauss1d[6] = gauss1d8; gauss1d[7] = gauss1d9; gauss1d[8] = gauss1d10;
		gauss1d[9] = gauss1d11;
	};
	SRgaussPoint* GetPoints(int n){ return gauss1d[n - 2]; };
};

class SRface;
class SRelement;

class SRGaussData
{
public:
	SRvector <SRgaussPoint3D> gp3d;
	SRvector <SRgaussPoint2D> gp2d;
};

class SRmath
{
public:
	SRmath();
	double Vector2Norm(int n, double a[]);
	double VectorDot(int n, double a[], double b[]);
	void VectorCopy(int n, double a[], double b[]);
	void GetGP3d(int i, double& r, double& s, double& t, double& w, int thread = 0);
	void PutGP3d(int i, double r, double s, double t, double w, int thread = 0);
	void GetGP2d(int i, double& r, double& s, double& w, int thread = 0);
	void PutGP2d(int i, double r, double s, double w, int thread = 0);
	int CountGaussPoints(SRface* face);
	int CountGaussPoints(SRelement* elem);
	int Round(double d);
	double Distance(SRvec3& p1, SRvec3& p2);
	int Sign(double x){ return (x>0.0) ? 1 : -1; };
	double Max(double x, double y){ return (x>y) ? x : y; };
	double Min(double x, double y){ return (x<y) ? x : y; };
	inline int Max(int x, int y){ return (x>y) ? x : y; };
	inline int Min(int x, int y){ return (x<y) ? x : y; };
	inline bool Even(int x){ return (((x / 2) * 2) == x) ? true : false; };
	inline bool Odd(int x){ return (((x / 2) * 2) != x) ? true : false; };
	inline bool Equal(double x, double y){ return (fabs(fabs(x) - fabs(y))<TINY) ? true : false; };
	int FillFaceGaussPoints(int maxPorder, int numSides, int thread = 0);
	int FillGaussPoints(SRface* face, int pIn = -1, int thread = 0);
	int FillGaussPoints(SRelement* elem);
	int FillTetGaussPoints(SRelement* elem);
	int FillWedgeGaussPoints(SRelement* elem);
	int FillBrickGaussPoints(SRelement* elem);
	SRgaussPoint* GetGaussPoints1d(int nint){ return gp1d.GetPoints(nint); };
	void GetTraction(double norm[], double stress[], double tract[]);
	double GetSvm(double stress[]);
	void GetPrinStress(double stress[], double& sp1, double& sp2);
	void GetPrinStress(double stress[], double& sp1, double& sp2, double& sp3);
	SRmat33& getIdentityMatrix(){ return identityMatrix; };
	void Setup();
private:
	SRmat33 identityMatrix;
	SRvector <SRGaussData> gpData;
	SRgauss1d gp1d;
};


#endif //!defined(SRMATH_INCLUDED)
