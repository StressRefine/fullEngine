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
// SRutil.h: interface for the SRutil class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRUTIL_INCLUDED)
#define SRUTIL_INCLUDED

#include <memory.h>
#include <vector>

using namespace std;

//temporary definitions to allow all new's and deletes to be
//easily found if later we decide to use a memory pool
#define ALLOCATEMEMORY new
#define DELETEMEMORY delete

class SRutil  
{
public:
	static void ErrorExit(const char* file, int line);
	SRutil(){};
	~SRutil(){};
};

class SRintVector
{
public:
	void Sort();
	int Find(int intIn);
	bool isEmpty(){ return (d == NULL); };
	int GetNum(){ return num; };
	void Free()
	{
		if(num == 0)
			return;
		DELETEMEMORY d;
		d = NULL;
		num = 0;
	};
	void Allocate(int nt)
	{
		if(nt == 0)
			return;
		if(d != NULL)
			Free();
		num = nt;
		d = ALLOCATEMEMORY int[num];
		Zero();
	};
	void PushBack(int v);
	void Zero()
	{
		for(int i = 0; i < num; i++)
			d[i] = 0;
	};
	void Set(int val)
	{
		for (int i = 0; i < num; i++)
			d[i] = val;
	};
	void Copy(SRintVector& v2, int numt = 0)
	{
		int len;
		if(numt == 0)
			len = v2.GetNum();
		else
			len = numt;
		if(num < len)
			Allocate(len);
		for(int i = 0; i < len; i++)
			d[i] = v2.d[i];
	};
	int* GetVector(){ return d; };
	inline int Get(int i){ return d[i];};
	inline void Put(int i,int di){ d[i] = di; };
	inline void PlusAssign(int i, int di){ d[i] += di; };
	inline void PlusPlus(int i){ d[i] ++; };
	int operator [] (int i) { return Get(i); };

	SRintVector(int nt){ Allocate(nt); };
	SRintVector(){ num = 0; d = NULL; };
	~SRintVector()
    {
		Free();
	};

	int* d;
	int num;
};

class SRdoubleVector
{
public:
	bool isEmpty(){ return (d == NULL); };
	void EquateVector(double* v){ d = v; };
	int GetNum(){ return num; };
	void Copy(SRdoubleVector& v2){ Copy(v2.d, v2.GetNum()); };
	void PushBack(double v);
	void Free()
	{
		try
		{
			if (num == 0)
				return;
			DELETEMEMORY d;
			d = NULL;
			num = 0;
		}
		catch (...)
		{
			return;
		}
	};
	void Allocate(int nt)
	{
		if(nt == 0)
			return;
		if(d != NULL)
			Free();
		num = nt;
		d = ALLOCATEMEMORY double[num];
		Zero();
	};
	void Zero()
	{
		for (int i = 0; i < num; i++)
			d[i] = 0;
	};
	void Set(double val)
	{
		for (int i = 0; i < num; i++)
			d[i] = val;
	};
	void Copy(double* v2, int len)
	{
		if(num < len)
			Allocate(len);
		for(int i = 0; i < len; i++)
			d[i] = v2[i];
	};
	double* GetVector(){ return d; };
	inline double Get(int i){ return d[i]; };
	inline void Put(int i, double di){ d[i] = di; };
	inline void PlusAssign(int i, double di){ d[i] += di; };
	double operator [] (int i) { return Get(i); };

	SRdoubleVector(int nt){ Allocate(nt); };
	SRdoubleVector(){ num = 0; d = NULL; };
	~SRdoubleVector(){ Free(); };

	double* d;
	int num;
};

//note: it is tricky to layer an SRvector over another vector class such as SRintVector
//e.g. for SRvector <SRintVector> tmp,
//don't work with tmp[i]
//safer to do SRintVector *oneTmp = tmp.GetPointer(i) then work with oneTmp.
template <class gen>
class SRvector
{
public:
	bool isEmpty(){ return (d.size() == 0); };
	void Free()
	{
		if (!isEmpty())
			d.clear();
	};
	void Allocate(int nt)
	{
		d.resize(nt);
	};
	int GetNum(){ return d.size(); };
	inline gen* GetVector(){ return &d[0]; };
	inline gen *GetPointer(int i){ return &d[i]; };
	gen& Get(int i){ return d[i]; };
	inline void Put(int i, gen &di){
		d[i] = di;
	};
	inline void PlusAssign(int i, gen &di){ d[i] += di; };
	void pushBack(gen dt){ d.push_back(dt); };
	gen operator [] (int i) { return Get(i); };
	SRvector(int nt){ d.resize(0); Allocate(nt); };
	SRvector(){ d.resize(0); };
	~SRvector(){ Free(); };
	vector <gen> d;
};

//use SRpointerVector instead of SRvector for cases like edges and faces where 
//have to conservatively allocate because only some extra pointers are wasted,
//not storage of extra class (e.g. sizeof(SRface)=100)):

template <class gen>
class SRpointerVector
{
public:
	bool isEmpty(){ return (num == 0); };

	void Free()
	{
		if (isEmpty())
			return;
		for (int i = 0; i < num; i++)
			DELETEMEMORY d[i];
		d.clear();
		num = 0;
	};
	void Free(int i)
	{
		DELETEMEMORY d[i];
		d[i] = NULL;
	};

	void Allocate(int nt)
	{
		d.resize(nt);
	};

	int GetNum(){ return num; };
	int GetNumAllocated(){ return d.size(); };
	gen* Add()
	{
		if (d[num] == NULL)
			d[num] = ALLOCATEMEMORY gen();
		num++;
		return d[num - 1];
	};
	inline gen *GetPointer(int i){return d[i]; };

	void packNulls()
	{
		if (num == 0)
			return;
		int npacked = 0;
		int i = 0;
		while (1)
		{
			if (i >= num)
				break;
			if (d[i] != NULL)
			{
				d[npacked] = d[i];
				npacked++;
				i++;
				continue;
			}
			while (1)
			{
				if (i >= num)
					break;
				if (d[i] == NULL)
				{
					i++;
					continue;
				}
				d[npacked] = d[i];
				npacked++;
				i++;
				break;
			}
		}
		num = npacked;
		d.resize(num);

	};


	SRpointerVector(){ num = 0; d.resize(0); };
	~SRpointerVector(){ Free(); };
	vector <gen *> d;
	int num;
};

class SRintMatrix
{
public:
	bool isEmpty(){ return (d == NULL); };
	void Zero()
	{
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < m; j++)
				d[i][j] = 0;
		}
	}
	void Allocate(int nt, int mt)
	{
		if(nt == 0 || mt == 0)
			return;
		if(d != NULL)
			Free();
		n = nt;
		m = mt;
		d = ALLOCATEMEMORY int* [n];
		for(int i = 0; i < n; i++)
		{
			d[i] = ALLOCATEMEMORY int[m];
		}
		Zero();
	};
	void Free()
	{
		if(d == NULL)
			return;
		for(int i = 0; i < n; i++)
			DELETEMEMORY d[i];
		DELETEMEMORY d;
		d = NULL;
		n = m = 0;
	};
	inline int Get(int i,int j)
	{
		return d[i][j];
	};
	inline void Put(int i,int j, int v)
	{
		d[i][j] = v;
	};
	inline void PlusAssign(int i, int j, int v)
	{
		d[i][j] += v;
	};
	void getSize(int& nt, int &mt){ nt = n; mt = m; };
	int getNumCols(){ return n; };
	int getNumRows(){ return m; };

	SRintMatrix(){ d = NULL; };
	SRintMatrix(int n, int m){ Allocate(n, m); };
	~SRintMatrix(){ Free(); };
private:
	int n,m;
	int** d;
};

class SRdoubleMatrix
{
public:
	bool isEmpty(){ return (d == NULL); };
	void Allocate(int nt, int mt)
	{
		if (nt == 0 || mt == 0)
			return;
		if (d != NULL)
			Free();
		n = nt;
		m = mt;
		d = ALLOCATEMEMORY double*[n];
		for (int i = 0; i < n; i++)
		{
			d[i] = ALLOCATEMEMORY double[m];
		}
		Zero();
	};
	void Free()
	{
		if (d == NULL)
			return;
		for (int i = 0; i < n; i++)
			DELETEMEMORY d[i];
		DELETEMEMORY d;
		d = NULL;
		n = m = 0;
	};
	void Zero()
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
				d[i][j] = 0.0;
		}
	};
	void Scale(int i, double s)
	{
		for (int j = 0; j < m; j++)
			d[i][j] *= s;
	}

	void Scale(int i, int j, double s)
	{
		d[i][j] *= s;
	}

	inline double Get(int i,int j)
	{
		return d[i][j];
	};
	inline void Put(int i, int j ,double v)
	{
		d[i][j] = v;
	};
	inline void PlusAssign(int i, int j, double v)
	{
		d[i][j] += v;
	};
	void getSize(int& nt, int &mt){ nt = n; mt = m; };
	int getNumCols(){ return n; };
	int getNumRows(){ return m; };
	void Copy(SRdoubleMatrix& that);
	void PlusAssign(SRdoubleMatrix& that);
	double *GetRow(int i) {	return d[i]; };

	SRdoubleMatrix(){ n = 0; m = 0; d = NULL; };
	SRdoubleMatrix(int n, int m){ Allocate(n, m); };
	~SRdoubleMatrix(){ Free(); };
private:
	int n,m;
	double** d;
};

template <class gen>
class SRmatrix
{
public:
	bool isEmpty(){ return (d == NULL); };
	void Allocate(int nt, int mt)
	{
		n = nt;
		m = mt;
		if(d != NULL)
			Free();
		d = ALLOCATEMEMORY gen* [n];
		for(int i = 0; i < n; i++)
		{
			d[i] = ALLOCATEMEMORY gen[m];
		}
	};
	void Free()
	{
		if(d == NULL)
			return;
		for (int i = 0; i < n; i++)
			DELETEMEMORY d[i];
		DELETEMEMORY d;
		d = NULL;
		n = m = 0;
	};

	SRmatrix(){ n = 0; m = 0; d = NULL; };
	SRmatrix(int n, int m){ Allocate(n, m); };
	~SRmatrix(){ Free(); };

	int n, m;
	gen** d;
};

#endif //!defined(SRUTIL_INCLUDED)


