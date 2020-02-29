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
// SRcoord.h: interface for the SRcoord class: Coordinate system definitions
//

#if !defined(SRCOORD_INCLUDED)
#define SRCOORD_INCLUDED

enum SRcoordType {cartesian,spherical,cylindrical};

class SRstring;
class SRcoord  
{
	friend class SRinput;

public:
	void Create(double x0, double y0, double z0, SRvec3 p13, SRvec3 p3);
	void Create(double x0, double y0, double z0);
	void CalculateBasisVectors(SRvec3& p, SRvec3 &e1l, SRvec3 &e2l, SRvec3 &e3l);
	void CalculateLocalDirection(SRvec3& p, int dof, SRvec3 &el);
	void GetRotationMatrix(bool toLocal, SRvec3& p, SRmat33& R);
	void Copy(SRcoord& c2);
	void operator =(SRcoord& c2){ Copy(c2); };
	SRcoordType GetType(){ return type; };
	const char* GetName();
	void PrintToFile(SRfile& f);
	void GetPos(double &x, double &y, double &z, SRvec3& pos);
	void VecTransform(SRvec3 p, SRvec3 &v);
	int checkParallelToGcs(int dof);

	SRcoord();
	//data:
private:
	SRstring name;
	SRcoordType type;
	SRvec3 origin;
	SRvec3 e1, e2, e3;
	bool gcsAligned;
	int uid;
	int otherCoordid;
	bool gcsaligned;
};

#endif // !defined(SRCOORD_INCLUDED)
