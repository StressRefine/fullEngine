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
// SRmaterial.h: interface for the SRmaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRMATERIAL_INCLUDED)
#define SRMATERIAL_INCLUDED

enum SRmaterialType {iso, ortho, genAniso};

struct SRcij
{
	double c11;
	double c12;
	double c13;
	double c22;
	double c23;
	double c33;
	double c44;
	double c55;
	double c66;
};

class SRgenAnisoCij
{
public:
	double c[6][6];
    bool symCheck();
};
class SRstring;

class SRmaterial
{

	friend class SRinput;
	friend class SRelement;
	friend class SRthermalForce;

public:
	//functions:
	void IsoCreate(double et, double nut);
	const char* GetName();
	SRmaterialType GetType(){ return type; };
	double GetRho(){ return rho; };
	double MatScale();
	double GetAllowableStress(){ return allowableStress; };
	bool isAllowableAssigned(){ return allowableAssigned; };
	void GetIso(double& et, double& nut){ et = E; nut = nu; };
	double GetMaxSvm(){ return maxSvm; };
	void PutMaxSvm(double s);
	double GetMaxStress(){ return maxStress; };
	void PutMaxStress(double s){ maxStress = s; };
	double GetMaxStrain(){ return maxStrain; };
	void PutMaxStrain(double e){ maxStrain = e; };
	bool GetHighStressWarned() { return highStressWarned; };
	void SetHighStressWarned(bool tf = true) { highStressWarned = tf; };
	double GetVolPercentYielded();
	void SetVolPercentYielded(double p) { volPercentYielded = p; };
	void AddToVolPercentYielded(double p) { volPercentYielded += p; };
	void printToFile(SRfile &f);
	bool isActive(){ return (numElements > 0); };
	bool diffElast(SRmaterial* that);
	int getId(){ return id; };
	SRmaterial();
	~SRmaterial(){};

	int id;
	int uid;
	SRstring name;
	SRmaterialType type;
	//density:
	double rho;
	//coefficient of thermal expansion:
	double alphax;
	double alphay;
	double alphaz;
	//elastic properties:
	double E;
	double nu;
	double lambda;
	double G;
	double c11; //=lambda+2G for iso
	double tref;
	double allowableStress;
	SRcij orthoCij;
	SRgenAnisoCij genAnisoCij;
	bool allowableAssigned;
	bool highStressWarned;
	double maxSvm;
	double maxStress; //max vm or in any component
	double maxStrain;
	double svmavg;
	double volPercentYielded;
	int numElements;
	int numelyielded;
};

class SRElProperty
{
public:
	int uid;
	int matid;
	int matuid;
};


#endif //!defined(SRMATERIAL_INCLUDED)
