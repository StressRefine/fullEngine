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
// SRmaterial.cpp: implementation of the SRmaterial class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SRutil.h"
#include "SRmodel.h"

SRmaterial::SRmaterial()
{
	type = iso;
	E = 0.0;
	nu = 0.0;
	rho = 0.0;
	tref = 0.0;
	allowableStress = 0.0;
	allowableAssigned = false;
	highStressWarned = false;
	maxSvm = 0.0;
	maxStress = 0.0;
	maxStrain = 0.0;
	volPercentYielded = 0.0;
	numElements = 0;
	numelyielded = 0;
};

double SRmaterial::GetVolPercentYielded()
{
	if (numElements == 0)
		return 0.0;
	else
		return volPercentYielded / numElements;
}


void SRmaterial::PutMaxSvm(double s)
{
	maxSvm = s; 
	if (s > maxStress)
		maxStress = s;
};


void SRmaterial::IsoCreate(double et, double nut)
{
	//create material of type iso
	//input:
		//et = Young's modulus
		//nut = Poisson's ratio
	//note:
		//fills class variables E, nu, G, lambda, abd c11
	type = iso;
	E = et;
	nu = nut;
	G = E / (2.0*(1.0 + nu));
	lambda = 2.0*G*nu / (1.0 - 2.0*nu);
	c11 = lambda + 2.0*G;
};

double SRmaterial::MatScale()
{
	//determine a characteristic scale for this material
	//return:
		//scale = Young's modulus for isotropic element or equivalent value for anisotropic
	if (type == iso)
		return E;
	else if (type == ortho)
	{
		double pseudolam, pseudoG, PseudoGnu, lg;
		SRcij& cij = orthoCij;
		pseudolam = ONETHIRD*(cij.c12 + cij.c13 + cij.c23);
		pseudoG = ONETHIRD*(cij.c44 + cij.c55 + cij.c66);
		lg = pseudolam / pseudoG;
		PseudoGnu = 0.5*lg / (1.0 + lg);
		return pseudoG*2.0*(1.0 + PseudoGnu);
	}
	else if (type == genAniso)
	{
		double pseudolam, pseudoG, PseudoGnu, lg;
		SRgenAnisoCij& gcij = genAnisoCij;
		pseudolam = (gcij.c[0][1] + gcij.c[0][2] + gcij.c[1][2]);
		pseudoG = ONETHIRD*(gcij.c[3][3] + gcij.c[4][4] + gcij.c[5][5]);
		lg = pseudolam / pseudoG;
		PseudoGnu = 0.5*lg / (1.0 + lg);
		return pseudoG*2.0*(1.0 + PseudoGnu);
	}
	ERROREXIT;
	return 0;
}

void SRmaterial::printToFile(SRfile &f)
{
	f.Print(" %s", name.str);
	if (type == iso)
	{
		f.Print(" isotropic");
		f.PrintReturn();
		f.PrintLine(" %lg %lg %lg %lg //rho, alpha tref allowable", rho, alphax, tref, allowableStress);
		f.PrintLine(" %lg %lg //E, nu", E, nu);
	}
	else if (type == ortho)
	{
		f.Print(" orthotropic");
		f.PrintReturn();
		f.PrintLine(" %lg %lg %lg %lg %lg %lg //rho, alphax, alphay, alphaz tref allowable", rho, alphax, alphay, alphaz);
		f.PrintLine(" %lg %lg %lg %lg %lg %lg %lg %lg %lg //c11,c12,c13,c22,c23,c33,c44,c55,c6",
			orthoCij.c11, orthoCij.c12, orthoCij.c13,
			orthoCij.c22, orthoCij.c23, orthoCij.c33,
			orthoCij.c44, orthoCij.c55, orthoCij.c66);
	}
	else if (type == ortho)
	{
		f.Print(" general");
		f.PrintReturn();
		f.PrintLine(" %lg %lg %lg %lg %lg %lg //rho, alphax, alphay, alphaz tref allowable", rho, alphax, alphay, alphaz);
		for (int row = 0; row < 6; row++)
		{
			for (int col = 0; col < 6; col++)
				f.Print(" %lg", genAnisoCij.c[row][col]);
			f.PrintReturn();
		}
	}
}

bool SRmaterial::diffElast(SRmaterial* that)
{
	if (type == iso)
	{
		if (that->type != iso)
			return true;
		if (fabs(E - that->E) > TINY*E)
			return true;
		if (fabs(nu - that->nu) > TINY*E)
			return true;
	}
	else if (type == ortho)
	{
		if (that->type != ortho)
			return true;
		double scale = orthoCij.c11;
		if (fabs(orthoCij.c11 - that->orthoCij.c11) > TINY*scale)
			return true;
		if (fabs(orthoCij.c12 - that->orthoCij.c12) > TINY*scale)
			return true;
		if (fabs(orthoCij.c13 - that->orthoCij.c13) > TINY*scale)
			return true;
		if (fabs(orthoCij.c22 - that->orthoCij.c22) > TINY*scale)
			return true;
		if (fabs(orthoCij.c23 - that->orthoCij.c23) > TINY*scale)
			return true;
		if (fabs(orthoCij.c44 - that->orthoCij.c44) > TINY*scale)
			return true;
		if (fabs(orthoCij.c55 - that->orthoCij.c55) > TINY*scale)
			return true;
		if (fabs(orthoCij.c66 - that->orthoCij.c66) > TINY*scale)
			return true;
	}
	else if (type == genAniso)
	{
		if (that->type != genAniso)
			return true;
		double scale = genAnisoCij.c[0][0];
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				if (fabs(genAnisoCij.c[i][j] - that->genAnisoCij.c[i][j]) > TINY*scale)
					return true;
			}
		}
	}
	return false;
}

bool SRgenAnisoCij::symCheck()
{
	//check if this material is symmetric
    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < i; j++)
        {
            double diff = fabs(c[i][j] - c[j][i]);
            if(diff > TINY)
                return false;
        }
    }
    return true;
}

