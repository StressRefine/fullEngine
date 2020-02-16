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
// SRpostProcess.h: interface for the SRpostProcess class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRPOSTPROCESS_INCLUDED)
#define SRPOSTPROCESS_INCLUDED

class SRfile;
class SRvec3;
class SRBreakoutData;

enum SRentityType;


enum SRstressComponent { xxComponent, yyComponent, zzComponent, xyComponent, xzComponent, yzComponent };

class SRpostStressVector
{
public:
	SRpostStressVector(){ for (int i = 0; i < 6; i++) stress[i] = 0.0; };
	void PlusAssign(double stresst[])
	{
		for (int i = 0; i < 6; i++)
			stress[i] += stresst[i];
	};
	void Assign(double stresst[])
	{
		for (int i = 0; i < 6; i++)
			stress[i] = stresst[i];
	};
	void Scale(double s)
	{
		for (int i = 0; i < 6; i++)
			stress[i] *= s;
	};
	double* Getstress(){ return stress; };
private:
	double stress[6];
};

struct SRElemWithRadius
{
	int elId;
	double Rad;
};

class SRsaveBreakoutSphere;
class SRpostProcess
{
	friend class SRinput;
public:
	SRpostProcess();
	void PostProcess();
	void CalculateMaxStress();
	double StressMaxCheck(SRelement* elem, SRvec3& pos, double stress[]);
	void StressMaxCheckvsAllowable(SRelement* elem, double svmMax);
	void CalculateElementStresses(SRelement* elem, bool checkMaxOnly = false);
	void fillSacricialElementNodalStress(SRelement* elem, bool doMaxClipping);
	void PostProcessElementStresses();
	void PostProcessForBreakout();
	void PostProcessForBreakoutPartial();
	double ElementStrainSmoothByMaterial();
	double GlobalStrainSmooth();
	void Cleanup();
	double *getSmoothedStrainVec(int i){ return smoothedStrains[i].GetVector(); };
	int getElSmooth(int i){ return elSmooth.Get(i); };
	void readPreviousDisplacementForBreakout();
	bool findSacrificialElementsAtKinks();
	int findSaveBreakoutElementsTopo(int nel, int& nelAllNodesDisp, bool& topoIsRagged);
	void saveBreakoutModel();
	void printNodalForce(int nodeUid, SRforce* force, SRfile&f);
	double GetMaxFaceRotation();
	int MeshToFrd();
	void findElidAtRad0andSetSB(SRBreakoutData *bdat);
	void findElidAtRad0AndMaxRadNearOrigin(SRBreakoutData *bdat);
	int autoBreakoutSphere(SRBreakoutData *bdat, int NumElBreakoutCandidate);
	void autoLocalAdapt();
	void PFringeOut();
	void writeFrd();
	void outputBreakout(int nelBreakout, int nbrcon);
	int fillBreakoutModelOriginFromDispCentroid(SRBreakoutData *bdat);
	int fillBreakoutModelOriginFromHCentroidStresses(SRBreakoutData *bdat, bool SBonlyNoBsurfNoSacr = false);
	int fillBreakoutModelOriginFromHCentroidStressesElemRadList(SRBreakoutData *bdat, int nelToCheck);
	void PlotSBElems(char* inName = NULL, int filenum = -1);
	void PlotElems(int numel, int* elems, char* inName = NULL, int filenum = -1);
	void PlotElemsSBOnly(int numel, int* elems, char* inName = NULL, int filenum = -1);
	void PlotFaces(int numFace, SRface* facev[]);
	void PlotFaces(int numFace, int fid[]);
	void PlotEdges(int numEdge, int edgeId[]);
	void PlotSacr();
	int topoFilterSaveBreakoutElems(bool &isRagged);
	void checkElemFacesForPartBdry(SRelement* elem, int& nextOutEl, SRintVector& outchecklist, int numLevels, int& numElOnTopoList);
	void SortElemsByRadius();
	void PlotTopoLevelElems(int numLevels);
	void PlotElemsByRadius(double Rad, int filenum);
	void writeClippedPrevFrd(int nelClipped);

	SRvector <SRElemWithRadius> elemRadiusList;
	SRdoubleVector smoothedStrains[6];
	SRintVector elemOnTopoList;
	SRintVector elSmooth;
	SRintVector nodalStressCount;
	SRdoubleMatrix nodalStress;
	bool outputCgx;
	SRintVector nodeElementOwnerSave;
	SRvector <SRvec3> nodeDisps;
	double svmAvg;
	SRintVector nodalPOrders;
	double maxRadNearOrigin;
	double maxRadSaveBreakout;
	int nelBreakoutMax;
	int minElForBreakout;
	int elIdAtRad0;
	SRpointerVector <SRconstraint> breakoutConstraints;
	SRintVector breakoutElemUids;
	bool needTopoFilter;
	bool noTopoFilter;
	int nelAllNodesDisp;
	int nelNearOrigin;
	SRintVector bktElSacrList;

	//f06 input:
	bool readPreviousDisplacementF06();
	void readPreviousResultsF06();
	void readPreviousStressF06();
	bool getF06DispOrStressLineCheckForSkip(SRfile& f, SRstring& line);

	//f06 output:
	void OutputF06();
	void OutputF06UseElementGetStress();
	void pageCheckF06(int& nlines, bool disp = false);
	void writeHeaderF06();
	void writeDispHeaderF06();
	void clipBreakoutSacrStreses();
};

#endif // !defined(SRPOSTPROCESS_INCLUDED)
