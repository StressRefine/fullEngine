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
// SRmodel.h: interface for the SRmodel class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRMODEL_INCLUDED)
#define SRMODEL_INCLUDED

#define MAXADAPTLOOPS 3
#define SKIPCONTINATION true
#define MAXNODEFACEOWNERS 101
#define ERRTOL 0.02
#define MINP 2
#define MAXP 8
#define THINEDGEMAXP 2


#include "SRmachDep.h"
#include "SRutil.h"
#include "SRmath.h"
#include "SRstring.h"
#include "SRfile.h"
#include "SRsolver.h"
#include "SRcoord.h"
#include "SRmaterial.h"
#include "SRconstraint.h"
#include "SRforce.h"
#include "SRmap.h"
#include "SRnode.h"
#include "SRedge.h"
#include "SRface.h"
#include "SRbasis.h"
#include "SRelement.h"
#include "SRerrorCheck.h"
#include "SRinput.h"
#include "SRpostProcess.h"
#include "SRmklUtil.h"
#include "SRpoundDefines.h"

#define ERROREXIT SRutil::ErrorExit(__FILE__,__LINE__)

enum SRentityType{ nodeType, edgeType, faceType, elementType };


class SRElementData
{
	friend class SRmodel;

public:
	SRdoubleMatrix& Getdbdx(){ return dbdx; }
	SRdoubleMatrix& Getdbdy(){ return dbdy; }
	SRdoubleMatrix& Getdbdz(){ return dbdz; }
	double* GetIntwt(){ return intWt.d; };
	double* GetBasisVec(){ return basisVec.d; };
	double* Getdbasisdr(){ return dbasisdr.d; };
	double* Getdbasisds(){ return dbasisds.d; };
	double* Getdbasisdt(){ return dbasisdt.d; };
	double* GetEnfdForceVec(){ return enfdForceVec.d; };
private:
	SRdoubleMatrix dbdx, dbdy, dbdz;
	SRdoubleVector basisVec;
	SRdoubleVector dbasisdr;
	SRdoubleVector dbasisds;
	SRdoubleVector dbasisdt;
	SRdoubleVector intWt;
	SRdoubleVector enfdForceVec;
	SRdoubleVector elementStiffness;
};

class SRBreakoutData
{
public:
	SRBreakoutData()
	{
		rad = 0.0;
		atMax = true;
		bktNode = -1;
		fromPartialDispFile = false;
	};
	SRfile saveBreakoutMshFile;
	double rad;
	SRvec3 origin;
	int bktNode;
	bool atMax;
	bool fromPartialDispFile;
};

class SRPassData
{
public:
	SRPassData()
	{
		maxp = 0;
		neq = 0;
		err = 0.0;
		maxsvm = 0.0;
		maxsp1 = 0.0;
		maxCustom = 0;
	};
	int maxp;
	int neq;
	double err;
	double maxsvm;
	double maxsp1;
	double maxCustom;
};

class SRmodel
{
	friend class SRinput;

public:

	SRmodel();
	void Initialize();

	void CleanUp(bool partial = false);
	void NumberEquationsSmooth();
	bool AdaptUniform(int pIteration, bool checkErrorOnly);
	double GetEnforcedDisplacementCoeff(int gfun, int dof);
	double GetDisplacementCoeff(int gfun, int dof);
	double GetDisplacementCoeff(int eq) { return solution.Get(eq); };
	void NumberEquations();
	void AllocateDofVectors();
	void Create();
	void CalculateElementStiffnesses(bool anyLcsEnfd);
	bool CalculateElementStiffnessesOmp();
	void NumberGlobalFunctions();
	bool Adapt(int pIteration, bool checkErrorOnly = false);
	void Run();
	bool PreProcessPenaltyConstraints();
	void ProcessConstraints();
	void EnforcedDisplacementAssemble(bool doingPrevSolution = false);
	void ProcessForces();
	void ProcessVolumeForces(SRvec3& ResF);
	void ProcessElemVolumeForce(SRvolumeForce* vf, SRelement*elem, SRvec3& ResF, SRdoubleVector& basvec, double* globalForce);
	void NodalMaxDisp(int& nodeUidAtMaxDisp, SRvec3& disp);
	void PrintLoadVector();
	void openElStiff();
	void closeElStiff();
	void WriteElementStiffness(SRelement* elem, int len, double* stiff);
	double* ReadElementStiffness(SRelement* elem, bool readPrevP = false, int processorNum = 0);
	SRfile* GetElementStiffnessFile(bool readPrevP = false);
	SRnode* addNode(){ return nodes.Add(); };
	SRconstraint* addConstraint(){ return constraints.Add(); };
	void resizeConstraints(int n){ constraints.Allocate(n); };
	void AddSoftSprings();
	void UpdateCustomCriterion(double stress[]);


	int GetNumNodes(){ return nodes.GetNum(); };
	SRnode* GetNode(int i){ return nodes.GetPointer(i); };

	int GetNumEdges(){ return edges.GetNum(); };
	SRedge* GetEdge(int i){ return edges.GetPointer(i); };

	int GetNumFaces(){ return faces.GetNum(); };
	SRface* GetFace(int i){ return faces.GetPointer(i); };

	int GetNumElements() { return elements.GetNum(); };
	SRelement* GetElement(int i){ return elements.GetPointer(i); };

	int GetNumMaterials() { return materials.GetNum(); };
	SRmaterial* GetMaterial(int i){ return materials.GetPointer(i); };

	int GetNumConstraints(){ return constraints.GetNum(); };
	SRconstraint* GetConstraint(int i){ return constraints.GetPointer(i); };

	int GetNumForces(){ return forces.GetNum(); };
	SRforce* GetForce(int i){ return forces.GetPointer(i); };

	int GetNumCoords(){ return Coords.GetNum(); };
	SRcoord* GetCoord(int i){ return Coords.GetPointer(i); };

	SRthermalForce* GetThermalForce() { return thermalForce; };

	double* GetElementStiffnessVector(int processorNum = 0);

	double* GetSolutionVector() { return solution.GetVector(); };
	void AddToSolutionVector(int i, double d){ solution.PlusAssign(i, d); };
	void CopyToSolutionVector(int i, double d){ solution.Put(i, d); };
	void CopyToSolutionVector(SRdoubleVector& v){ solution.Copy(v); };

	void PutFunctionEquation(int fun, int dof, int id){ functionEquations.Put(fun, dof, id); };
	int GetFunctionEquation(int fun, int dof){ return functionEquations.Get(fun, dof); };

	int GetSmoothFunctionEquation(int fun){ return smoothFunctionEquations.Get(fun); };
	void AllocateSmoothFunctionEquations(int n){ smoothFunctionEquations.Allocate(n); };

	void PutEnforcedDisplacement(int fun, int dof, double val){ enforcedDisplacements.Put(fun, dof, val); };
	double GetEnforcedDisplacement(int fun, int dof) { return enforcedDisplacements.Get(fun, dof); };

	void AllocateSkipFun(int n){ skipFun.Allocate(n); };

	void CreateElemEdges(SRelement *elem, int nnodes, int inputNodes[]);
	void CreateElemEdgesSBO(SRelement *elem);
	void FillGlobalFaces(bool needGlobalNodeOrderForBreakout = false);
	void FindElemsAdjacentToBreakout();
	void setElemsAllNodesHaveDisp();

	//"read" routines for private data:
	double GetSize(){ return size; };
	double GetErrorTolerance(){ return ErrorTolerance; };
	double GetStressMax(){ return stressMax; };
	double GetStressMaxComp(int i){ return stressMaxComp[i]; };
	int GetNumFunctions(){ return numFunctions; }
	int GetNumEquations(){ return numEquations; }
	int GetNumSmoothEquations(){ return numSmoothEquations; }
	bool IsAnyEnforcedDisplacement(){ return anyEnforcedDisplacement; };
	int GetMaxPorder(){ return maxPorder; };
	void SetMaxPorderLowStress(int p){ maxPorderLowStress = p; };
	int GetMaxPorderLowStress(){ return maxPorderLowStress; };
	int GetMinPorder(){ return minPorder; };
	int GetThinEdgeMaxPorder(){ return thinEdgeMaxPorder; };
	int GetmaxNumElementFunctions(){ return maxNumElementFunctions; };
	int GetnumEquations(){ return numEquations; };
	int GetnumSmoothEquations(){ return numSmoothEquations; };
	int GetMaxNumThreads() { return maxNumThreads; };
	bool isBreakout(){ return breakout; };
	bool areElementsInMemory(){ return elementsInMemory; };
	int GetMaxNodeUid(){ return maxNodeUid; };
	int GetMaxElemUid(){ return maxElemUid; };
	bool GetFunUncon(int fun){ return (funUncon.Get(fun) == 1); };
	SRvec3* GetMaxStressPos(){ return &maxStressPos; };
	bool checkOrphanNode(int id){ return GetNode(id)->isOrphan(); };
	SRnode* GetNodeFromUid(int i);
	void SetNumThreads();

	//scalar set routines:
	void SetErrorMax(double error, int eluid, double errorSmoothRaw, double errorFaceJumps);
	void SetStressMax(SRelement* elem, SRvec3& pos, double svm);
	void SetStressMaxComp(double *stressComp);
	double GetStrainMax(){ return strainMax; };
	void SetStrainMax(double e){ strainMax = e; };
	
	void updateVolume(double elVol){ volume += elVol; };
	void zeroStressMax();
	void setMaxPorder(int maxp){ maxPorderFinalAdapt = maxp; };
	void setMinPorder(int minp);
	void setErrorTolerance(double tol){ ErrorTolerance = tol; };
	void allocateElementData();
	void FreeElementData();
	void allocateSmallElementData(bool anyLcsEnfd = false);
	void FreeSmallElementData();
	void PrintStresses(double *stress);
	void SetDetectSacr(bool sacr){ detectSacrificialElements = sacr; };
	void SetFreezePorder(int p){ freezeSacrificalP = p; };
	int GetFreezePorder(){ return freezeSacrificalP; };
	void setAdaptLoopMax(int nits);
	void PStats();
	void checkElementMapping();
	bool checkFlattenedElementHighStress();
	bool FreezeSacrificialEdges(){ return freezeSacrEdges; };
	void SetFreezeSacrificialEdges(){ freezeSacrEdges = true; };
	void SetsimpleElements(){ simpleElements = true; };
	bool UseSimpleElements(){ return simpleElements; };
	bool useParallelElements(){ return maxNumThreads > 1; };
	void SetsimplePardiso(){ simplePardiso = true; };
	bool UseSimplePardiso(){ return simplePardiso; };
	void SetDetectThinElements(){ detectThinElements = true; };
	void CheckElementsFitInMemory();
	void WriteAllElements();
	void dumpFileOpen(bool append = false);
	void dumpFileClose(){ dumpFile.Close(); };
	void SetMaxNodeUid(int uid){ maxNodeUid = uid; };
	int GetMaxPinModel() {return maxPinModel;};
	void SetEdgesToPorder(int p);
	bool GetAnyWedges(){ return anywedges; };
	bool GetAnyBricks(){ return anybricks; };
	void GetFileNameFromExtension(char* ext, SRstring& name);
	bool SingStressCheck();
	bool readResultsSrr();
	double findMaxStressNearBreakoutOrigin();
	void checkForHotSpotElemsNearBreakoutBoundary();

	SRElementData* GetElementData(int i);

	//these should be set during input or from SRmodel routines and only accessed readonly elsewhere through query functions
	double size;
	double ErrorTolerance;
	int maxPinModel;
	int adaptLoopMax;
	bool pOrderUniform;
	bool anyEnforcedDisplacement;
	int numFunctions;
	int maxNumElementFunctions;
	int numEquations;
	int numSmoothEquations;
	int numprocessors;
	bool detectSacrificialElements;
	bool detectThinElements;
	int freezeSacrificalP;
	bool elementsInMemory;
	int maxNumThreads;
	bool freezeSacrEdges;
	bool simpleElements;
	bool simplePardiso;
	bool breakout;
	bool breakoutByMat;
	bool breakoutAtNode;
	double feasvmmax;
	SRvec3 feavmmaxpos;
	int numShellOrBeamNodes;
	bool partialFlatten;
	int maxStressElid;
	bool needSoftSprings;
	double MaxElementMem;
	double errorMax;
	double errorSmoothRawAtMax;
	double errorFaceJumpAtMax;
	int elUidAtMaxError;
	double volume;
	double strainMax;
	double stressMaxComp[6];
	int minPorder;
	int maxPorder;
	int maxPorderFinalAdapt;
	int maxPorderLowStress;
	int thinEdgeMaxPorder;
	int maxNodeUid;
	int maxElemUid;
	int adaptIt;
	SRvec3 maxStressPos;
	int nodeidAtMaxStressPos;
	int maxPJump;
	double stressUnitConversion;
	SRstring stressUnitstr;
	double lengthUnitConversion;
	SRstring lengthUnitstr;
	bool anyUnsupportedFace;
	bool ReadDispStressSRR;
	bool useUnits;
	SRintVector faceBktConNodes;

	SRpointerVector <SRnode> nodes;
	SRpointerVector <SRedge> edges;
	SRpointerVector <SRface> faces;
	SRpointerVector <SRconstraint> constraints;
	SRpointerVector <SRcoord> Coords;
	SRpointerVector <SRmaterial> materials;
	SRpointerVector <SRforce> forces;
	SRpointerVector <SRvolumeForce> volumeForces;
	SRvector <SRelement> elements;
	SRthermalForce* thermalForce;
	SRpointerVector <SRFaceForceGroup> faceForceGroups;
	SRpointerVector <SRElProperty> elProps;

	SRdoubleVector solution;
	SRintMatrix functionEquations;
	SRintVector smoothFunctionEquations;
	SRdoubleMatrix enforcedDisplacements;
	SRintVector skipFun;
	SRintVector funUncon;
	SRvector <SRPassData> passData;

	//scratch space needed by elements:
	SRvector <SRElementData> elData;

	char stringBuffer[MAXLINELENGTH];
	SRfile elementFileEven;
	SRfile elementFileOdd;
	bool anybricks;
	bool anytets;
	bool anywedges;
	bool allConstraintsAsPenalty;
	SRvec3 umaxInModel;
	int nodeuidAtumax;
	bool doEnergySmooth;

	//1 instance of each utility class:
	SRbasis basis;
	SRmath math;
	SRmap map;
	SRsolver solver;
	SRpostProcess post;
	SRerrorCheck errorChecker;
	SRinput input;

	SRstring wkdir;
	SRstring outdir;
	SRfile outputFile;
	SRfile inputFile;
	SRfile repFile;
	SRfile logFile;
	SRfile sumFile;
	SRstring fileNameTail;
	SRfile cgxFrdFile;
	SRfile settingsFile;
	SRfile f06File;
	SRfile f06outFile;
	SRfile dumpFile; //file for dumping diffs for debugging

	double maxFlattened;
	int maxFlattenedElUId;

	SRBreakoutData saveBreakoutData;

	bool saveBreakout;
	bool allMatsHaveAllowable;
	double maxAllowableAnyActiveMat;
	int numFlattenedElHighStress;
	bool allMatsEquivElast;
	bool anyMatYielded;
	double prevStressMax;
	double stressMax;
	double stressMaxP2;
	SRstring breakoutMat;
	double maxsp1;
	double minsp2;
	double maxCustom;
	bool checkForHotSpotatMax;
	bool anyHotSpotElFound;
	double breakoutRadSacr;
	bool localAdapt;
	bool allFacesFlat;
	bool anyUnsupNode;
	bool anyBsurf;
};

class SRFaceForceGroup
{
public:
	SRFaceForceGroup(){ nfaceFunMax = 0; fromNodalForces = false; };
	void addFace(SRface* face)
	{
		faceIds.PushBack(face->id);
		if (face->GetNumNodesTotal() > nfaceFunMax)
			nfaceFunMax = face->GetNumNodesTotal();
	};

	SRintVector faceIds;
	SRintVector nodeIds;
	int nfaceFunMax;
	SRdoubleMatrix ForceDof;
	SRdoubleVector smoothStiff;
	SRintMatrix faceFunLoc;
	bool fromNodalForces;
};


#endif //!defined(SRMODEL_INCLUDED)
