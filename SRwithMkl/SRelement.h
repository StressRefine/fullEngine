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
// SRelement.h: interface for the SRelement class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRELEMENT_INCLUDED)
#define SRELEMENT_INCLUDED

#include "SRbasis.h"

enum SRelementType { tet, wedge, brick};

class SRlocalEdge;

class SRlocalFace
{
	friend class SRelement;
public:
	SRface* GetFace();

	SRlocalFace()
	{ 
		globalFaceId = -1;
		for (int i = 0; i < 4; i++)
			globalNodeOrder[i] = -1;
	};
	int GetGlobalFaceId(){ return globalFaceId; };
	void PutGlobalFaceId(int f) { globalFaceId = f; };
	void PutGlobalNodeOrder(int gn, int ln){ globalNodeOrder[gn] = ln; };
private:
	int globalFaceId;
	int globalNodeOrder[4];
};

class SRBTC
{
public:
	SRBTC();
	double BTC11, BTC12, BTC13, BTC14, BTC15, BTC16;
	double BTC21, BTC22, BTC23, BTC24, BTC25, BTC26;
	double BTC31, BTC32, BTC33, BTC34, BTC35, BTC36;
};


class SRelement
{
	friend class SRinput;

public:
	SRelement();
	~SRelement(){ Cleanup(); };
	void FillMappingNodes();
	void FillMappingNodesBreakoutOnly(bool needShapeDerivs = false);
	void DownloadSmoothedStrains();
	void CalibratePenalty();
	void GetStress(double r, double s, double t, double stress[], bool checkMaxStrain = false);
	void GetSmoothedStrain(double r, double s, double t, double strain[]);
	int MapFaceToElFuns(int lface, int FaceToElFuns[]);
	int MapEdgeToElFuns(int ledge, int EdgeToElFuns[]);
	int MapEdgeP2FnToElFun(int ledge);
	double CalculateRawStrain(double r, double s, double t, double strain[], double& etx, double& ety, double& etz);
	int GetMaxPorder();
	void DownloadDisplacement();
	void GetBrickNodePos(int lnode, double& r, double& s, double& t);
	void GetWedgeNodePos(int lnode, double& r, double& s, double& t);
	void GetBrickFaceNodes(int lface, int& n1, int& n2, int& n3, int& n4);
	void GetWedgeFaceNodes(int lface, int& n1, int& n2, int& n3, int& n4);
	void Clear();
	void NodeNaturalCoords(int lnode, double& r, double& s, double& t);
	void StraintoStress(double r, double s, double t, double strain[], double stress[]);
	void AddPenaltyToStiffnessandEnfDisp(double* stiff);
	SRnode* GetNode(int localnodenum);
	int localEdgeMatch(int gId);
	SRnode* GetLocalEdgeNode(int lej, int localnodenum){ return localEdges.Get(lej).GetEdge()->GetNode(localnodenum); };
	void FillKel33GenAnisoRowWise(SRBTC& btc, int rowfun, int colfun, double dbdx, double dbdy, double dbdz, double kel33[3][3]);
	void FillBTC(double intwt, double dbdx, double dbdy, double dbdz, SRBTC& btc);
	double* CalculateStiffnessMatrix(int processorNum, int &len);
	void colFunLoopIsoOrtho(double* dbdxv, double* dbdyv, double* dbdzv, double w, int rowfun, SRBTC& btc, double *stiff);
	void colFunLoopGenAniso(double* dbdxv, double* dbdyv, double* dbdzv, double w, int rowfun, SRBTC& btc, double *stiff);
	SRface* GetFace(int localfacenum);
	SRedge* GetEdge(int localedgenum);
	void GetFaceNodes(int lface, int& n1, int& n2, int& n3, int& n4);
	void GetEdgeLocalNodeIds(int lej, int& ln0, int& ln1);
	void GetEdgeNodeIds(int lej, int& n1, int& n2);
	void GetQuadFaceGlobalCoords(int lface, double rfl, double sfl, double& rf, double& sf);
	void GetQuadFaceDerivs(int lface, double rfl, double sfl, double& drfdrfl, double& drfdsfl, double& dsfdrfl, double& dsfdsfl);
	void Create(int userid, int nnodes, int nodest[], SRmaterial* mat);
	void Create(SRelementType typet, int userid, int nnodes, int nodest[], SRmaterial* mat);
	int GetUserid(){ return userId; };
	int GetId(){ return id; };
	SRelementType GetType(){ return type; };
	int GetNumNodes(){ return nodeIds.GetNum(); };
	int GetNumNodesTotal(){ return nodeIds.GetNum() + localEdges.GetNum(); };
	int GetNodeId(int i){ return nodeIds.Get(i); };
	int GetNumLocalEdges(){ return localEdges.GetNum(); };
	int GetLocalEdgeGlobalId(int i){ return localEdges.Get(i).GetGlobalEdgeId(); };
	int GetLocalEdgePOrder(int i){ return localEdges.Get(i).GetPOrder(); };
	int GetLocalEdgeMidNodeId(int i){ return localEdges.Get(i).GetMidNodeId(); };
	int GetLocalEdgeNodeId(int lej, int lnode);
	int GetLocalEdgeDirection(int lej){ return localEdges.Get(lej).GetDirection(); };
	int GetNumLocalFaces(){ return localFaces.GetNum(); };
	int GetLocalFaceGlobalId(int i){ return localFaces.Get(i).GetGlobalFaceId(); };
	int GetLocalFaceGlobalNodeOrder(int lface, int lnode){ return localFaces.Get(lface).globalNodeOrder[lnode]; };
	int GetMaterialId(){ return elMat->id; };
	int GetNumFunctions() { return globalFunctionNumbers.GetNum(); };
	int GetFunctionNumber(int i){ return globalFunctionNumbers.Get(i); };
	void AllocateRawStrains(int n){ rawStrains.Allocate(n, 6); };
	void FreeRawStrains();
	void PutRawStrain(int i, double *strain);
	double GetRawStrain(int i, int e){ return rawStrains.Get(i, e); };
	SRvec3 GetDisp(double r, double s, double t);
	bool isSacrificial(){ return (sacrificial != 0); };
	void SetSacrificial(int pf = 1);
	int GetSacrificialPorder(){ return sacrificial; };
	double GetMaxStress(){ return maxStress; };
	double GetSvmMax(){ return svmMax; };
	void SetSvmMax(double s){ svmMax = s; };
	double GetPenaltyStiffness(){ return penaltyStiffness; };
	void NodeLocalCopy(int lnode, int gnode);
	double FillMapping(double r, double s, double t, bool detJonly = false);
	void FillJacobian(double r, double s, double t, bool FillElJac = true);
	void FillTetJacobian(double r, double s, double t, bool FillElJac = true);
	void FillWedgeJacobian(double r, double s, double t, bool FillElJac = true);
	void FillBrickJacobian(double r, double s, double t, bool FillElJac = true);
	void BrickShapeDerivs(double r, double s, double t);
	void FillJac();
	void Position(double r, double s, double t, SRvec3 &p);
	void PositionLinear(double r, double s, double t, SRvec3 &p);
	int ShapeFunctions(double r, double s, double t, double N[]);
	int ShapeFunctionsLinear(double r, double s, double t, double N[]);
	int BrickShapeFns(double r, double s, double t, double N[]);
	int WedgeShapeFns(double r, double s, double t, double N[]);
	void XyzDerivatives(double dfdr, double dfds, double dfdt, double& dfdx, double& dfdy, double& dfdz);
	bool hasLcsConstraint(){ return (nodeLcSConstraints.GetNum() != 0 || faceLCSConstraints.GetNum() != 0); };
	void AddfaceLCSConstraint(int lface){ faceLCSConstraints.pushBack(lface); };
	void AllocateGlobalFunctionNumbers(int n){ globalFunctionNumbers.Allocate(n); };
	void PutGlobalFunctionNumbers(int lfun, int gfun){ globalFunctionNumbers.Put(lfun, gfun); };
	int GetNewPorder(){ return newPorder; };
	void PutNewPorder(int p){ newPorder = p; };
	void AssignLocalEdge(int lej, int gej, int direction){ localEdges.GetPointer(lej)->Assign(gej, direction); };
	SRlocalFace* GetLocalFace(int lface){ return localFaces.GetPointer(lface); };
	void AddNodeLcSConstraints(int l){ nodeLcSConstraints.pushBack(l); };
	double* FillBasisFuncs(double r, double s, double t, SRBasisCallType calltype);
	void Cleanup();
	int GetStiffLength(){ return stiffLength; };
	int GetFilePos(){ return filePos; };
	void SetFilePos(int pos){ filePos = pos; };
	double *GetStiffnessMatrix(){ return stiffnessMatrix.GetVector(); };
	void SetPChanged(bool tf){ pChanged = tf; };
	bool GetPChanged(){ return pChanged; };
	bool DetectThinEdges();
	void FillStiffDiag(int neq);
	int GetStiffnessLocation(int row, int col);
	int GetStiffnessLocationAboveDiag(int row, int col);
	SRmaterial* GetMaterial(){ return elMat; };
	void SetMaterial(SRmaterial* mat){ elMat = mat; };
	bool testMapping(bool okToFail = true);

	bool checkForStraightenedEdges();
	void rescaleFlattenFraction();
	double GetSize(){ return size; };
	double GetVolume(){ return elVol; };
	double getdNdr(int i){ return dNdr[i]; };
	double getdNds(int i){ return dNds[i]; };
	double getdNdt(int i){ return dNdt[i]; };
	void SetBasisData(int processorNum = 0);
	bool isConstrained(){ return constrained; };
	void SetConstrained(){ constrained = true; };
	void deleteStiffnessMaxtrix(){ stiffnessMatrix.Free(); };
	SRmat33* GetElJac(){ return &elJac; };
	void PutThread(int t){ thread = t; };
	int GetThread(){ return thread; };
	void PutError(double e){ error = e; };
	double GetError(){ return error; };
	bool nodeDistCheck(SRvec3& pos, double radius);
	bool InsideBoundingBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
	void dumpData();
	void checkMaterialStressMax(double s);
	void checkMaterialStrainMax(double e);
	int GetNumCorners();
	int GetLocalCornerIdFromGlobalFaceNodeId(int faceCornerId);
	double minCornerDistance(SRvec3& pos);
	double centroidDistance(SRvec3& pos);
	int getNodeOrMidNodeId(int n);
	SRnode* getNodeOrMidNode(int n);
	void approxCentroid(SRvec3& elcen);
	bool checkSlopeKink();
	void checkAllNodesHaveDisp();
	bool checkNearBktBdry();
	void scaleStress(double stress[], double& svm, double scaleClip);
	bool checkAnyNodeBreakout();
	bool checkForNodalHotSpot();
	void CalculateHStress(double r, double s, double t, double stress[]);
	double CalculateSvmHCentroid();

	int thread;
	double error;
	int userId;
	int id;
	SRelementType type;
	SRintVector nodeIds;
	SRvector <SRlocalEdge> localEdges; //num. of edges=localedges.num
	SRvector <SRlocalFace> localFaces; //num. of faces=localFaces.num
	SRmaterial* elMat;
	SRintVector globalFunctionNumbers;
	SRdoubleMatrix rawStrains;
	SRdoubleVector smoothedStrains[6];
	int sacrificial;
	int newPorder;
	double svmMax;
	int stiffLength;
	SRvector <int> nodeLcSConstraints;
	SRvector <int> faceLCSConstraints;
	int maxPorder;
	int filePos;
	double penaltyStiffness;
	SRdoubleVector stiffnessMatrix;
	SRintVector stiffDiag;

	double* xnode;
	double* ynode;
	double* znode;
	double* dNdr;
	double* dNds;
	double* dNdt;
	SRdoubleVector xnodev, ynodev, znodev;
	SRdoubleVector dNdrv, dNdsv, dNdtv;

	SRmat33 elJac;
	double size;
	int numNodesTotal;
	bool pChanged;
	double elVol;
	double approxVol;
	double maxStress; //this is max stress for any component, not svm
	double *basisVec;
	double *dbasisdr;
	double *dbasisds;
	double *dbasisdt;
	SRdoubleVector dispEl;
	bool constrained;
	bool saveForBreakout;
	bool flattenedHighStress;
	bool flattened;
	double flattenFraction;
	bool allNodesHaveDisp;
	double minJac;
};

#endif //!defined(SRELEMENT_INCLUDED)
