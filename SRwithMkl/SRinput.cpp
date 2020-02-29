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
// SRinput.cpp: implementation of the SRinput class.
//
//////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <search.h>
#include "SRmodel.h"
#include "SRmachDep.h"
#include "SRinput.h"

extern SRmodel model;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif


void SRinput::ReadModel()
{
	//read the model from input file in model.inputFile
	//also read any run settings specific to this model

	SCREENPRINT("reading input model\n");

	SRstring filename, line, basename, tail, tok;

	//assign output, results, and combined results names; delete existing
	//files in case future opens will be for "append":
	model.outdir.Right('\\', tail);
	basename = model.outdir;
	basename += "\\";
	basename += tail;

	filename = basename;
	filename += ".out";
	model.outputFile.filename = filename;
	model.outputFile.Delete();

	filename = basename;
	filename += ".sum";
	model.sumFile.filename = filename;
	model.sumFile.Delete();

	filename = basename;
	filename += ".frd";
	model.cgxFrdFile.filename = filename;

	filename = basename;
	filename += "_SR.f06";
	model.f06outFile.filename = filename;

	GetBinaryFileName("elstiffEven", filename);
	model.elementFileEven.filename = filename;
	if(SRfile::Existcheck(filename))
		model.elementFileEven.Delete();
	GetBinaryFileName("elstiffOdd", filename);
	model.elementFileOdd.filename = filename;
	if(SRfile::Existcheck(filename))
		model.elementFileOdd.Delete();

	LOGPRINT("\nReading Settings...");
	readSettings();

	if (!model.inputFile.Open(SRinputMode))
	{
		const char *tmp = filename.LastChar('\\', true);
		LOGPRINT(" file not found: %s", tmp);
		ERROREXIT;
		return;
	}

	LOGPRINT("\nReading mesh...");

	/* read model: */
	numUnSupportedNode = 0;
	nnode = nelem = ncon = nforce = nvolforce = nfacePressure = nfaceTraction = 0;
	nBreakoutCon = nFaceCon = 0;
	nmat = ncoord = numBreakoutCon = 0;
	numNodalBreakoutCons = 0;
	numFaceFromNodalLoad = 0;
	anyFaceHasMultipleLoads = false;

	model.materials.Allocate(MAXNMAT);
	model.Coords.Allocate(MAXNCOORD);
	//create default GCS coordinate system (needed for lcs constraints)
	SRcoord* coord = model.Coords.Add();
	coord->Create(0.0, 0.0, 0.0);
	coord->name.Copy("SR_DEFAULT_GCS");

	model.inputFile.GetLine(line);
	model.inputFile.ToTop();

	//count entities: nodes, elements, forces, local output nodes:
	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;

		if (line.isCommentOrBlank())
			continue;

		if (line == "nodes")
			CountEntities(nnode);
		else if (line == "elements")
			CountElements(nelem);
		else if (line == "forces")
			CountEntities(nforce);
		else if (line == "facePressures")
			CountEntities(nfacePressure);
		else if (line == "faceTractions")
			CountEntities(nfaceTraction);
		else if (line == "volumeForce")
			CountEntities(nvolforce);
		else if (line == "constraints")
			CountEntities(ncon);
		else if (line == "breakoutConstraints")
			CountEntities(nBreakoutCon);
		else if (line == "faceConstraints")
			CountEntities(nFaceCon);
		else if (line == "nodalbreakoutconstraints")
			CountEntities(numNodalBreakoutCons);
		else if (line.CompareUseLength("coord"))
		{
			//input coordinates and materials on 1st pass because
			//they will be referred to by name:
			InputCoordinates();
		}
		else if (line == "materials")
			InputMaterials();
	}

	if (nnode == 0 || nelem == 0)
	{
		LOGPRINT("Mesh definition incomplete. Missing nodes or elements");
		ERROREXIT;
	}

	if (!userPSettings && !model.saveBreakout && nelem > 40000)
	{
		//auto econsolve:
		model.maxPorder = model.maxPorderFinalAdapt = 5;
		model.adaptLoopMax = 2;
	}

	OUTPRINT("Number of elements: %d", nelem);
	OUTPRINT("Number of nodes: %d", nnode);
	OUTPRINT("\n");

	model.elements.Allocate(nelem);

	model.volumeForces.Allocate(nvolforce);
	//worst case of model with all bricks. Not wasteful because
	//these are just pointers
	int nedge = 12 * nelem;
	model.edges.Allocate(nedge);
	nnode += nedge;

	nforce += nfacePressure;
	nforce += nfaceTraction;
	model.forces.Allocate(nforce);

	model.nodes.Allocate(nnode);
	model.inputFile.ToTop();

	while (1)
	{
		if (!model.inputFile.GetLine(line))
		{
			OUTPRINT("mesh input file incomplete");
			ERROREXIT;
		}

		if (line.isCommentOrBlank())
			continue;

		if (line == "nodes")
		{
			InputNodes();
			break;
		}
	}

	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank())
			continue;

		if (line == "elements")
		{
			InputElements();
			break;
		}
	}

	model.constraints.Allocate(ncon + numUnSupportedNode);

	bool doMeshOnlyFrdPlot = false;
	if (doMeshOnlyFrdPlot)
		MeshOnlyFrdPlot();


	int nel = model.elements.GetNum();

	model.SetEdgesToPorder(model.minPorder);

	//worst case for number of faces:
	int nface;
	if (model.anybricks)
		nface = 6 * nel;
	else if (model.anywedges)
		nface = 5 * nel;
	else
		nface = 4 * nel;

	model.faces.Allocate(nface);

	nnode = model.GetNumNodes();
	numNodeFaces.Allocate(nnode);
	nodeFaces.Allocate(nnode, MAXNODEFACEOWNERS);

	//ttd!! rethink if full model unsup has regs
	//if (!model.saveBreakout || numUnSupportedNode > 0)
	if (!model.saveBreakout)
	{
		model.FillGlobalFaces();
		nodalToFaceUnsup();
		model.map.Setup();
	}

	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;

		curLine = line;

		if (line.isCommentOrBlank())
			continue;
		if (line == "forces")
			InputNodalForces();
		else if (line == "facePressures")
			InputFacePressures();
		else if (line == "faceTractions")
			InputFaceTractions();
		else if (line == "volumeForce")
			InputVolumeForces();
		else if (line == "thermal")
			InputThermal();
		else if (line == "constraints")
			InputNodalConstraints();
		else if (line == "breakoutConstraints")
			InputBreakoutConstraints(nBreakoutCon);
		else if (line == "nodalbreakoutconstraints")
			InputnodalBreakoutConstraints(numNodalBreakoutCons);
		else if (line == "faceConstraints")
			InputFaceConstraints(nFaceCon);
		else if (line == "bktElSacrList")
		{
			InputBreakoutSacrElems();
			PlotBktElSacr();
		}
	}

	model.inputFile.Close();


	dumpModel();

	numNodeEdges.Free();
	nodeEdges.Free();

	if (numUnSupportedNode != 0 && !model.saveBreakout)
	{
		OUTPRINT("Non solid elements encountered. See limitations in user's guide");
		LOGPRINT("Non solid elements encountered. See limitations in user's guide");
	}

	if (numUnSupportedNode != 0)
		applyNodalBreakoutCons();

	fillMultiFaceForceGroups();
	EnergySmoothNodalForces();
}

void SRinput::InputNodes()
{
	//input nodes

	//Input Spec:
		//userId x y z (followed by 1 if needs breakout else 0 or omitted)
	//return:
		//no. of nodes that need breakout constraints. These nodes are marked as sacrificial

	numUnSupportedNode = 0;
	nodeUidOffset = 0;
	SRstring line, tok;
	SRnode* node;
	int nnode = 0, uid;
	double x, y, z, xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = BIG;
	xmax = -BIG;
	ymin = BIG;
	ymax = -BIG;
	zmin = BIG;
	zmax= -BIG;
	while(1)
	{
		if(!model.inputFile.GetLine(line))
			break;
        if(line.isCommentOrBlank())
            continue;
		if(line == "end")
			break;
		node = model.addNode();
		line.TokRead(uid);
		if (nnode == 0)
			nodeUidOffset = uid;
		else if (nodeUidOffset != -1)
		{
			if (uid - nodeUidOffset != nnode)
				nodeUidOffset = -1;
		}
		line.TokRead(x);
		line.TokRead(y);
		line.TokRead(z);
		tok = line.Token();
		if (tok.CompareUseLength("unsup"))
		{
			//translator marked this node as connected to an unsupported element
			//in the original FEA mesh
			node->sacrificialType = unsupported;
			model.anyUnsupNode = true;
			numUnSupportedNode++;
		}
		else if (tok.CompareUseLength("onBsurf"))
		{
			//translator marked this node as connected to an unsupported element
			//in the original FEA mesh
			node->sacrificialType = onBsurf;
			model.anyUnsupNode = true;
			model.anyBsurf = true;
			numUnSupportedNode++;
		}
		else if (tok.CompareUseLength("shellOrBeamNode"))
		{
			//translator marked this node as connected to a shell or beam element
			//in the original FEA mesh
			node->sacrificialType = shellOrBeamNode;
			model.anyUnsupNode = true;
			model.numShellOrBeamNodes++;
		}

		if (x < xmin)
			xmin = x;
		if(x > xmax)
			xmax = x;
		if(y < ymin)
			ymin = y;
		if(y > ymax)
			ymax = y;
		if(z < zmin)
			zmin = z;
		if(z > zmax)
			zmax = z;
		node->Create(uid, x, y, z);
		node->id = nnode;
		if (uid > model.maxNodeUid)
			model.maxNodeUid = uid;
		nnode++;
	}
	double dx, dy, dz;
	dx = xmax - xmin;
	dy = ymax - ymin;
	dz = zmax - zmin;
	model.size = sqrt(dx*dx + dy*dy + dz*dz);

#if 0
	OUTPRINT(" model bounding box:");
	OUTPRINT(" xmin: %lg xmax: %lg", xmin, xmax);
	OUTPRINT(" ymin: %lg ymax: %lg", ymin, ymax);
	OUTPRINT(" zmin: %lg zmax: %lg", zmin, zmax);
#endif

	FillAndSortNodeUids();

	numNodeEdges.Allocate(nnode);
	nodeEdges.Allocate(nnode, 20);
}

void SRinput::InputElements()
{
	//Input Spec:
		//userId Material-name nodeuserId,i=1 to # of nodes
			//# of nodes=
			//10 for a tet, 15 for a wedge, 20 for a brick
			//for linear mesh # of nodes=
			//4 for a tet, 6 for a wedge, 8 for a brick

	SRstring line, tok;

	SRelement* elem;
	int nelem = 0;
	int uid, nodeuid, mid, nodev[20], nt, id;
	SRnode* node;
	elemUidOffset = 0;
	while (1)
	{
		if(!model.inputFile.GetLine(line))
			break;
        if(line.isCommentOrBlank())
            continue;
		if(line == "end")
			break;
		elem = model.GetElement(nelem);
		elem->id = nelem;
		line.TokRead(uid);
		if (nelem == 0)
			elemUidOffset = uid;
		else if (elemUidOffset != -1)
		{
			if (uid - elemUidOffset != nelem)
				elemUidOffset = -1;
		}

		if (uid > model.maxElemUid)
			model.maxElemUid = uid;

		tok = line.Token();
		mid = GetMaterialId(tok);
		SRmaterial* mat = model.materials.GetPointer(mid);
		mat->numElements++;
		nt = 0;
		while(1)
		{
			if(!line.TokRead(nodeuid))
				break;
			id = NodeFind(nodeuid);
			nodev[nt] = id;
			node = model.GetNode(id);
			if (node->firstElementOwner == -1)
				node->firstElementOwner = elem->id;
			nt++;
		}
		if (nt == 4 || nt == 6 || nt == 8)
			ERROREXIT; //linear mesh not supported

		if (nt == 10)
		{
			elem->Create(tet, uid, nt, nodev, mat);
			model.anytets = true;
		}
		else if (nt == 15)
		{
			elem->Create(wedge, uid, nt, nodev, mat);
			model.anywedges = true;
		}
		else if (nt == 20)
		{
			elem->Create(brick, uid, nt, nodev, mat);
			model.anybricks = true;
		}
		else
		{
			OUTPRINT("incorrect element definition in mesh file");
			ERROREXIT;
		}

		if (!model.saveBreakout || (numUnSupportedNode != 0))
			model.CreateElemEdges(elem, nt, nodev);
		nelem++;
	}

	bool anyorphan = false;
	for (int i = 0; i < model.nodes.GetNum(); i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
		{
			anyorphan = true;
			break;
		}
	}
	if (anyorphan)
		OUTPRINT("Ignoring Nodes that are not Referenced by elements");

	if (elemUidOffset == -1)
	{
		//fill elem uid vector and sort in ascending order of uid for faster
		//elem-finding:
		elemUids.Allocate(nelem);
		SRuidData *euid;
		for (int i = 0; i < nelem; i++)
		{
			elem = model.GetElement(i);
			euid = elemUids.GetPointer(i);
			euid->id = i;
			euid->uid = elem->userId;
		}
		SortElems();
	}

}

void SRinput::InputMaterials()
{
	//input materials

	//Input Spec:
		//name type
			//type = isotropic, orthotropic or general
			//if iso: rho, alpha then E,nu on second line
			//if ortho: rho, alphax, alphay, alphaz, then c11,c12,c13,c22,c23,c33,c44,c55,c66
			//if general: rho, alphax, alphay, alphaz, then full cij matrix, 36 constants,6 per line
				//(must be symmetric)

	SRstring line;
	SRmaterial* mat;
	SRstring type;
	double E,nu;

	SRfile* f;
	f = &model.inputFile;

	while(1)
	{
		double alpha;
		if(!f->GetLine(line))
			break;
        if(line.isCommentOrBlank())
            continue;
		if(line == "end")
			break;
		if(model.GetNumMaterials() == MAXNMAT)
		{
			LOGPRINT("Maximum Number of Materials allowed is %d\n", MAXNMAT);
			ERROREXIT;
		}
		mat = model.materials.Add();
		mat->id = model.materials.GetNum() - 1;
		mat->name = line.Token();
		type = line.Token();
		f->GetLine(line);
		line.TokRead(mat->rho);
		if(type == "iso")
		{ 
			mat->type = iso;
			line.TokRead(alpha);
			mat->alphax = alpha;
			mat->alphay = alpha;
			mat->alphaz = alpha;
			line.TokRead(mat->tref, CHECKFORTRAILINGCOMMENT);
			if (line.TokRead(mat->allowableStress, CHECKFORTRAILINGCOMMENT))
			{
				if (mat->allowableStress > SMALL)
					mat->allowableAssigned = true;
			}
			f->GetLine(line);
			line.TokRead(E);
			line.TokRead(nu);
			mat->IsoCreate(E,nu);
		}
		else if(type == "ortho")
		{
			mat->type = ortho;
			line.TokRead(mat->alphax);
			line.TokRead(mat->alphay);
			line.TokRead(mat->alphaz);
			line.TokRead(mat->tref);
			if (line.TokRead(mat->allowableStress))
			{
				if (mat->allowableStress > SMALL)
					mat->allowableAssigned = true;
			}
			f->GetLine(line);
			line.TokRead(mat->orthoCij.c11);
			line.TokRead(mat->orthoCij.c12);
			line.TokRead(mat->orthoCij.c13);
			line.TokRead(mat->orthoCij.c22);
			line.TokRead(mat->orthoCij.c23);
			line.TokRead(mat->orthoCij.c33);
			line.TokRead(mat->orthoCij.c44);
			line.TokRead(mat->orthoCij.c55);
			line.TokRead(mat->orthoCij.c66);
		}
		else if(type == "gen")
		{
			mat->type = genAniso;
			line.TokRead(mat->alphax);
			line.TokRead(mat->alphay);
			line.TokRead(mat->alphaz);
			line.TokRead(mat->tref);
			if (line.TokRead(mat->allowableStress))
			{
				if (mat->allowableStress > SMALL)
					mat->allowableAssigned = true;
			}
			for(int i = 0; i < 6; i++)
			{
				f->GetLine(line);
				for(int j = 0; j < 6; j++)
					line.TokRead(mat->genAnisoCij.c[i][j]);
			}
			if (!mat->genAnisoCij.symCheck())
			{
				OUTPRINT("improper definition of general anisotropic material %s", mat->name.getStr());
				OUTPRINT("material is not symmetric");
				REPPRINT("improper definition of general anisotropic material %s", mat->name.getStr());
				REPPRINT("material is not symmetric");
				ERROREXIT;
			}

		}
		else
		{
			OUTPRINT("improper type in definition of material");
			ERROREXIT;
		}
	}

	equivMatTest();
}

void SRinput::InputCoordinates()
{
	//input local coordinate systems
	//Input Spec: name type "NotGcsAligned"
	//type = cartesian,spherical,cylindrical
	//x0,y0,z0 (origin)
	//if NotGcsAligned:
	//p1, p3 are points along local e1 and e3 axes

	SRstring line, type, tok;
	SRcoord* coord;
	int ncoord = 0;
	double x0, y0, z0;
	bool gcsaligned;

	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank())
			continue;
		if (line == "end")
			break;
		if (ncoord > MAXNCOORD)
		{
			OUTPRINT("Maximum Number of Coordinate Systems allowed is %d\n", MAXNCOORD);
			LOGPRINT("Maximum Number of Coordinate Systems allowed is %d\n", MAXNCOORD);
			ERROREXIT;
		}
		coord = model.Coords.Add();
		coord->name = line.Token();
		type = line.Token();
		tok = line.Token();
		if (tok == "NotGcsAligned")
			gcsaligned = false;
		else
			gcsaligned = true;
		model.inputFile.GetLine(line);
		line.TokRead(x0);
		line.TokRead(y0);
		line.TokRead(z0);
		if (type == "cartesian")
			coord->type = cartesian;
		else if (type == "spherical")
			coord->type = spherical;
		else if (type == "cylindrical")
			coord->type = cylindrical;
		else
			ERROREXIT;

		if (gcsaligned)
			coord->Create(x0, y0, z0);
		else
		{
			SRvec3 p13, p3;
			model.inputFile.GetLine(line);
			line.TokRead(p13.d[0]);
			line.TokRead(p13.d[1]);
			line.TokRead(p13.d[2]);
			line.TokRead(p3.d[0]);
			line.TokRead(p3.d[1]);
			line.TokRead(p3.d[2]);
			coord->Create(x0, y0, z0, p13, p3);
		}
		ncoord++;
	}
}

void SRinput::InputNodalForces()
{
	//input nodal forces

	//input spec (all on one line):
	//node-Id
	//pressure "pressure" if pressure load (only valid for nodal loads that will become face loads) OR
	//or "coord" coord name
	//or "gcs"
	//if "pressure", 
	// pressure value
	// else
	//force vals for 3 dofs; 0 if not forced
	//notes:
	// edge and face forces will be automatically created if all nodes are loaded compatibly
	// and if doNodalToFace

	SRstring tok, line;

	SRforce* force;
	int nodeuid;
	SRforce forceTmp;

	int nnode = model.GetNumNodes();
	currentNodeForceStore.Allocate(nnode);

	currentForces.Allocate(nnode);
	if (nodeForceStore.GetNum() == 0)
		nodeForceStore.Allocate(nnode);

	bool rotateNodalLcs = true;

	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank())
			continue;
		if (line == "end")
			break;
		line.TokRead(nodeuid);
		int nId = NodeFind(nodeuid);
		if (nId == -1)
			continue;//node doesn't exist in model, probably is on an unsupported entity e.g. shell or beam

		//read in the force, then check for duplicate forces on this node:
		forceTmp.Clear();
		forceTmp.entityId = nId;
		tok = line.Token();
		if (tok == "pressure")
		{
			forceTmp.pressure = true;
			forceTmp.forceVals.Allocate(1, 1);
			double p;
			line.TokRead(p);
			forceTmp.forceVals.PlusAssign(0, 0, p);
		}
		else if (tok == "coord")
		{
			tok = line.Token();
			int cId = GetCoordId(tok);
			forceTmp.coordId = cId;
		}
		else if (!tok.Compare("gcs"))
		{
			OUTPRINT("incorrect force type on node %d. must be 'pressure', 'coord Name' or 'gcs'", nodeuid);
			ERROREXIT;
		}

		if (!forceTmp.pressure)
		{
			//for non pressure forces there are 3 force values. if the user omits the 2nd or 3rd they are set to 0
			forceTmp.forceVals.Allocate(1, 3);
			double ft;
			for (int dof = 0; dof < 3; dof++)
			{
				if (!line.TokRead(ft))
					ft = 0.0;
				forceTmp.forceVals.PlusAssign(0, dof, ft);
			}
			if (forceTmp.coordId != -1 && rotateNodalLcs)
			{
				SRvec3 fl, fg;
				SRmat33 R;
				for (int dof = 0; dof < 3; dof++)
					fl.d[dof] = forceTmp.forceVals.Get(0, dof);
				SRcoord* coord = model.GetCoord(forceTmp.coordId);
				SRnode* node = model.GetNode(nId);
				coord->GetRotationMatrix(false, node->pos, R);

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
						fg.d[i] += (R.rows[i].d[j] * fl.d[j]);
				}
				for (int dof = 0; dof < 3; dof++)
					forceTmp.forceVals.Put(0, dof, fg.d[dof]);
				forceTmp.coordId = -1;
			}
		}

		int nnodalForce = currentNodeForceStore.d[nId].forces.GetNum();

		//check for duplicate force, same node, could happen e.g. at a corner
		//if two edges sharing a face loaded.
		bool dupe = false;
		for (int iprev = 0; iprev < nnodalForce; iprev++)
		{
			int prevId = currentNodeForceStore.d[nId].GetForceId(iprev);
			force = currentForces.GetPointer(prevId);
			if (force->coordId == forceTmp.coordId && force->pressure == forceTmp.pressure)
			{
				//duplicate and compatible force
				//add the force vectors
				dupe = true;
				force->AddNodalForce(forceTmp);
				break;
			}
		}
		if (!dupe)
		{
			//new force
			int forceid = currentForces.GetNum();
			currentNodeForceStore.d[nId].forces.PushBack(forceid);
			force = currentForces.Add();
			force->Copy(forceTmp);
		}
	}

	//AND currentForces and model.forces and update nodeForceStore:
	for (int i = 0; i < nnode; i++)
	{
		int nnodalForce = nodeForceStore.d[i].forces.GetNum();
		int curnumnodalForce = currentNodeForceStore.d[i].forces.GetNum();
		for (int curforcenum = 0; curforcenum < curnumnodalForce; curforcenum++)
		{
			int curforceid = currentNodeForceStore.d[i].GetForceId(curforcenum);
			SRforce* curforce = currentForces.d[curforceid];
			SRforce* force;
			bool dupe = false;
			for (int forcenum = 0; forcenum < nnodalForce; forcenum++)
			{
				force = nodeForceStore.d[i].GetForce(forcenum);
				if (force->coordId == curforce->coordId && force->pressure == curforce->pressure)
				{
					force->AddNodalForce(*curforce, true); //summing sets is true
					dupe = true;
				}
			}
			if (!dupe)
			{
				//new force:
				int forceid = model.forces.GetNum();
				force = model.forces.Add();
				nodeForceStore.d[i].forces.PushBack(forceid);
				force->Copy(*curforce);
			}
		}
	}

	//catch case that only midside nodes were loaded. This can occur for constant force on a face,
	//equiv nodal load at corners is 0 so femap doesn't put it out
	int numNewCornerForces = 0;
	for (int i = 0; i < nnode; i++)
	{
		int nnodalForce = nodeForceStore.d[i].forces.GetNum();
		SRnode* node = model.GetNode(i);
		if (nnodalForce == 0)
			continue;
		int eid = node->midSideEdgeOwner;
		if (eid == -1)
			continue;
		SRedge* edge = model.GetEdge(eid);
		for (int n = 0; n < 2; n++)
		{
			int nid = edge->GetNodeId(n);
			if (nodeForceStore.d[nid].forces.GetNum() == 0)
				numNewCornerForces++;
		}
	}
	if (numNewCornerForces != 0)
	{
		int nforce = model.forces.GetNumAllocated();
		nforce += numNewCornerForces;
		nforce += nfacePressure;
		nforce += nfaceTraction;
		model.forces.Allocate(nforce);
		for (int i = 0; i < nnode; i++)
		{
			int nnodalForce = nodeForceStore.d[i].forces.GetNum();
			SRnode* node = model.GetNode(i);
			if (nnodalForce == 0)
				continue;
			int eid = node->midSideEdgeOwner;
			if (eid == -1)
				continue;
			SRedge* edge = model.GetEdge(eid);
			for (int n = 0; n < 2; n++)
			{
				int nid = edge->GetNodeId(n);
				if (nodeForceStore.d[nid].forces.GetNum() == 0)
				{
					for (int f = 0; f < nnodalForce; f++)
					{
						force = nodeForceStore.d[i].GetForce(f);
						int forceid = model.forces.GetNum();
						SRforce* cornerforce = model.forces.Add();
						cornerforce->Copy(*force);
						cornerforce->entityId = nid;
						cornerforce->forceVals.Zero();
						nodeForceStore.d[nid].forces.PushBack(forceid);
					}
				}
			}
		}
	}
	doNodalToFaceForces();
}

void SRinput::InputFacePressures()
{
	//input pressures on faces

	//input spec (all on one line):
	//eluid, n1, n2, n3, n4, p1, p2, p3, p4
	//eluid = uid of element that owns the face
	//n1, n2, n3, n4 = nodes at corner of face, n4 = 1 for tri;
	//p1, p2, p3, p4 = pressures at corner of face, p4 omitted for tri face

	SRstring tok, line;

	SRforce* force;
	int gno[4];
	int i;
	int nfacep = 0;
	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank())
			continue;
		if (line == "end")
			break;
		int forceId = model.forces.GetNum();
		force = model.forces.Add();
		force->type = faceForce;
		force->pressure = true;
		int eluid;
		line.TokRead(eluid);
		int elId = elemFind(eluid);
		SRelement* elem = model.GetElement(elId);
		int nv[4], nuid[4];
		for (i = 0; i < 4; i++)
		{
			line.TokRead(nuid[i]);
			if (nuid[i] != -1)
				nv[i] = NodeFind(nuid[i]);
			else
				nv[i] = -1;
		}
		int fId = -1;
		if (!model.saveBreakout)
		{
			fId = elemFaceFind(elem, nv, gno);
			if (fId == -1)
			{
				OUTPRINT(" improperly defined pressure load on face. element user id: %d", eluid);
				OUTPRINT(" face node user ids: %d %d %d %d", nuid[0], nuid[1], nuid[2], nuid[3]);
				OUTPRINT(" face not found in the mesh");
				ERROREXIT;
			}
			force->entityId = fId;
			SRface *face = model.GetFace(fId);
			face->forceIds.PushBack(forceId);
		}
		else
		{
			for (int i = 0; i < 4; i++)
			{
				force->nv[i] = nv[i];
				gno[i] = i;
			}
		}
		force->elemId = elId;
		double pv[8];
		int nn = 4;
		if (nv[3] == -1)
			nn = 3;
		for (i = 0; i < nn; i++)
			line.TokRead(pv[i]);

		force->forceVals.Allocate(2 * nn, 1);
		for (i = 0; i < nn; i++)
		{
			int idg = gno[i];
			force->forceVals.Put(idg, 0, pv[i]);
		}
		//interpolate on face to get pv of midnodes:
		if (nn == 3)
		{
			//tri face
			pv[3] = 0.5*(force->forceVals.Get(0, 0) + force->forceVals.Get(1, 0));
			pv[4] = 0.5*(force->forceVals.Get(1, 0) + force->forceVals.Get(2, 0));
			pv[5] = 0.5*(force->forceVals.Get(0, 0) + force->forceVals.Get(2, 0));
		}
		else
		{
			//quad face
			pv[4] = 0.5*(force->forceVals.Get(0, 0) + force->forceVals.Get(1, 0));
			pv[5] = 0.5*(force->forceVals.Get(1, 0) + force->forceVals.Get(2, 0));
			pv[6] = 0.5*(force->forceVals.Get(3, 0) + force->forceVals.Get(2, 0));
			pv[7] = 0.5*(force->forceVals.Get(0, 0) + force->forceVals.Get(3, 0));
		}

		for (int l = 0; l < nn; l++)
		{
			int m = l + nn;
			//linearly interpolate from corners to get midside values:
			force->forceVals.Put(m, 0, pv[m]);
		}
	}
	nfacep++;
}


void SRinput::InputFaceTractions()
{
	//input tractions on faces

	//input spec (all on one line):
	//eluid, n1, n2, n3, n4
	//eluid = uid of element that owns the face
	//n1, n2, n3, n4 = uids of nodes at corner of face, n4 = -1 for tri
	//then 3 continuation lines of t1, t2, t3, t4
	// (tractions at corner of face, t4 omitted for tri face); 1 for each dof

	SRstring tok, line;

	SRforce* force = NULL;
	int gno[4];
	int i;

	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank())
			continue;
		if (line == "end")
			break;
		int eluid;
		line.TokRead(eluid);
		int elId = elemFind(eluid);
		SRelement* elem = model.GetElement(elId);
		int nv[4], nuid[4];
		for (i = 0; i < 4; i++)
		{
			line.TokRead(nuid[i]);
			nv[i] = NodeFind(nuid[i]);
		}

		bool found = false;
		int fId = -1;
		if (!model.saveBreakout)
		{
			fId = elemFaceFind(elem, nv, gno);
			if (fId == -1)
			{
				OUTPRINT(" improperly defined traction load on face. element user id: %d", eluid);
				OUTPRINT(" face node user ids: %d %d %d %d", nuid[0], nuid[1], nuid[2], nuid[3]);
				OUTPRINT(" face not found in the mesh");
				ERROREXIT;
			}

			//see if a force with fId exists:
			for (i = 0; i < model.GetNumForces(); i++)
			{
				force = model.GetForce(i);
				if (force->GetEntityId() == fId)
				{
					if (force->isPressure())
					{
						OUTPRINT(" duplicate pressure and traction load on face. element user id: %d", eluid);
						OUTPRINT(" face node user ids: %d %d %d %d", nuid[0], nuid[1], nuid[2], nuid[3]);
						ERROREXIT;
					}

					found = true;
					break;
				}
			}
		}

		if (!found)
		{
			int forceId = model.forces.GetNum();
			force = model.forces.Add();
			force->entityId = fId;
			force->type = faceForce;
			if (!model.saveBreakout)
			{
				SRface* face = model.GetFace(fId);
				face->forceIds.PushBack(forceId);
			}
			else
			{
				for (int i = 0; i < 4; i++)
				{
					force->nv[i] = nv[i];
					gno[i] = i;
				}
			}
		}
		force->elemId = elId;

		force->type = faceForce;
		for (int dof = 0; dof < 3; dof++)
		{

			if (!model.inputFile.GetLine(line) || !line.CompareUseLength("#"))
			{
				if (dof < 2)
				{
					OUTPRINT(" not enough degrees of freedom specified for traction load on face. element user id: %d", eluid);
					OUTPRINT(" face node user ids: %d %d %d %d", nuid[0], nuid[1], nuid[2], nuid[3]);
					ERROREXIT;
				}
			}
			tok = line.Token(); //skip the # character
			double tv[8];
			int nn = 4;
			if (nv[3] == -1)
				nn = 3;
			line.TokRead(tv[0]);
			for (i = 1; i < nn; i++)
				tv[i] = tv[0];
			for (i = 1; i < nn; i++)
			{
				double tt;
				if (!line.TokRead(tt))
					break;
				tv[i] = tt;
			}
			if (dof == 0)
				force->forceVals.Allocate(2 * nn, 3);
			for (i = 0; i < nn; i++)
			{
				int idg = gno[i];
				force->forceVals.Put(idg, dof, tv[i]);
			}
			//interpolate on face to get tv of midnodes:
			if (nn == 3)
			{
				//tri face
				tv[3] = 0.5*(force->forceVals.Get(0, dof) + force->forceVals.Get(1, dof));
				tv[4] = 0.5*(force->forceVals.Get(1, dof) + force->forceVals.Get(2, dof));
				tv[5] = 0.5*(force->forceVals.Get(0, dof) + force->forceVals.Get(2, dof));
			}
			else
			{
				//quad face
				tv[4] = 0.5*(force->forceVals.Get(0, dof) + force->forceVals.Get(1, dof));
				tv[5] = 0.5*(force->forceVals.Get(1, dof) + force->forceVals.Get(2, dof));
				tv[6] = 0.5*(force->forceVals.Get(3, dof) + force->forceVals.Get(2, dof));
				tv[7] = 0.5*(force->forceVals.Get(0, dof) + force->forceVals.Get(3, dof));
			}

			for (int l = 0; l < nn; l++)
			{
				int m = l + nn;
				//linearly interpolate from corners to get midside values:
				force->forceVals.Put(m, dof, tv[m]);
			}
		} // for (int dof = 0; dof < 3; dof++)

	} // while (1)
}

void SRinput::InputVolumeForces()
{
	//input volume force (gravity or centrifugal)
	//input spec:
		//type (gravity,centrifugal)
		//g1,g2,g3 for gravity or
		//omega, axis, origin for centrifugal

	SRstring tok, line;

	SRvolumeForce* force;
	double omega;
	bool nextLineWasRead = false;
	while(1)
	{
		if (!nextLineWasRead)
		{
			if (!model.inputFile.GetLine(line))
				break;
		}
        if(line.isCommentOrBlank())
            continue;
		if (line == "end")
			break;
		force = model.volumeForces.Add();
		if (line == "gravity")
		{
			force->type = gravity;
			tok = line.Token();
			line.TokRead(force->g1);
			line.TokRead(force->g2);
			line.TokRead(force->g3);
		}
		else if (line == "centrifugal")
		{
			force->type = centrifugal;
			tok = line.Token();
			line.TokRead(omega);
			force->omega2 = omega*omega;
			line.TokRead(force->axis.d[0]);
			line.TokRead(force->axis.d[1]);
			line.TokRead(force->axis.d[2]);
			line.TokRead(force->origin.d[0]);
			line.TokRead(force->origin.d[1]);
			line.TokRead(force->origin.d[2]);
			double alpha = 0.0;
			line.TokRead(alpha);
			force->alpha = alpha;
		}
		else
		{
			OUTPRINT("improper volume force type. should be 'gravity' or 'centrifugal'");
			ERROREXIT;
		}
		if (!model.inputFile.GetLine(line))
			break;
		//optional continuation lines of element list
		if (line.CompareUseLength("#"))
		{
			nextLineWasRead = false;
			tok = line.Token(); //skip the # character
			tok = line.Token(); //skip the elements keyword
			int numel;
			line.TokRead(numel);
			force->elList.Allocate(numel);
			int nread = 0;
			int eluid, eid;
			while (nread < numel)
			{
				if (!model.inputFile.GetLine(line))
					break;
				tok = line.Token(); // skip the # character
				while (1)
				{
					if (!line.TokRead(eluid))
						break;
					eid = elemFind(eluid);
					force->elList.d[nread] = eid;
					nread++;
					if (nread == numel)
						break;
				}
			}
		}
		else
			nextLineWasRead = true;
	}
}

void SRinput::InputThermal()
{
	//input thermal loads
	//input spec:
		//"constant" temp
		//else "variable", then multiple lines of
		//nodeid, temperature, 1 for each node in model, one per line,
		//then "end"

	SRstring tok, line;

	model.thermalForce = ALLOCATEMEMORY SRthermalForce;
	SRthermalForce* therm = model.thermalForce;
	double temp;
	model.inputFile.GetLine(line);
	int uid;
	if (line == "constant")
	{
		tok = line.Token();
		therm->constantTemp = true;
		line.TokRead(therm->temp);
	}
	else if (line == "variable")
	{
		therm->constantTemp = false;
		int n = model.nodes.GetNum();
		therm->nodalTemp.Allocate(n);

		while(1)
		{
			model.inputFile.GetLine(line);
			if (line == "end")
				break;
			line.TokRead(uid);
			int i = NodeFind(uid);
			line.TokRead(temp);
			if (i != -1)
				therm->nodalTemp.Put(i, temp);
		}
	}
	else
		ERROREXIT;
}

void SRinput::InputNodalConstraints()
{
	//input nodal constraints

	//input spec: the following are all on one line:
	//node user id
	//enforced disp vals for each constrained dof, 0 if not enforced, "-" if not constrained
	//"coord" coord-name
	//notes:
		//edge and face constraints will be automatically created
		//if all nodes of the edge or face are constrained compatibly

	SRstring tok, line;
	SRconstraint* constraint;
	int dof;
	int nodeuid;

	SRconstraint conTmp;

	int nnode = model.GetNumNodes();

	currentNodeConStore.Allocate(nnode);

	currentConstraints.Allocate(nnode);

	if (nodeConStore.GetNum() == 0)
		nodeConStore.Allocate(nnode);


	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (line == "end")
			break;
		line.TokRead(nodeuid);
		int nId = NodeFind(nodeuid);
		if (nId == -1)
			continue;//node doesn't exist in model, probably is on an unsupported entity e.g. shell or beam

		SRnode* node = model.GetNode(nId);
		if (node->checkUnsup())
		{
			//any applied constraints on this node are overridden by the breakout constraint,
			//handled below
			continue;
		}
		//read in the constraint, then check for duplicate constraints on this node:
		conTmp.Clear();
		conTmp.type = nodalCon;
		conTmp.entityId = nId;
		double enfd;
		for (dof = 0; dof < 3; dof++)
		{
			conTmp.constrainedDof[dof] = false;
			tok = line.Token();
			if (tok.getLength() == 1 && tok == "-")
				continue;
			conTmp.constrainedDof[dof] = true;
			enfd = tok.RealRead();
			if (fabs(enfd) > TINY)
				model.anyEnforcedDisplacement = true;
			if (conTmp.enforcedDisplacementData.isEmpty())
				conTmp.enforcedDisplacementData.Allocate(1, 3);
			conTmp.enforcedDisplacementData.Put(0, dof, enfd);
		}
		tok = line.Token();
		if (tok == "coord")
		{
			tok = line.Token();
			conTmp.coordId = GetCoordId(tok);
		}

		int nnodalCon = currentNodeConStore.d[nId].cons.GetNum();

		//check for duplicate constraint, same node, could happen e.g. at a corner
		//if two edges sharing a face constrained.
		bool dupe = false;
		for (int iprev = 0; iprev < nnodalCon; iprev++)
		{
			int prevId = currentNodeConStore.d[nId].GetConId(iprev);
			constraint = currentConstraints.GetPointer(prevId);
			if (constraint->coordId == conTmp.coordId)
			{
				//duplicate and compatible constraint
				//add the constrainedDof flags and the enforced displacement vectors
				dupe = true;
				constraint->AddNodalConstraint(conTmp);
				break;
			}
		}
		if (!dupe)
		{
			//new constraint:
			int conid = currentConstraints.GetNum();
			currentNodeConStore.d[nId].cons.PushBack(conid);
			constraint = currentConstraints.Add();
			constraint->Copy(conTmp);
		}
	}

	//AND currentConstraints and model.constraints and update nodeConStore:
	for (int i = 0; i < nnode; i++)
	{
		int nnodalCon = nodeConStore.d[i].cons.GetNum();
		int curnumnodalCon = currentNodeConStore.d[i].cons.GetNum();
		for (int curconnum = 0; curconnum < curnumnodalCon; curconnum++)
		{
			int curconid = currentNodeConStore.d[i].GetConId(curconnum);
			SRconstraint* curcon = currentConstraints.GetPointer(curconid);
			SRconstraint* con;
			bool dupe = false;
			for (int connum = 0; connum < nnodalCon; connum++)
			{
				con = nodeConStore.d[i].GetCon(connum);
				if (con->coordId == curcon->coordId)
				{
					con->AddNodalConstraint(*curcon, true); //summing sets is true
					dupe = true;
				}
			}
			if (!dupe)
			{
				//new constraint:
				int conid = model.constraints.GetNum();
				con = model.constraints.Add();
				con->id = conid;
				nodeConStore.d[i].cons.PushBack(conid);
				con->Copy(*curcon);
				SRnode* node = model.GetNode(con->entityId);
				node->constraintId = con->id;
			}
		}
	}

	currentNodeConStore.Free();
	currentConstraints.Free();

	doNodalToFaceConstraints();
}



void SRinput::InputFaceConstraints(int nfc)
{
	//input face constraints
	//input spec:
	//1st line: elem uid of the face owner, then nodeuids of the face corners (4th nodeuid =-1 for triangular face), then 3 constrained dof flags (0 or 1)
	//then optional "coord" coord-name. constraint is gcs is coord flag omitted
	//Then for each constrained dof, continuation lines of enforced disp vals for each corner (or just one value if constant)
	//note:
	//the continuation lines must start with "#"

	SRstring tok, line;

	if (nfc == 0)
		return;

	int ncon = model.constraints.GetNum();
	ncon += nfc;
	model.constraints.Allocate(ncon);

	int gno[4];
	int i;
	bool condof[3];

	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (line == "end")
			break;
		int euid;
		line.TokRead(euid);
		int eid = elemFind(euid);
		if (eid == -1)
			ERROREXIT;
		SRelement *elem = model.GetElement(eid);
		int nv[4];
		nv[3] = -1;
		for (i = 0; i < 4; i++)
		{
			if (!line.TokRead(nv[i]))
				break;
			nv[i] = NodeFind(nv[i]);
		}
		int nnodes = nv[3] == -1 ? 3 : 4;
		if (nv[3] == -1)
			nnodes = 3;
		for (i = 0; i < 3; i++)
		{
			condof[i] = false;
			int ct;
			if (!line.TokRead(ct))
				break;
			if (ct == 1)
				condof[i] = true;
		}

		int coordId = -1;
		tok = line.Token();
		if (tok.CompareUseLength("coord"))
		{
			tok = line.Token();
			coordId = GetCoordId(tok);
		}

		int fId = elemFaceFind(elem, nv, gno);
		if (fId == -1)
			ERROREXIT;
		int conId = model.constraints.GetNum();
		SRconstraint* constraint = model.constraints.Add();
		constraint->type = faceCon;
		constraint->entityId = fId;
		SRface* face = model.GetFace(fId);
		constraint->coordId = coordId;
		face->constraintId = conId;
		int nnodesTotal = face->GetNumNodesTotal();
		double enfd[8];
		double enfdm[8][3];
		bool allEndZero = true;
		bool enfdZero[3];
		for (int dof = 0; dof < 3; dof++)
		{
			enfdZero[dof] = true;
			if (!condof[dof])
				continue;
			constraint->SetConstrainedDof(dof);
			if (!model.inputFile.GetLine(line))
				ERROREXIT;
			tok = line.Token();//skip continuation character ("#")
			int nread = 0;
			for (i = 0; i < nnodesTotal; i++)
			{
				if (!line.TokRead(enfd[i]))
					break;
				if (fabs(enfd[i]) > TINY)
				{
					enfdZero[dof] = false;
					allEndZero = false;
				}
				nread++;
			}
			if (nread == 1)
			{
				for (i = 1; i < nnodes; i++)
					enfd[i] = enfd[0];

			}
			else if (nread > nnodesTotal)
				ERROREXIT;
			if (enfdZero[dof])
				continue;

			if (nread < nnodesTotal)
			{
				//interpolate on face to get enfd of midnodes:
				if (nnodes == 3)
				{
					//tri face
					enfd[3] = 0.5*(enfdm[0][dof] + enfdm[1][dof]);
					enfd[4] = 0.5*(enfdm[1][dof] + enfdm[2][dof]);
					enfd[5] = 0.5*(enfdm[0][dof] + enfdm[2][dof]);
				}
				else
				{
					//quad face
					enfd[4] = 0.5*(enfdm[0][dof] + enfdm[1][dof]);
					enfd[5] = 0.5*(enfdm[1][dof] + enfdm[2][dof]);
					enfd[6] = 0.5*(enfdm[3][dof] + enfdm[2][dof]);
					enfd[7] = 0.5*(enfdm[0][dof] + enfdm[3][dof]);
				}
			}
			for (i = 0; i < nnodesTotal; i++)
				enfdm[i][dof] = enfd[i];
		}
		if (!allEndZero)
		{
			model.anyEnforcedDisplacement = true;
			constraint->enforcedDisplacementData.Allocate(nnodesTotal, 3);
			for (int dof = 0; dof < 3; dof++)
			{
				if (!enfdZero[dof])
				{
					for (i = 0; i < nnodesTotal; i++)
						constraint->enforcedDisplacementData.Put(i, dof, enfdm[i][dof]);
				}
			}
		}
	}
}

void SRinput::InputBreakoutConstraints(int nbr)
{
	//input breakout constraints
	//input spec:
		//1st line: elem uid of the face owner, then local face number
		//Then n continuation lines of enforced disp vals for 3 constrained dofs where n is number of nodes and midnodes of face
	//note:
		//in input file. all breakout constraints are on faces.
		//the continuation lines must start with "#"

	if (nbr == 0)
		return;

	int ncon = model.constraints.GetNum();
	ncon += nbr;
	model.constraints.Allocate(ncon);

	SRstring tok, line;

	int gno[4];
	int i;

	while (1)
	{
		if (!model.inputFile.GetLine(line))
			break;
		if (line == "end")
			break;
		int euid;
		SRstring linesav = line;
		line.TokRead(euid);
		int eid = elemFind(euid);
		if (eid == -1)
			ERROREXIT;
		SRelement *elem = model.GetElement(eid);
		int tnv[8];
		int nv[4];
		int midv[4];
		nv[3] = -1;
		midv[3] = -1;
		i = 0;
		while (1)
		{
			if (!line.TokRead(tnv[i]))
				break;
			tnv[i] = NodeFind(tnv[i]);
			i++;
		}
		int nnodes = 4;
		if (i == 6)
			nnodes = 3;
		for (int j = 0; j < nnodes; j++)
			nv[j] = tnv[j];
		for (int j = 0; j < nnodes; j++)
			midv[j] = tnv[j + nnodes];

		int fId = elemFaceFind(elem, nv, gno);
		if (fId == -1)
			ERROREXIT;
		int conId = model.constraints.GetNum();
		SRconstraint* constraint = model.constraints.Add();
		//assign a coord id so breakout constraints are forced to be processed as penalty constraints, to avoid conflicts:
		constraint->coordId = 0; //0 is reserved for default gcs coord system
		constraint->breakout = true;
		constraint->type = faceCon;
		constraint->entityId = fId;
		SRface* face = model.GetFace(fId);
		face->constraintId = conId;
		int nnodesTotal = face->GetNumNodesTotal();
		constraint->enforcedDisplacementData.Allocate(nnodesTotal, 3);
		for (int i = 0; i < nnodesTotal; i++)
		{
			if (!model.inputFile.GetLine(line))
				ERROREXIT;
			tok = line.Token();//skip continuation character ("#")
			int idg;
			if (i < nnodes)
				idg = gno[i];
			else
				idg = face->midNodeMatch(midv[i - nnodes]);
			if (idg == -1)
				ERROREXIT;
			for (int dof = 0; dof < 3; dof++)
			{
				constraint->SetConstrainedDof(dof);
				double enfd;
				if (!line.TokRead(enfd))
					ERROREXIT;
				constraint->enforcedDisplacementData.Put(idg, dof, enfd);
			}
		}
	}
}

void SRinput::readSettings()
{
	//check for settings in SRsettings file:
	SRstring filename, line, basename, tail, tok;
	bool econSolve = false;
	bool highAcc = false;

	SRfile& settingsFile = model.settingsFile;
	bool settingsOpened = settingsFile.Open(SRinputMode);
	if (!settingsOpened)
		return;

	OUTPRINT("Stress Refine\n");
	OUTPRINT("\nCustom Settings for this run:");
	bool anyCustom = false;
	userPSettings = false;
	customPmax = 0;
	numSettingsStrings = 0;

	bool lineWasRead = false;

	SRvec3 breakoutOrigin;
	bool breakoutOriginAssigned = false;

	OUTPRINT("Stress Refine\n");
	OUTPRINT("\nCustom Settings for this run:");

	while (settingsOpened)
	{
		if (!lineWasRead)
		{
			if (!settingsFile.GetLine(line))
				break;
		}
		lineWasRead = false;
		if (line.isCommentOrBlank())
			continue;
		else if (line.CompareUseLength("parallel"))
		{
			model.SetNumThreads();
		}
		else if (line.CompareUseLength("NOUNITS"))
			model.useUnits = false;
		else if (line.CompareUseLength("stress_conversion"))
		{
			line.Token();
			line.TokRead(model.stressUnitConversion);
			model.stressUnitstr = line.Token();
		}
		else if (line.CompareUseLength("length_conversion"))
		{
			line.Token();
			line.TokRead(model.lengthUnitConversion);
			model.lengthUnitstr = line.Token();
		}
		else if (line.CompareUseLength("feasvmmax"))
		{
			line.Token();
			line.TokRead(model.feasvmmax);
			line.TokRead(model.feavmmaxpos.d[0]);
			line.TokRead(model.feavmmaxpos.d[1]);
			line.TokRead(model.feavmmaxpos.d[2]);
		}
		else if (line.CompareUseLength("localAdapt"))
		{
			model.localAdapt = true;
			//11/13/2017: testing shows this doesn't help much so ui for it is turned off in ui
			//revisit if ever try p2-itsolv
		}
		else if (line.CompareUseLength("originOfBreakout"))
		{
			line.Token();
			line.TokRead(breakoutOrigin.d[0]);
			line.TokRead(breakoutOrigin.d[1]);
			line.TokRead(breakoutOrigin.d[2]);
			breakoutOriginAssigned = true;
			model.breakoutByMat = false;
		}
		else if (line.CompareUseLength("breakout"))
		{
			model.breakout = true;
			line.Token();
			line.TokRead(model.post.maxRadNearOrigin);
			int atnode;
			line.TokRead(atnode);
			if (atnode != -1)
			{
				model.breakoutAtNode = true;
			}
		}
		else if (line.CompareUseLength("nbktElSacr"))
		{
			tok = line.Token();
			int nbktElSacr;
			line.TokRead(nbktElSacr);
			model.post.bktElSacrList.Allocate(nbktElSacr);
		}
		else if (line.CompareUseLength("needSoftSprings"))
		{
			model.needSoftSprings = true;
			LOGPRINT("use soft springs for stabilization");
		}
		else if (line == "PLIMITS")
		{
			tok = line.Token();
			int minp, maxp;
			line.TokRead(minp);
			if (minp < 2)
				minp = 2;
			line.TokRead(maxp);
			model.setMinPorder(minp);
			model.setMaxPorder(maxp);
			if (maxp == minp)
				model.maxPorder = maxp;
			OUTPRINT(" min p order: %d max p order %d\n", model.minPorder, model.maxPorderFinalAdapt);
			anyCustom = true;
			customPmax = true;
		}
		else if (line == "PMAX")
		{
			tok = line.Token();
			int maxp;
			line.TokRead(maxp);
			model.setMaxPorder(maxp);
			OUTPRINT(" min p order: %d max p order %d\n", model.minPorder, model.maxPorderFinalAdapt);
			anyCustom = true;
		}
		else if (line.CompareUseLength("MAXITERATIONS"))
		{
			tok = line.Token();
			int nits;
			line.TokRead(nits);
			model.setAdaptLoopMax(nits);
			OUTPRINT(" maximum adaptive loop iterations: %d\n", model.adaptLoopMax);
		}
		else if (line.CompareUseLength("maxPJump"))
		{
			tok = line.Token();
			line.TokRead(model.maxPJump);
		}
		else if (line.CompareUseLength("uniformP"))
		{
			model.pOrderUniform = true;
			OUTPRINT(" Uniform P-adaptivity\n");
			anyCustom = true;
		}
		else if (line.CompareUseLength("ERRORTOL"))
		{
			tok = line.Token();
			double errTol;
			line.TokRead(errTol);
			if (errTol < 0 || errTol > 25)
			{
				OUTPRINT(" error Tolerance must be in range 0 to 25 percent");
				ERROREXIT;
			}
			model.setErrorTolerance(errTol / 100.0);
			OUTPRINT(" Error Tolerance: %lg\n", model.GetErrorTolerance());
			anyCustom = true;
		}
		else if (line.CompareUseLength("FREEZESACRIFICIAL"))
		{
			tok = line.Token();
			int freezep;
			line.TokRead(freezep);
			model.SetFreezePorder(freezep);
			OUTPRINT(" Freeze Sacrificial Element Freeze P order to: %d\n", model.freezeSacrificalP);
			anyCustom = true;
		}
		else if (line.CompareUseLength("FREEZESAFEELEMENTS"))
		{
			tok = line.Token();
			int freezep;
			line.TokRead(freezep);
			model.SetMaxPorderLowStress(freezep);
			OUTPRINT("Freeze Low stress elements to p-order: %d\n", model.maxPorderLowStress);
			anyCustom = true;
		}
		else if (line.CompareUseLength("NOSACRIFICIAL"))
		{
			OUTPRINT("Sacrificial Element detection disabled\n");
			anyCustom = true;
			model.SetDetectSacr(false);
		}
		else if (line.CompareUseLength("DETECTTHINEDGES"))
		{
			model.detectThinElements = true;
			OUTPRINT("detect Thin Elements");
			anyCustom = true;
		}
		else if (line.CompareUseLength("LOWSTRESSTOL"))
		{
			double tol;
			tok = line.Token();
			line.TokRead(tol);
			model.errorChecker.SetLowStressTol(tol);
			OUTPRINT("Low Stress Tolerance: %lg\n", model.errorChecker.GetLowStressTol());
			anyCustom = true;
		}
		else if (line.CompareUseLength("FREEZESACREDGES"))
		{
			model.SetFreezeSacrificialEdges();
			OUTPRINT("Freeze p-order of all sacrifical elements");
			anyCustom = true;
		}
		else if (line.CompareUseLength("NOCHECKREENTRANT"))
		{
			model.errorChecker.checkReentrant = false;
			OUTPRINT("do not check for reentrant corners during sacrifical element detection");
			anyCustom = true;
		}
		else if (line.CompareUseLength("ALLCONSTRAINTSPENALTY"))
		{
			model.allConstraintsAsPenalty = true;
			OUTPRINT("calculate all constraints using penalty method\n");
			anyCustom = true;
		}
		else if (line.CompareUseLength("econSolve"))
		{
			econSolve = true;
			OUTPRINT(" economy solution setting");
			anyCustom = true;
		}
		else if (line.CompareUseLength("highAccuracy"))
		{
			highAcc = true;
			OUTPRINT(" high accuracy setting");
			anyCustom = true;
		}
		else if (line.CompareUseLength("noEnergySmooth"))
		{
			model.doEnergySmooth = false;
			OUTPRINT(" no energy smoothing of nodal loads");
			anyCustom = true;
		}
		else if (line.CompareUseLength("outputF06"))
		{
			model.outputf06 = true;
			OUTPRINT(" Output Ascii Nastran Stresses (.f06)");
			anyCustom = true;
		}
		else if (line.CompareUseLength("savebreakout"))
		{
			SRstring lineSav;
			lineSav = line;
			model.saveBreakout = true;
			SRBreakoutData *bdat = &model.saveBreakoutData;
			bdat->atMax = false;
			anyCustom = true;
			model.setMaxPorder(2);//don't run p-loop on the model being saved for breakout. It will be run on the breakout model
			tok = line.Token(); //skip saveBreakout keyword
			tok = line.Token();
			bdat->atMax = true;
			if (tok == "node")
			{
				line.TokRead(bdat->bktNode);
				if (bdat->bktNode != -1)
					bdat->atMax = false;
			}

			//create a folder for the saved breakout by appending _breakout to the filename
			char buf[20];
			SPRINTF(buf, "_breakout");
			SRstring breakoutTail;
			breakoutTail.Copy(model.fileNameTail);
			breakoutTail.Cat(buf);
			model.outdir.Left('\\', line);
			line.Cat("\\");
			line.Cat(breakoutTail);
			SRmachDep::CreateDir(line.getStr());
			//create .msh file for output in the new folder:
			filename.Copy(line);
			filename.Cat("\\");
			filename.Cat(breakoutTail);
			filename.Cat(".msh");
			bdat->saveBreakoutMshFile.filename = filename;
			bdat->saveBreakoutMshFile.Delete(); //clean up the folder if it was preexisting
		}
		else if (line.CompareUseLength("NoTopoFilter"))
			model.post.noTopoFilter = true;
		else if (line.CompareUseLength("ReadF06Results"))
		{
			tok = line.Token();
			model.f06File.filename = line.Token();
		}
		else if (line.CompareUseLength("partialDisplacements"))
			model.saveBreakoutData.fromPartialDispFile = true;
		else
			OUTPRINT("unrecognized setting: %s", line.getStr());
	}
	if (!anyCustom)
		OUTPRINT("   --None");

	settingsFile.Close();

	if (model.breakout && breakoutOriginAssigned)
		model.saveBreakoutData.origin.Copy(breakoutOrigin);

	if (econSolve)
	{
		model.maxPorder = model.maxPorderFinalAdapt = 5;
		model.adaptLoopMax = 2;
	}
	if (highAcc)
	{
		model.ErrorTolerance = 0.03;
		model.maxPorder = model.maxPorderFinalAdapt = 8;
		model.errorChecker.lowStressTolFinalAdapt = 0.5;
		model.errorChecker.lowStressTol = 0.5;
		model.SetFreezePorder(3);
		model.freezeSacrEdges = false;
	}
}

void SRinput::doNodalToFaceForces()
{
	//convert nodal forces to face forces
	int nfaceForce = nodalToFaceForces(countOnlyFlag);
	if (nfaceForce != 0)
	{
		int nforce0 = model.forces.GetNum();
		int nforceTotal = nforce0 + nfaceForce;
		model.forces.Allocate(nforceTotal);
		nodalToFaceForces(addTheForcesFlag);//false to add the face forces
	}

	if (nfaceForce == 0)
		return;

	//free the nodal forces which becaem inactive:
	for (int f = 0; f < model.GetNumForces(); f++)
	{
		SRforce* force = model.GetForce(f);
		if (force->type == inactiveForce)
			model.forces.Free(f);
	}
	model.forces.packNulls();
	//face tractions and face pressures may not have been read yet. reallocate for them:
	model.forces.Allocate(model.forces.GetNum() + nfacePressure + nfaceTraction);

	//correct force ids of face forces
	//(may have changed because of removal of inactive nodal forces):
	for (int i = 0; i < model.GetNumForces(); i++)
	{
		SRforce* f = model.GetForce(i);
		if (f == NULL)
			continue;
		int eid = f->GetEntityId();
		if (f->GetType() == faceForce)
			model.GetFace(eid)->forceIds.PushBack(i);
	}

	//error check. there should be no "orphan" nodal forces with pressure loads, or orphan midside nodes
	for (int f = 0; f < model.forces.GetNum(); f++)
	{
		SRforce* force = model.GetForce(f);
		bool okforce = true;
		if ((force->GetType() == nodalForce))
		{
			if (force->isPressure())
			{
				okforce = false;
				break;
			}
			int nid = force->entityId;
			SRnode* node = model.GetNode(nid);
			if (node->isMidSide())
			{
				okforce = false;
				break;
			}
		}
		if (!okforce)
		{
			SRnode* node = model.GetNode(force->entityId);
			OUTPRINT(" incorrect nodal force definition for node %d", node->userId);
			OUTPRINT(" pressure loads may not be applied to individual nodes, only to all nodes of faces");
			ERROREXIT;
		}
	}

	nodeForceStore.Free();

}

void SRinput::doNodalToFaceConstraints()
{
	if (model.GetNumFaces() == 0)
		return;
	//create face constraints from nodal constraints:
	int nfaceCon = nodalToFaceConstraints(countOnlyFlag);
	if (nfaceCon == 0)
		return;
	//reallocate model.constraints with room for the face constraints, then add the face constraints.
	int ncon0 = model.constraints.GetNum();
	int nconTotal = ncon0 + nfaceCon;
	model.constraints.Allocate(nconTotal);
	nodalToFaceConstraints(addTheConstraintsFlag);

	//fatal error if any nodes remain with more than one active constraint;
	for (int i = 0; i < nodeConStore.GetNum(); i++)
	{
		int numactive = 0;
		int n = nodeConStore.d[i].cons.GetNum();
		for (int j = 0; j < n; j++)
		{
			SRconstraint* con = nodeConStore.Get(i).GetCon(j);
			if (con->type == nodalCon)
				numactive++;
		}
		if (numactive > 1)
			ERROREXIT;
	}

	//get rid of any remaining constraints that aren't sensible:
	//nodal constraint on orphan node, or constrained midside node whose edge is not constrained.
	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		SRconstraint* con = model.GetConstraint(i);
		if (con->GetType() == nodalCon)
		{
			int nodeId = con->GetEntityId();
			SRnode* node = model.GetNode(nodeId);
			if (node->isOrphan())
				con->type = inactiveCon;
			else if (node->isMidSide())
			{
				int eid = model.GetNode(nodeId)->GetMidSideEdgeOwner();
				SRedge* edge = model.GetEdge(eid);
				if (edge->GetNode(0)->constraintId == -1 || edge->GetNode(1)->constraintId == -1)
					con->type = inactiveCon;
			}
		}
	}

	//delete the inactive nodal constraints:
	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		SRconstraint* con = model.GetConstraint(i);
		if (con->GetType() == inactiveCon)
		{
			int nodeId = con->GetEntityId();
			SRnode* node = model.GetNode(nodeId);
			node->constraintId = -1;
			model.constraints.Free(i);
		}
	}

	model.constraints.packNulls();

	//correct constraint ids of entities involved in constraints
	//(may have changed because of removal of inactive nodal constraints):
	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		SRconstraint* con = model.GetConstraint(i);
		if (con == NULL)
			continue;
		int eid = con->GetEntityId();
		if (con->GetType() == nodalCon)
			model.GetNode(eid)->SetConstraintId(i);
		if (con->GetType() == faceCon)
			model.GetFace(eid)->SetConstraintId(i);
		if (model.allConstraintsAsPenalty && con->isGcs())
			con->coordId = 0; //0 is reserved for default gcs coord system
	}

	nodeConStore.Free();

}

void SRinput::dumpModel(bool meshonly)
{
	bool dodump = false;
	if (!dodump)
		return;
	model.dumpFileOpen();
	model.dumpFile.PrintLine("nodes");
	for (int n = 0; n < model.GetNumNodes(); n++)
	{
		SRnode* node = model.GetNode(n);
		node->dumpData();
	}
	model.dumpFile.PrintLine("edges");
	for (int n = 0; n < model.GetNumEdges(); n++)
	{
		SRedge* edge = model.GetEdge(n);
		edge->dumpData();
	}
	if (!meshonly)
	{
		model.dumpFile.PrintLine("faces");
		for (int f = 0; f < model.faces.GetNum(); f++)
		{
			SRface* face = model.GetFace(f);
			//if (face->IsBoundaryFace())
				face->dumpData();
		}
	}
	model.dumpFile.PrintLine("elements");
	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		elem->dumpData();
	}

	if (meshonly)
		return;
	model.dumpFile.PrintLine("forces");
	for (int f = 0; f < model.forces.GetNum(); f++)
	{
		SRforce* force = model.GetForce(f);
		force->dumpData();
	}
	model.dumpFile.PrintLine("constraints");
	for (int c = 0; c < model.GetNumConstraints(); c++)
	{
		SRconstraint* con = model.GetConstraint(c);
		con->dumpData();
	}
	model.dumpFileClose();
}



