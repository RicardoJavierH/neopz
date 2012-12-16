//$Id: pzpostprocanalysis.cpp,v 1.10 2010-11-23 18:58:35 diogo Exp $
#include "pzanalysis.h"
#include "pzpostprocanalysis.h"
#include "pzpostprocmat.h"
#include "pzcompelpostproc.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzvec.h"
#include "tpzautopointer.h"
#include "tpzcompmeshreferred.h"
#include "pzstring.h"
#include "pzelastoplasticanalysis.h"
#include "pzcreateapproxspace.h"

#include <map>
#include <set>
#include <stdio.h>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr PPAnalysisLogger(Logger::getLogger("analysis.postproc"));
#endif

using namespace std;

TPZPostProcAnalysis::TPZPostProcAnalysis() : fpMainAnalysis(NULL)
{	
}

TPZPostProcAnalysis::TPZPostProcAnalysis(TPZAnalysis * pRef):TPZAnalysis(), fpMainAnalysis(pRef)
{
//    
//    TPZCompMesh* pcMainMesh = fpMainAnalysis->Mesh();
//    //TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(pcMainMesh);
//    
//    TPZGeoMesh * pgmesh = pcMainMesh->Reference();
//
//    
//    TPZCompMeshReferred * pcPostProcMesh = new TPZCompMeshReferred(pgmesh);
//    
//    
//    //TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(pcPostProcMesh);
//    
//    //pcPostProcMesh->ApproxSpace().SetAllCreateFunctionsDiscontinuousReferred();
//    //SetAllCreateFunctionsDiscontiunous(pcPostProcMesh);
//    //pcPostProcMesh->LoadReferred(pcMainMesh);
//    
//    pcPostProcMesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
//    TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(pcPostProcMesh);
//
//    fCompMesh = pcPostProcMesh;
    
    

    TPZCompMesh* pcMainMesh = fpMainAnalysis->Mesh();
    
    TPZGeoMesh * pgmesh = pcMainMesh->Reference();

   // TPZPostProcAnalysis::SetAllCreateFunctionsPostProc();
    
    TPZCompMeshReferred * pcPostProcMesh = new TPZCompMeshReferred(pgmesh);
    
    fCompMesh = pcPostProcMesh;

    TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(pcPostProcMesh);
    
    
    
}

TPZPostProcAnalysis::~TPZPostProcAnalysis()
{
}

void TPZPostProcAnalysis::SetPostProcessVariables(TPZVec<int> & matIds, TPZVec<std::string> &varNames)
{
	int i, j, nMat, matNumber;


	TPZCompMesh * pcMainMesh = fpMainAnalysis->Mesh();
	
	TPZGeoMesh * pgmesh = pcMainMesh->Reference();

	TPZCompMeshReferred * pcPostProcMesh = dynamic_cast<TPZCompMeshReferred *>(this->Mesh());
	
	TPZStack<int> avlMatIds;
	int nel = pgmesh->NElements();	
	for(i = 0; i < nel; i++)
	{
		int matId = pgmesh->ElementVec()[i]->MaterialId();
		int isMatPostProc = 0;
		int isMatAvl = 0;
		j = 0;
		nMat = matIds.NElements();
		while(j < nMat && !isMatPostProc)
		{
			if(matId == matIds[j])isMatPostProc = 1;
			j++;
		}
		
		if(!isMatPostProc)
		{
			nMat = avlMatIds.NElements();
			j = 0;
			while(j < nMat && !isMatAvl)
			{
				if(matId == avlMatIds[j])isMatAvl = 1;
				j++;
			}
			
			if(!isMatAvl)
			{
				avlMatIds.Push(matId);
			}
		}
	}
	
	nMat = matIds.NElements();
	for(i = 0; i < nMat; i++)
	{
		TPZMaterial * pmat = pcMainMesh->FindMaterial(matIds[i]);
		if(!pmat)
		{
			PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::SetPostProcessVariables() material Id " << matIds[i] << " not found in original mesh!\n";
			continue;
		}
		
		TPZPostProcMat * pPostProcMat = new TPZPostProcMat(matIds[i]);
		
		pPostProcMat->SetPostProcessVarIndexList(varNames,pmat);
		
		matNumber = pcPostProcMesh->InsertMaterialObject(pPostProcMat);

	}
	
	//pcPostProcMesh->AutoBuild();
	//pcPostProcMesh->AutoBuildDisc();
	AutoBuildDisc();
	
	pcPostProcMesh->LoadReferred(pcMainMesh);
}

void TPZPostProcAnalysis::AutoBuildDisc() 
{
	TPZAdmChunkVector<TPZGeoEl *> &elvec = Mesh()->Reference()->ElementVec();
	int i, nelem = elvec.NElements();
	int neltocreate = 0;
	int index;
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			neltocreate++;
		}
	}
	std::set<int> matnotfound;
	int nbl = Mesh()->Block().NBlocks();
	if(neltocreate > nbl) Mesh()->Block().SetNBlocks(neltocreate);
	Mesh()->Block().SetNBlocks(nbl);
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			int matid = gel->MaterialId();
			TPZMaterial * mat = Mesh()->FindMaterial(matid);
			if(!mat)
			{
				matnotfound.insert(matid);
				continue;
			}
			int printing = 0;
			if (printing) {
				gel->Print(cout);
			}
			
			//if(!gel->Reference() && gel->NumInterfaces() == 0)
			//gel->CreateCompEl(*Mesh(),index);
			//gel->ResetReference();
			Mesh()->CreateCompEl(gel,index);
            gel->ResetReference();
			
		}
	}
	
	Mesh()->InitializeBlock();
#ifdef LOG4CXX
    if(PPAnalysisLogger->isDebugEnabled())
    {
        std::stringstream sout;
        Mesh()->Print(sout);
        LOGPZ_DEBUG(PPAnalysisLogger, sout.str())
    }
#endif
	if(matnotfound.size())
	{
		std::cout << "Malha post proc was created without these materials ";
		std::set<int>::iterator it;
		for(it = matnotfound.begin(); it!= matnotfound.end(); it++)
		{
			std::cout << *it << " ";
		}
		std::cout << std::endl;
	}
	
}

void TPZPostProcAnalysis::Assemble()
{
   PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::Assemble() should never be called\n";
}

void TPZPostProcAnalysis::Solve(){
   PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::Solve() should never be called\n";
}

void TPZPostProcAnalysis::TransferSolution()
{

    
    TPZAnalysis::AssembleResidual();
    fSolution = Rhs();
    TPZAnalysis::LoadSolution();
    

	
}


#include "pzintel.h"

#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "tpzpoint.h"

#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "tpzline.h"

#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "tpztriangle.h"

#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"

#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "tpzprism.h"

#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "tpztetrahedron.h"

#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "tpzpyramid.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "tpzcube.h"

#include "pzelctemp.h"


/*void TPZPostProcAnalysis::SetAllCreateFunctionsPostProc()
{
	pzgeom::TPZGeoPoint::fp = TPZPostProcAnalysis::CreatePointEl;
	pzgeom::TPZGeoLinear::fp = TPZPostProcAnalysis::CreateLinearEl;
	pzgeom::TPZGeoQuad::fp = TPZPostProcAnalysis::CreateQuadEl;
	pzgeom::TPZGeoTriangle::fp = TPZPostProcAnalysis::CreateTriangleEl;
	pzgeom::TPZGeoPrism::fp = TPZPostProcAnalysis::CreatePrismEl;
	pzgeom::TPZGeoTetrahedra::fp = TPZPostProcAnalysis::CreateTetraEl;
	pzgeom::TPZGeoPyramid::fp = TPZPostProcAnalysis::CreatePyramEl;
	pzgeom::TPZGeoCube::fp = TPZPostProcAnalysis::CreateCubeEl;


}*/


void TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(TPZCompMesh *cmesh)
{
    
    TPZManVector<TCreateFunction,10> functions(8);
    
    functions[EPoint] = &TPZPostProcAnalysis::CreatePointEl;
    functions[EOned] = TPZPostProcAnalysis::CreateLinearEl;
    functions[EQuadrilateral] = TPZPostProcAnalysis::CreateQuadEl;
    functions[ETriangle] = TPZPostProcAnalysis::CreateTriangleEl;
    functions[EPrisma] = TPZPostProcAnalysis::CreatePrismEl;
    functions[ETetraedro] = TPZPostProcAnalysis::CreateTetraEl;
    functions[EPiramide] = TPZPostProcAnalysis::CreatePyramEl;
    functions[ECube] = TPZPostProcAnalysis::CreateCubeEl;
    cmesh->ApproxSpace().SetCreateFunctions(functions);
  /*
    functions[EPoint] = TPZCompElDisc::CreateDisc;
    functions[EOned] = TPZCompElDisc::CreateDisc;
    functions[ETriangle] = TPZCompElDisc::CreateDisc;
    functions[EQuadrilateral] = TPZCompElDisc::CreateDisc;
    functions[ETetraedro] = TPZCompElDisc::CreateDisc;
    functions[EPiramide] = TPZCompElDisc::CreateDisc;
    functions[EPrisma] = TPZCompElDisc::CreateDisc;
    functions[ECube] = TPZCompElDisc::CreateDisc;
    cmesh->ApproxSpace().SetCreateFunctions(functions);
   */ 
}



using namespace pzshape;

TPZCompEl *TPZPostProcAnalysis::CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc< TPZIntelGen<TPZShapePoint> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeLinear> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeQuad> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeTriang> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeCube> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc< TPZIntelGen<TPZShapePrism> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapePiram> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeTetra> >(mesh,gel,index);
	return NULL;
}


TPZCompEl * TPZPostProcAnalysis::CreatePostProcDisc(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElPostProc< TPZCompElDisc > (mesh,gel,index);
}
