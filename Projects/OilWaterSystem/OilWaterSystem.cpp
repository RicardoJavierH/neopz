#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgnode.h"
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzgeoel.h"
#include "pzmatrix.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "pzconvectionproblem.h"
#include "pzmultiphase.h"
#include "pzl2projection.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "fad.h"

#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzbndcond.h"
#include "pzlog.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"

#include "pzskylstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzpoisson3d.h"
#include <time.h>
#include <stdio.h>

// Using Log4cXX as logging tool
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.multiphase.data"));
#endif
//
// End Using Log4cXX as logging tool


TPZCompMesh *ComputationalMeshPseudopressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *ComputationalMeshBulkflux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *ComputationalMeshWaterSaturation(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *ComputationalMeshMultiphase(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
TPZCompMesh *L2ProjectionQ(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);
TPZCompMesh *L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);

void SolveSyst(TPZAnalysis &an, TPZCompMesh *cmesh);
void InitialFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void InitialPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void InitialSaturation(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void PosProcessBulkflux(TPZAnalysis &an, std::string plotfile);
void PosProcessL2(TPZAnalysis &an, std::string plotfile);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void UniformRefinement(TPZGeoMesh *gMesh, int nh);
void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
void SolveSystemTransient(REAL deltaT,REAL maxTime, TPZAnalysis *an, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics);
void CheckConvergence(TPZFMatrix<STATE> SoltUattn,TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics);
void ComputeResidual(TPZFMatrix<STATE> SoltUattn,REAL alpha, TPZFMatrix<STATE> DeltaU, TPZFMatrix<STATE> &ResAlpha, TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics);

bool ftriang = false;

int main()
{
	
#ifdef LOG4CXX
	InitializePZLOG("../OilWaterLog4cxx.cfg");
#endif		
	
	
	#ifdef LOG4CXX
		InitializePZLOG();
	#endif		

	//	Reading mesh
	std::string GridFileName;	
	GridFileName = "OilWaterSystemUnit.dump";	
//	GridFileName = "OilWaterSystemUnitTwo.dump";
//	GridFileName = "OilWaterSystemUnitOneHRef.dump";	
//	GridFileName = "OilWaterSystemUnitTwoHRef.dump";
	
	
	bool twoMaterial = false;
	TPZReadGIDGrid GeometryInfo;
	GeometryInfo.SetfDimensionlessL(1.0);
	TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);	
	
	{
		//	Print Geometrical Base Mesh
		ofstream argument("GeometicMesh.txt");
		gmesh->Print(argument);
		ofstream Dummyfile("GeometricMesh.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
	}
	
	std::vector<double> dd(2,0);
	
	
	int Href = 3;
	int div = 0;		
	int POrderBulkFlux = 1;	
	int POrderPseudopressure = 1;
	int POrderWaterSaturation = 0;	
	
	UniformRefinement(gmesh, Href);	
	
	{
		//	Print Geometrical Base Mesh
		ofstream argument("RefGeometicMesh.txt");
		gmesh->Print(argument);
		ofstream Dummyfile("RefGeometricMesh.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
//		gmesh->BuildConnectivity();
	}	

	
	 // Computational meshes	
	
	//	First computational mesh
	TPZCompMesh * CMeshBulkflux = ComputationalMeshBulkflux(gmesh, POrderBulkFlux);
	//	Print Second computational mesh
	ofstream ArgumentBulkflux("ComputationalMeshForBulkflux.txt");
	CMeshBulkflux->Print(ArgumentBulkflux);	
	
	//	Second computational mesh
	TPZCompMesh * CMeshPseudopressure = ComputationalMeshPseudopressure(gmesh, POrderPseudopressure);
	//	Print First computational mesh
	ofstream ArgumentPseudopressure("ComputationalMeshForPseudopressure.txt");
	CMeshPseudopressure->Print(ArgumentPseudopressure);
	
	//	Third computational mesh
	TPZCompMesh * CMeshWaterSaturation = ComputationalMeshWaterSaturation(gmesh, POrderWaterSaturation);
	//	Print Second computational mesh
	ofstream ArgumentWaterSaturation("ComputationalMeshForWaterSaturation.txt");
	CMeshWaterSaturation->Print(ArgumentWaterSaturation);
	
	TPZAnalysis Anbulkflux(CMeshBulkflux);
    std::string outputfile1;
	outputfile1 = "SolutionBulkflux";
    std::stringstream outputfiletemp1;
    outputfiletemp1 << outputfile1 << ".vtk";
    std::string plotfilebuklflux = outputfiletemp1.str();
	TPZFMatrix<STATE> InitialQSolution = Anbulkflux.Solution();	
	int rwosQ= InitialQSolution.Rows();
	for (int i=0; i < InitialQSolution.Rows(); i++) 
	{
		InitialQSolution(i)=0.0;
	}
	
	
	REAL initqx = 0.1;
	
	
//	InitialQSolution(8)=1.0;	
//	InitialQSolution(9)=0.0;		
//	InitialQSolution(10)=0.0;
//	InitialQSolution(11)=0.0;	
//
//	InitialQSolution(0)=initqx;
//	InitialQSolution(1)=initqx;
//	InitialQSolution(2)=0.0;
//	InitialQSolution(3)=0.0;
//	InitialQSolution(4)=0.0;
//	InitialQSolution(5)=0.0;
//	InitialQSolution(6)=-initqx;
//	InitialQSolution(7)=-initqx;
//	InitialQSolution(8)=0.0;	
//	InitialQSolution(9)=0.0;		
//	InitialQSolution(10)=0.0;
//	InitialQSolution(11)=0.0;		
//	
//	InitialQSolution(0)=0.0;
//	InitialQSolution(1)=0.0;
//	InitialQSolution(2)=initqx;
//	InitialQSolution(3)=initqx;
//	InitialQSolution(4)=0.0;
//	InitialQSolution(5)=0.0;
//	InitialQSolution(6)=0.0;
//	InitialQSolution(7)=0.0;
//	InitialQSolution(8)=-initqx;//	
//	InitialQSolution(9)=-initqx;//		
//	InitialQSolution(10)=0.0;
//	InitialQSolution(11)=0.0;
//	InitialQSolution(12)=0.0;
//	InitialQSolution(13)=0.0;
//	InitialQSolution(14)=-initqx;//
//	InitialQSolution(15)=-initqx;//
//	InitialQSolution(16)=0.0;
//	InitialQSolution(17)=0.0;
//	InitialQSolution(18)=0.0;
//	InitialQSolution(19)=0.0;
//	InitialQSolution(20)=0.0;
//	InitialQSolution(21)=0.0;
//	InitialQSolution(22)=0.0;	
//	InitialQSolution(23)=0.0;
//	InitialQSolution(24)=initqx;
//	InitialQSolution(25)=initqx;	
//
//	InitialQSolution(26)=0.0;
//	InitialQSolution(27)=0.0;
//	InitialQSolution(28)=0.0;
//	InitialQSolution(29)=0.0;	
//	InitialQSolution(30)=initqx;
//	InitialQSolution(31)=initqx;
//
//	InitialQSolution(32)=0.0;
//	InitialQSolution(33)=0.0;
//	InitialQSolution(34)=-initqx;//
//	InitialQSolution(35)=-initqx;//
//	InitialQSolution(36)=0.0;
//	InitialQSolution(37)=0.0;	
//	InitialQSolution(38)=0.0;
//	InitialQSolution(39)=0.0;
	
	Anbulkflux.LoadSolution(InitialQSolution);
	int num= InitialQSolution.Rows();


	TPZVec<STATE> soliniQ(InitialQSolution.Rows(),0.0);
    TPZCompMesh  * cmeshQL2 = L2ProjectionQ(gmesh, POrderBulkFlux, soliniQ);
    
    TPZAnalysis anQL2(cmeshQL2);
    SolveSyst(anQL2, cmeshQL2);
//    Anbulkflux.LoadSolution(anQL2.Solution());	
	PosProcessBulkflux(Anbulkflux,plotfilebuklflux);
	


	TPZAnalysis AnPressure(CMeshPseudopressure);
	
    std::string outputfile2;
	outputfile2 = "SolutionPressure";
    std::stringstream outputfiletemp2;
    outputfiletemp2 << outputfile2 << ".vtk";
    std::string plotfilePressure = outputfiletemp2.str();
	TPZFMatrix<STATE> InitialPSolution = AnPressure.Solution();	
	for (int i=0; i < InitialPSolution.Rows(); i++) {
		InitialPSolution(i)=-0.0;
	}
	
	
	TPZVec<STATE> solini(InitialPSolution.Rows(),0.0);
    TPZCompMesh  * cmeshL2 = L2ProjectionP(gmesh, POrderPseudopressure, solini);
    
    TPZAnalysis anL2(cmeshL2);
    SolveSyst(anL2, cmeshL2);
	
//	InitialPSolution(0)=0.0;
//	InitialPSolution(1)=0.0;
//	InitialPSolution(2)=-6.36396103;
//	InitialPSolution(3)=5.50000000;
		
	
//	InitialPSolution(0)=5.5;
//	InitialPSolution(1)=14.5;
//	InitialPSolution(2)=14.5;
//	InitialPSolution(3)=5.5;	
    AnPressure.LoadSolution(anL2.Solution());
	

//	int TotalComputationalElements	= CMeshPseudopressure->NElements();		
//	for ( int el = 0; el < TotalComputationalElements ; el++ )
//	{
//		
//		TPZCompEl * ComputationalElement = CMeshPseudopressure->ElementVec()[el];
//		TPZGeoEl * GeoemtricElement;
//		GeoemtricElement = ComputationalElement->Reference();
//		int NumberNodeCorners = ComputationalElement->Reference()->NCornerNodes();
//		
//		for (int ConnectIndex = 0; ConnectIndex < NumberNodeCorners ; ConnectIndex++) 
//		{
//			if (NumberNodeCorners == 2) {
//				continue;
//			}
//			else 
//			{
//				TPZConnect & ThisConnect = ComputationalElement->Connect(ConnectIndex);
//				TPZVec <REAL> Coordenates(3,0.0);
//				REAL InitialValue;			
//				Coordenates[0] = GeoemtricElement->NodePtr(ConnectIndex)->Coord(0);
//				Coordenates[1] = GeoemtricElement->NodePtr(ConnectIndex)->Coord(1);
//				Coordenates[2] = GeoemtricElement->NodePtr(ConnectIndex)->Coord(2);						
//				
//				int Sequence = ThisConnect.SequenceNumber();
//				int BlockPos = CMeshPseudopressure->Block().Position(Sequence);
//				
//				InitialPSolution(BlockPos)= 1.0e6;				
//			}
//
//
//		}
//	}
	
//	AnPressure.LoadSolution(InitialPSolution);	
	PosProcessL2(AnPressure,plotfilePressure);	
	

	TPZAnalysis AnSaturation(CMeshWaterSaturation);
    std::string outputfile3;
	outputfile3 = "SolutionSaturation";
    std::stringstream outputfiletemp3;
    outputfiletemp3 << outputfile3 << ".vtk";
    std::string plotfileSaturation = outputfiletemp3.str();
	TPZFMatrix<STATE> InitialSSolution = AnSaturation.Solution();	
	for (int i=0; i < InitialSSolution.Rows(); i++) {
		InitialSSolution(i)=0.0;
	}
	AnSaturation.LoadSolution(InitialSSolution);
	PosProcessL2(AnSaturation,plotfileSaturation);		
	
	//	This is so rare!!
	if (twoMaterial) 
	{
		gmesh->AddInterfaceMaterial(1,2, 1);
		gmesh->AddInterfaceMaterial(2,1, 1);
	}
//	
//	
//	RefinUniformElemComp(CMeshBulkflux, 2);
//	RefinUniformElemComp(CMeshPseudopressure, 2);
//	RefinUniformElemComp(CMeshWaterSaturation, 2);	
	
    //	Multiphysics Mesh
    TPZVec<TPZCompMesh *> meshvec(3);
    meshvec[0] = CMeshBulkflux;
    meshvec[1] = CMeshPseudopressure;	
    meshvec[2] = CMeshWaterSaturation;
	
    
    TPZCompMesh * MultiphysicsMesh = ComputationalMeshMultiphase(gmesh,meshvec);
    ofstream ArgumentMultiphysic("MultiphysicsMesh.txt");
    MultiphysicsMesh->Print(ArgumentMultiphysic);
	
//	TPZBuildMultiphysicsMesh::UniformRefineCompMesh(MultiphysicsMesh, Href, true)	
	
    TPZAnalysis *MultiphysicsAn = new TPZAnalysis(MultiphysicsMesh);
	int	Nthreads = 2;
	//    //TPZFStructMatrix matsk(mphysics);
//    TPZSkylineStructMatrix matsk(MultiphysicsMesh);	
	TPZSkylineNSymStructMatrix matsk(MultiphysicsMesh);	
//	TPZFStructMatrix matsk(MultiphysicsMesh);
	MultiphysicsAn->SetStructuralMatrix(matsk);
	MultiphysicsAn->StructMatrix()->SetNumThreads(Nthreads);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELU); 					
	MultiphysicsAn->SetSolver(step);
		
	
//	an.Run();	
//	an.Solve();	
	
//	TPZFMatrix<STATE> Initialsolution = MultiphysicsAn.Solution();
    std::string outputfile;
	outputfile = "TransientSolution";
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
	PosProcessMultphysics(meshvec,MultiphysicsMesh,*MultiphysicsAn,plotfile);		
	
	int datar = MultiphysicsAn->Solution().Rows();
	

    REAL deltaT = 0.8; //second
    REAL maxTime = 8.0;
    SolveSystemTransient(deltaT, maxTime, MultiphysicsAn, meshvec, MultiphysicsMesh);

	return 0;
	
}	
	
	TPZCompMesh *ComputationalMeshBulkflux(TPZGeoMesh *gmesh, int pOrder)
	{
		/// criar materiais
		int dim = 2;
		int matId1 = 1;
		int matId2 = 2;		
		
		TPZMatPoisson3d *material1;
		material1 = new TPZMatPoisson3d(matId1,dim);
		TPZMaterial * mat1(material1);
		material1->NStateVariables();

//		TPZMatPoisson3d *material2;
//		material2 = new TPZMatPoisson3d(matId2,dim);
//		TPZMaterial * mat2(material2);
//		material2->NStateVariables();		
		
		TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
		cmesh->SetDefaultOrder(pOrder);
		cmesh->SetDimModel(dim);
		cmesh->InsertMaterialObject(mat1);
//		cmesh->InsertMaterialObject(mat2);		
				
		TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
		TPZMaterial * BCond2 = material1->CreateBC(mat1,2,0, val1, val2);
		TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);		
//		TPZMaterial * BCond4 = material2->CreateBC(mat2,4,0, val1, val2);
		TPZMaterial * BCond4 = material1->CreateBC(mat1,4,0, val1, val2);
//		TPZMaterial * BCond6 = material2->CreateBC(mat2,6,0, val1, val2);
		TPZMaterial * BCond5 = material1->CreateBC(mat1,5,0, val1, val2);
//		TPZMaterial * BCond8 = material2->CreateBC(mat2,8,0, val1, val2);		
		
		cmesh->SetAllCreateFunctionsHDiv();
		cmesh->InsertMaterialObject(BCond2);
		cmesh->InsertMaterialObject(BCond3);		
//		cmesh->InsertMaterialObject(BCond4);
		cmesh->InsertMaterialObject(BCond4);
//		cmesh->InsertMaterialObject(BCond6);
		cmesh->InsertMaterialObject(BCond5);
//		cmesh->InsertMaterialObject(BCond8);
		
			
		cmesh->AutoBuild();
		cmesh->AdjustBoundaryElements();
		cmesh->CleanUpUnconnectedNodes();
		
		return cmesh;
	}

	TPZCompMesh *ComputationalMeshPseudopressure(TPZGeoMesh *gmesh, int pOrder)
	{
		/// criar materiais
		int dim = 2;
		int matId1 = 1;
		int matId2 = 2;
		
		TPZMatPoisson3d *material1;
		material1 = new TPZMatPoisson3d(matId1,dim);
		TPZMaterial * mat1(material1);
		material1->NStateVariables();
		
//		TPZMatPoisson3d *material2;
//		material2 = new TPZMatPoisson3d(matId2,dim);
//		TPZMaterial * mat2(material2);
//		material2->NStateVariables();
		
		TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
		cmesh->SetDefaultOrder(pOrder);
		cmesh->SetDimModel(dim);
		cmesh->InsertMaterialObject(mat1);
//		cmesh->InsertMaterialObject(mat2);		
		
		
		TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
		TPZMaterial * BCond2 = material1->CreateBC(mat1,2,0, val1, val2);
		TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);		
//		TPZMaterial * BCond4 = material2->CreateBC(mat2,4,0, val1, val2);
		TPZMaterial * BCond4 = material1->CreateBC(mat1,4,0, val1, val2);
//		TPZMaterial * BCond6 = material2->CreateBC(mat2,6,0, val1, val2);
		TPZMaterial * BCond5 = material1->CreateBC(mat1,5,0, val1, val2);
//		TPZMaterial * BCond8 = material2->CreateBC(mat2,8,0, val1, val2);		
		
		cmesh->SetAllCreateFunctionsDiscontinuous(); // L2 approximation space
		cmesh->InsertMaterialObject(BCond2);
		cmesh->InsertMaterialObject(BCond3);		
//		cmesh->InsertMaterialObject(BCond4);
		cmesh->InsertMaterialObject(BCond4);
//		cmesh->InsertMaterialObject(BCond6);
		cmesh->InsertMaterialObject(BCond5);
//		cmesh->InsertMaterialObject(BCond8);
		
		
		cmesh->SetDefaultOrder(pOrder);
		cmesh->SetDimModel(dim);
		cmesh->AutoBuild();			
		
		
		///inserir connect da pressao
		int ncon = cmesh->NConnects();
		for(int i=0; i<ncon; i++)
		{
			TPZConnect &newnod = cmesh->ConnectVec()[i];
			//newnod.SetPressure(true);
			newnod.SetLagrangeMultiplier(1);
		}
		
		///set order total da shape
		int nel = cmesh->NElements();
		for(int i=0; i<nel; i++){
			TPZCompEl *cel = cmesh->ElementVec()[i];
			TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
			celdisc->SetConstC(1.);
			celdisc->SetCenterPoint(0, 0.);
			celdisc->SetCenterPoint(1, 0.);
			celdisc->SetCenterPoint(2, 0.);
			celdisc->SetTrueUseQsiEta();
			if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
			{
				if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
				else celdisc->SetTensorialShape();
			}
		}
//		
//		
//#ifdef DEBUG
//		int ncel = cmesh->NElements();
//		for(int i =0; i<ncel; i++){
//			TPZCompEl * compEl = cmesh->ElementVec()[i];
//			if(!compEl) continue;
//			TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
//			if(facel)DebugStop();
//			
//		}
//#endif		
		
		cmesh->AdjustBoundaryElements();
		cmesh->CleanUpUnconnectedNodes();
		
		
		return cmesh;
	}
	

	TPZCompMesh *ComputationalMeshWaterSaturation(TPZGeoMesh * gmesh, int pOrder)
	{
		/// criar materiais
		int dim = 2;
		int matId1 = 1;
		int matId2 = 2;	
		
		TPZMatConvectionProblem *material1 = new TPZMatConvectionProblem(matId1,dim);
		TPZMaterial * mat1(material1);
		
//		TPZMatConvectionProblem *material2 = new TPZMatConvectionProblem(matId2,dim);
//		TPZMaterial * mat2(material2);		
		
		TPZVec<REAL> convdir(dim,0.);
		convdir[0]=1.;
		REAL flux = 0.;
		REAL rho = 1.;
		
		material1->SetParameters(rho,convdir);
		material1->SetInternalFlux(flux);
		material1->NStateVariables();
		
//		material2->SetParameters(rho,convdir);
//		material2->SetInternalFlux(flux);
//		material2->NStateVariables();		
		
		TPZCompEl::SetgOrder(pOrder);
		TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
		cmesh->SetDimModel(dim);
		cmesh->InsertMaterialObject(mat1);
//		cmesh->InsertMaterialObject(mat2);		
		
		TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
		TPZMaterial * BCond2 = material1->CreateBC(mat1,2,0, val1, val2);
		TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);		
//		TPZMaterial * BCond4 = material2->CreateBC(mat2,4,0, val1, val2);
		TPZMaterial * BCond4 = material1->CreateBC(mat1,4,0, val1, val2);
//		TPZMaterial * BCond6 = material2->CreateBC(mat2,6,0, val1, val2);
		TPZMaterial * BCond5 = material1->CreateBC(mat1,5,0, val1, val2);
//		TPZMaterial * BCond8 = material2->CreateBC(mat2,8,0, val1, val2);
		
		cmesh->SetAllCreateFunctionsDiscontinuous();//  L2 approximation space
		cmesh->InsertMaterialObject(BCond2);
		cmesh->InsertMaterialObject(BCond3);		
//		cmesh->InsertMaterialObject(BCond4);
		cmesh->InsertMaterialObject(BCond4);
//		cmesh->InsertMaterialObject(BCond6);
		cmesh->InsertMaterialObject(BCond5);
//		cmesh->InsertMaterialObject(BCond8);
				
		///set order total da shape
		int nel = cmesh->NElements();
		for(int i=0; i<nel; i++){
			TPZCompEl *cel = cmesh->ElementVec()[i];
			TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
			if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
			{
				if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
				else celdisc->SetTensorialShape();
			}
		}
				
		cmesh->AdjustBoundaryElements();
		cmesh->CleanUpUnconnectedNodes();		
		cmesh->AutoBuild();
		
		
		return cmesh;
	}
	
	TPZCompMesh *ComputationalMeshMultiphase(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
		
		
		//Creating computational mesh for multiphysic elements
		gmesh->ResetReference();
		TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
//		mphysics->SetAllCreateFunctionsMultiphysicElem();
		
		int dim =2;
		int matId1 = 1;
		int matId2 = 2;
		// Setting data
		
		TPZMultiphase *material1 = new TPZMultiphase(matId1,dim);
//		TPZMultiphase *material2 = new TPZMultiphase(matId2,dim);
		
		REAL deltaT = 0.1;
		REAL maxTime = 0.1;
		
		material1->SetTimeStep(1.0);
		material1->fnewWS=true;
		material1->SetTheta(0.5);
		material1->LoadKMap("Permeabilities.txt");
		material1->SetYorN(false);
//		material2->SetTimeStep(1.0);
//		material2->SetTheta(0.5);

//		material1->SetCurrentState();
//		material2->SetCurrentState();
//		
//		material1->SetLastState();
//		material2->SetLastState();		
		
//		material1->SetPermeability , densitys, etc ...

	
		TPZMaterial *mat1(material1);
		mphysics->InsertMaterialObject(mat1);
		mphysics->SetDimModel(dim);
		
//		TPZMaterial *mat2(material2);
//		mphysics->InsertMaterialObject(mat2);
//		mphysics->SetDimModel(dim);
		
		
		TPZFMatrix<STATE> val1(4,2,0.), val2(4,1,0.);

		val2(0,0)=0.0;// qx
		val2(1,0)=0.0;// qy
		val2(2,0)=0.1;// P
		val2(3,0)=1.0;// S	
		TPZMaterial * BCond5 = material1->CreateBC(mat1,5,1, val1, val2);
		
		val2(0,0)=0.0;// qx
		val2(1,0)=0.0;// qy
		val2(2,0)=0.0;// P
		val2(3,0)=0.0;// S			
		TPZMaterial * BCond3 = material1->CreateBC(mat1,2,3, val1, val2);		
		//		TPZMaterial * BCond4 = material2->CreateBC(mat2,4,0, val1, val2);
		
		val2(0,0)=0.0;// qx
		val2(1,0)=0.0;// qy
		val2(2,0)=0.0;// P
		val2(3,0)=0.0;// S			
		TPZMaterial * BCond2 = material1->CreateBC(mat1,4,3, val1, val2);
		//		TPZMaterial * BCond4 = material1->CreateBC(mat1,4,2, val1, val2);		
		//		TPZMaterial * BCond6 = material2->CreateBC(mat2,6,3, val1, val2);
		
		val2(0,0)=0.0;// qx
		val2(1,0)=0.0;// qy
		val2(2,0)=0.0;// P
		val2(3,0)=0.0;// S			
		TPZMaterial * BCond4 = material1->CreateBC(mat1,3,2, val1, val2);				
		//		TPZMaterial * BCond3 = material1->CreateBC(mat1,3,1, val1, val2);
		//		TPZMaterial * BCond8 = material2->CreateBC(mat2,8,1, val1, val2);			
		
//		val2(0,0)=10.0;// qx
//		val2(1,0)=0.0;// qy
//		val2(2,0)=0.0;// P
//		val2(3,0)=0.0;// S	
//		TPZMaterial * BCond5 = material1->CreateBC(mat1,5,3, val1, val2);
//		
//		val2(0,0)=0.0;// qx
//		val2(1,0)=0.0;// qy
//		val2(2,0)=0.0;// P
//		val2(3,0)=0.0;// S			
//		TPZMaterial * BCond3 = material1->CreateBC(mat1,2,3, val1, val2);		
////		TPZMaterial * BCond4 = material2->CreateBC(mat2,4,0, val1, val2);
//		
//		val2(0,0)=0.0;// qx
//		val2(1,0)=0.0;// qy
//		val2(2,0)=0.0;// P
//		val2(3,0)=0.0;// S			
//		TPZMaterial * BCond2 = material1->CreateBC(mat1,4,3, val1, val2);
////		TPZMaterial * BCond4 = material1->CreateBC(mat1,4,2, val1, val2);		
////		TPZMaterial * BCond6 = material2->CreateBC(mat2,6,3, val1, val2);
//				
//		val2(0,0)=0.0;// qx
//		val2(1,0)=0.0;// qy
//		val2(2,0)=2.0;// P
//		val2(3,0)=0.0;// S			
//		TPZMaterial * BCond4 = material1->CreateBC(mat1,3,1, val1, val2);				
////		TPZMaterial * BCond3 = material1->CreateBC(mat1,3,1, val1, val2);
////		TPZMaterial * BCond8 = material2->CreateBC(mat2,8,1, val1, val2);		
		
		mphysics->SetAllCreateFunctionsMultiphysicElem();		
		mphysics->InsertMaterialObject(BCond2);
//		mphysics->InsertMaterialObject(BCond4);
		mphysics->InsertMaterialObject(BCond3);
//		mphysics->InsertMaterialObject(BCond6);
		mphysics->InsertMaterialObject(BCond4);
//		mphysics->InsertMaterialObject(BCond8);
		mphysics->InsertMaterialObject(BCond5);		
					

		mphysics->AutoBuild();
		mphysics->AdjustBoundaryElements();
		mphysics->CleanUpUnconnectedNodes();
		
		// Creating multiphysic elements into mphysics computational mesh
		TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
		TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
		TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);	
		
		mphysics->Reference()->ResetReference();
		mphysics->LoadReferences();		
		
		// Creation of interface elements
		int nel = mphysics->ElementVec().NElements();
		for(int el = 0; el < nel; el++)
		{
			TPZCompEl * compEl = mphysics->ElementVec()[el];
			if(!compEl) continue;
			int index = compEl ->Index();
			if(compEl->Dimension() == mphysics->Dimension())
			{
				TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(mphysics->ElementVec()[index]);
				if(!InterpEl) continue;
				InterpEl->CreateInterfaces();
			}
		}
		
		return mphysics;
	}

void PosProcessBulkflux(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(0), vecnames(1);
	vecnames[0]= "FluxL2";
    
	const int dim = 2;
	int div = 5;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malhaflux.txt");
	an.Print("nothing",out);
}

void PosProcessL2(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(1), vecnames(0);
	scalnames[0]= "Solution";
    
	const int dim = 2;
	int div = 2;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malhaflux.txt");
	an.Print("nothing",out);
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(3), vecnames(1);
	
    scalnames[0] = "WeightedPressure";
    scalnames[1] = "WaterSaturation";
    scalnames[2] = "OilSaturation";
    vecnames[0] = "BulkVelocity";
    
	const int dim = 2;
	int div =2;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
    
}

void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		long n = gMesh->NElements();
		for ( long i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
		}//for i
	}//ref
}

void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv)
{
	
	TPZVec<long > subindex;
	for (long iref = 0; iref < ndiv; iref++) {
		TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
		long nel = elvec.NElements(); 
		for(long el=0; el < nel; el++){
			TPZCompEl * compEl = elvec[el];
			if(!compEl) continue;
			long ind = compEl->Index();
			compEl->Divide(ind, subindex, 0);
		}
	}
	
}

void SolveSystemTransient(REAL deltaT,REAL maxTime, TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics){
    
	//	Hard code and hard core

	//	Newton analysis and loading  initial Solution
//	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);	
//	TPZAnalysis NonLinearAn(mphysics);
	int	Nthreads = 1;
//	TPZSkylineNSymStructMatrix
	
	
//	TPZSkylineStructMatrix  matsk(mphysics);
//	NonLinearAn.SetStructuralMatrix(matsk);
//	//	an.StructMatrix()->SetNumThreads(Nthreads);
//	TPZStepSolver<STATE> step;
//	step.SetDirect(ELDLt); 					
//	NonLinearAn.SetSolver(step);	
//    //TPZFStructMatrix matsk(mphysics);
//    TPZSkylineStructMatrix matsk(mphysics);
//	NonLinearAn.SetStructuralMatrix(matsk);
//    //matsk.SetNumThreads(30);
//	TPZStepSolver<STATE> step;virtual clonedrive 5.4.6.0.rar
//	step.SetDirect(ELDLt);
//	//step.SetDirect(ELU);
//	NonLinearAn.SetSolver(step);	
	
	
	
	
//	
    TPZMaterial *mat1 = mphysics->FindMaterial(1);
//    TPZMaterial *mat2 = mphysics->FindMaterial(2);	
	
    TPZMultiphase * material1 = dynamic_cast<TPZMultiphase *>(mat1);  
//    TPZMultiphase * material2 = dynamic_cast<TPZMultiphase *>(mat2);  	
    material1->SetTimeStep(deltaT);
//    material2->SetTimeStep(deltaT);
    material1->SetTheta(1.0);
//    material2->SetTheta(0.5);

//	Starting Newton Iterations
	
	
	TPZFMatrix<STATE> DeltaX = mphysics->Solution();
	TPZFMatrix<STATE> LastSolution = mphysics->Solution();
	

	

	REAL TimeValue = 0.0;
	double Tolerance = 1.0e-5;
	int cent = 0;
	int MaxIterations = 30;
	TimeValue = cent*deltaT;
	double NormValue =1.0;
	bool StopCriteria = false;
	TPZFMatrix<STATE> RhsAtn, RhsAtnPlusOne, Residual;
	
	while (TimeValue < maxTime)
	{	
	
		material1->SetLastState();
//		material1->SetLastState();
//		material2->SetLastState();	
		NonLinearAn->Assemble();
		RhsAtn = NonLinearAn->Rhs();	
		
		material1->SetCurrentState();
		//			material2->SetCurrentState();
		NonLinearAn->Assemble();
		RhsAtnPlusOne = NonLinearAn->Rhs();
		Residual= RhsAtn + RhsAtnPlusOne;		
		NormValue = Norm(Residual);
		
//		// Getting EK
//		TPZFStructMatrix matsp(mphysics);		
//		std::set< int > materialid;
//		int matid = material1->MatId();
//		materialid.insert(matid);
//		materialid.insert(2);
//		materialid.insert(3);
//		materialid.insert(4);
//		materialid.insert(5);			
//		matsp.SetMaterialIds (materialid);
//		TPZAutoPointer<TPZGuiInterface> guiInterface;
//		TPZFMatrix<STATE> Un;
//		TPZAutoPointer <TPZMatrix<STATE> > matK2 = matsp.CreateAssemble(Un,guiInterface);
//		TPZFMatrix<STATE> TangentMat = *(matK2.operator->());
		
//#ifdef LOG4CXX
//		if(logdata->isDebugEnabled())
//		{
//			std::stringstream sout;
//			RhsAtn.Print("RhsAtn = ", sout,EMathematicaInput);				
//			RhsAtnPlusOne.Print("RhsAtnPlusOne = ", sout,EMathematicaInput);
//			NonLinearAn->Solver().Matrix()->Print("Tangent2 = ", sout,EMathematicaInput);
//			matK2->Print("Tangent = ", sout,EMathematicaInput);					
//			LastSolution.Print("LastSolution = ", sout,EMathematicaInput);
//			LOGPZ_DEBUG(logdata,sout.str())
//		}
//#endif		
		
//#ifdef LOG4CXX
//		if(logdata->isDebugEnabled())
//		{
//			std::stringstream sout;
//			RhsAtn.Print("RhsAtn = ", sout,EMathematicaInput);				
//			RhsAtnPlusOne.Print("RhsAtnPlusOne = ", sout,EMathematicaInput);
//			Residual.Print("Residual = ", sout,EMathematicaInput);
//			LOGPZ_DEBUG(logdata,sout.str())
//		}
//#endif		

//		LastSolution.Print("LastSolution = :");
//		Residual.Print("Residual = :");	
		int iterations= 0;		
		while (NormValue > Tolerance)
		{		
			
//			// Getting EK
//			TPZFStructMatrix matsp(mphysics);		
//			std::set< int > materialid;
//			int matid = material1->MatId();
//			materialid.insert(matid);
//			materialid.insert(2);
//			materialid.insert(3);
//			materialid.insert(4);
//			materialid.insert(5);			
//			matsp.SetMaterialIds (materialid);
//			TPZAutoPointer<TPZGuiInterface> guiInterface;
//			TPZFMatrix<STATE> Un;
//			TPZAutoPointer <TPZMatrix<STATE> > matK2 = matsp.CreateAssemble(Un,guiInterface);				
//			
//#ifdef LOG4CXX
//			if(logdata->isDebugEnabled())
//			{
//				std::stringstream sout;
//				RhsAtn.Print("RhsAtn = ", sout,EMathematicaInput);				
//				RhsAtnPlusOne.Print("RhsAtnPlusOne = ", sout,EMathematicaInput);
//				Residual.Print("Residual = ", sout,EMathematicaInput);
//				matK2->Print("Tangent = ", sout,EMathematicaInput);					
//				LastSolution.Print("LastSolution = ", sout,EMathematicaInput);
//				LOGPZ_DEBUG(logdata,sout.str())
//			}
//#endif


			
		
//			
//#ifdef LOG4CXX
//			if(logdata->isDebugEnabled())
//			{
//				std::stringstream sout;
//				LastSolution.Print("LastSolution = ", sout,EMathematicaInput);
//				NonLinearAn->Solution().Print("NonLinearAn.Solution() = ", sout,EMathematicaInput);			
//				LOGPZ_DEBUG(logdata,sout.str())
//			}
//#endif				
						
			
			Residual*=-1.0;
			NonLinearAn->Rhs()=Residual;
			NonLinearAn->Solve();
//			LastSolution.Print("LastSolution = :");
//			NonLinearAn->Solution().Print("NonLinearAn->Solution() = :");				
			DeltaX = NonLinearAn->Solution();
//			mphysics->LoadSolution((LastSolution + NonLinearAn->Solution()));
			NonLinearAn->LoadSolution((LastSolution + NonLinearAn->Solution()));			
			TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
			
			LastSolution = NonLinearAn->Solution();	
			
//			DeltaX.Print("DeltaX = :");
//			NonLinearAn->Solution().Print("NonLinearAn->Solution() = :");		

			
			material1->SetCurrentState();
//			TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);			
			NonLinearAn->Assemble();
			RhsAtnPlusOne = NonLinearAn->Rhs();
			Residual= RhsAtn + RhsAtnPlusOne;		
			NormValue = Norm(Residual);
			
//			Residual.Print("Residual = :");
//			mphysics->Solution().Print("NonLinearAn->Solution() = :");				
			
			
//#ifdef LOG4CXX
//			if(logdata->isDebugEnabled())
//			{
//				std::stringstream sout;				
//				NonLinearAn.Solution().Print("LastSolution = ", sout,EMathematicaInput);
//				Residual.Print("Residual = ", sout,EMathematicaInput);
//				LOGPZ_DEBUG(logdata,sout.str())
//			}
//#endif					
			iterations++;
			std::cout << " Iteration number = : " << iterations  << "\n Corresponding L2 norm = : " << NormValue <<  std::endl;
			if (iterations == MaxIterations) 
			{
				StopCriteria = true;
				std::cout << " Time Step number = : " << iterations  << "\n Exceed max iterations numbers = : " << MaxIterations <<  std::endl;					
				break;
			}
		}	
		
//	NonLinearAn.LoadSolution(LastSolution);			
    std::string outputfile;
	outputfile = "TransientSolution";
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
	PosProcessMultphysics(meshvec,mphysics,*NonLinearAn,plotfile);		
		
	if (StopCriteria) {
		std::cout << " Iteration number = : " << iterations  << "\n Corresponding L2 norm = : " << NormValue <<  std::endl;		
		break;
	}
		
	cent++;
	TimeValue = cent*deltaT;
		
	std::cout << " Time Step number = : " << cent  << "\n Corresponding time = : " << TimeValue <<  std::endl;	
    
	}
	
	CheckConvergence(RhsAtn,NonLinearAn, meshvec, mphysics);	
	
	
}


// Setting up initial conditions

void SolveSyst(TPZAnalysis &an, TPZCompMesh *Cmesh)
{

	TPZSkylineStructMatrix full(Cmesh);
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt);
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
//	//Saida de Dados: solucao e  grafico no VT
//	ofstream file("Solutout");
//	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

TPZCompMesh *L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
	int dim = 2;
	TPZL2Projection *material;
	material = new TPZL2Projection(1, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialPressure);
    material->SetForcingFunction(forcef);
	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	cmesh->AutoBuild();
	return cmesh;
	
}

TPZCompMesh *L2ProjectionQ(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
	int dim = 2;
	TPZL2Projection *material;
	material = new TPZL2Projection(1, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialFlux);
    material->SetForcingFunction(forcef);
	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	cmesh->AutoBuild();
	
	return cmesh;
	
}

// It requires modfify L2 number os state variables
void InitialFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
	double x = pt[0];
	double y = pt[1];
    disp[0] = 0.0;
//    disp[1] = 0.0;	
    
}

void InitialPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
	double x = pt[0];
	double y = pt[1];
    disp[0] = 0.0*(10.0 - 10.0 * x);
    
}

void InitialSaturation(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
	double x = pt[0];
	double y = pt[1];
    disp[0] = 0.0;
    
}


void CheckConvergence(TPZFMatrix<STATE> SoltUattn,TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics)
{

	TPZFMatrix<REAL> U = NonLinearAn->Solution();
	long neq = mphysics->NEquations();
	int nsteps = 10;
	REAL alpha;
	TPZFMatrix<REAL> alphas(nsteps,1,0.0),ResNorm(nsteps,1,0.0),ConvergenceOrder(nsteps-1,1,0.0);

	TPZFMatrix<REAL> DeltaX(neq,1,0.0001),ResAlpha(neq,0.0);

    for(int i = 0; i < nsteps; i++)
    {
        alpha = (i+1.0)/100.0;
        alphas(i,0) = log(alpha);
		NonLinearAn->LoadSolution(U);
		TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
		ComputeResidual(SoltUattn,alpha, DeltaX, ResAlpha, NonLinearAn, meshvec, mphysics);
		
		ResNorm(i,0)=log(Norm(ResAlpha));
		
	}
	
    for(int i = 0; i < nsteps - 1; i++){ ConvergenceOrder(i,0) =  (ResNorm(i+1,0)-ResNorm(i,0))/(alphas(i+1,0)-alphas(i,0));}	
	
	ConvergenceOrder.Print("ConvergenceOrder = ");	


}

void ComputeResidual(TPZFMatrix<STATE> SoltUattn, REAL alpha, TPZFMatrix<STATE> DeltaU, TPZFMatrix<STATE> &ResAlpha, TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics)
{

	TPZFMatrix<STATE> TangentRes;	
    TPZMaterial *mat1 = mphysics->FindMaterial(1);
    TPZMultiphase * material1 = dynamic_cast<TPZMultiphase *>(mat1); 	
	
	
	// Computing the first part of the residual expresion.		
	material1->SetCurrentState();
	NonLinearAn->Assemble();
		
	TPZFMatrix<STATE> RhsAtnPlusOne = NonLinearAn->Rhs();
	TPZFMatrix<STATE> ResidualAtU = SoltUattn + RhsAtnPlusOne;		
	NonLinearAn->Solver().Matrix()->Multiply(alpha*DeltaU,TangentRes);
		
	NonLinearAn->LoadSolution(NonLinearAn->Solution()+alpha*DeltaU);			
	TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);		
	
	material1->SetCurrentState();
	NonLinearAn->Assemble();
	RhsAtnPlusOne = NonLinearAn->Rhs();
	TPZFMatrix<STATE> ResidualAtUplusAlphaDeltaU = SoltUattn + RhsAtnPlusOne;
	
	ResAlpha = (ResidualAtUplusAlphaDeltaU - (ResidualAtU + TangentRes));

}