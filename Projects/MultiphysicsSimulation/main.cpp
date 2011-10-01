#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "tpzgeolinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

#include "pzmultiphysicselement.h"

#include "pzlog.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
#endif

using namespace std;

const int matId = 1;
const int dirichlet = 0;
const int bc0 = -1;

TPZGeoMesh *MalhaGeom( );

TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh,int pOrder);
void RefinamentoUniforme(TPZGeoMesh  *gMesh, int nh);
void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh, int MatId, int indexEl);
void RefinElemComp(TPZCompMesh  *cMesh, int indexEl);
void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file);
void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file);
void GeoElMultiphysicVec(TPZManVector<TPZCompMesh  *> cmeshVec,std::set <int> &geoelVec);
void AddElements(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh);


int main(int argc, char *argv[])
{	
	std::string logs("log4cxx.doubleprojection1d");
	InitializePZLOG();
	
	int p =2;
	//primeira malha
	
      TPZGeoMesh * gmesh = MalhaGeom();
	ofstream arg1("gmesh1.txt");
	gmesh->Print(arg1);
	ofstream file1("malhageoInicial.vtk");
	PrintGMeshVTK(gmesh, file1);
	
	TPZCompMesh * cmesh1= MalhaComp(gmesh,  p);
	ofstream arg2("cmesh1.txt");
	cmesh1->Print(arg2);
	TPZCompMesh * cmesh2 = MalhaComp(gmesh, p);
	ofstream arg3("cmesh2.txt");
	cmesh2->Print(arg3);
	
	//Refinar as malhas
	gmesh->ResetReference();
	cmesh1->LoadReferences();
	RefinElemComp(cmesh1, 7);
	RefinElemComp(cmesh1, 10);
	cmesh1->AdjustBoundaryElements();
	cmesh1->CleanUpUnconnectedNodes();
	
	ofstream arg4("cmesh12.txt");
	cmesh1->Print(arg4);
	ofstream arg5("gmesh2.txt");
	gmesh->Print(arg5);
	ofstream file3("malhageo1.vtk");
	PrintGMeshVTK(gmesh, file3);
		
	
	gmesh->ResetReference();
	cmesh2->LoadReferences();
	RefinElemComp(cmesh2, 6);
	RefinElemComp(cmesh2, 7);
	cmesh2->AdjustBoundaryElements();
	cmesh2->CleanUpUnconnectedNodes();
	
	TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh1;
	meshvec[1] = cmesh2;
	
	ofstream arg6("cmesh22.txt");
	cmesh2->Print(arg6);
	ofstream arg7("gmesh3.txt");
	gmesh->Print(arg7);
	ofstream file5("malhageo2.vtk");
	PrintGMeshVTK(gmesh, file5);
	
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    mphysics->MaterialVec() = cmesh1->MaterialVec();
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->AutoBuild();
	
	AddElements(meshvec, mphysics);
	
#ifdef LOG4CXX
    {
        std::stringstream sout;
        mphysics->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	
	//-------------------------------------------------
	//int ngeo=gmesh->NElements();
//	TPZStack <TPZGeoEl*> geovec;
//	cout<<"Num element Geom = "<<ngeo<<endl;
//	int ind=0;
//	for(int iel=0; iel<ngeo; iel++){
//		TPZGeoEl * gEl = gmesh->ElementVec()[iel];
//		TPZStack<TPZCompElSide> ElSideVec;
//		int ns = gEl->NSides();
//		TPZGeoElSide *geoside = new TPZGeoElSide(gEl,ns-1);
//		//geoside->SetElement(gEl);
//		geoside->HigherLevelCompElementList2(ElSideVec, 1,1);
//		int nel = ElSideVec.NElements();
//		if (nel==0) {
//			//std::cout<<" Num Elem Comp higher level = "<< nel<<endl;
//			geovec.Push(gEl);
//			//cout << " ======= ======= "<<endl;
//			geovec[ind]->Print();
//			cout << "======= ====== "<<endl;
//			ind++;
//		}
//	}
//	cout<<"ind = "<<ind<<endl;
	
	
	std::set<int> geoelVec;
	TPZManVector<TPZCompMesh *,2> cmeshVec(2);
	cmeshVec[0]=cmesh1;
	cmeshVec[1]=cmesh2;
	GeoElMultiphysicVec(cmeshVec, geoelVec);
	
	set<int>::iterator it;
	cout << "myset contains:";
	for (it=geoelVec.begin() ; it != geoelVec.end(); it++ )
		cout << " " << *it;
	cout << endl;
	
	return EXIT_SUCCESS;
}


TPZGeoMesh *MalhaGeom()
{
	
	int Qnodes = 6;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolLine(2);
	
	//indice dos nos
	int id = 0;
	REAL valx;
	REAL dx=0.5;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = 1. - xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,1. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
	id++;
	
	TopolLine[0] = 1;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
	id++;
	
	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
	id++;
	
	TopolLine[0] = 3;
	TopolLine[1] = 4;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
	id++;
	
	TopolLine[0] = 4;
	TopolLine[1] = 5;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
	id++;
	
	TopolLine[0] = 5;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
	id++;
	
	TopolQuad[0] = 0;
	TopolQuad[1] = 1;
	TopolQuad[2] = 4;
	TopolQuad[3] = 5;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
	id++;
	
	TopolQuad[0] = 1;
	TopolQuad[1] = 2;
	TopolQuad[2] = 3;
	TopolQuad[3] = 4;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
	
	gmesh->BuildConnectivity();
			
	//ofstream arg("gmesh.txt");
//	gmesh->Print(arg);
		
	return gmesh;
	
}

TPZCompMesh*MalhaComp(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	TPZAutoPointer<TPZMaterial> mat(material);
	
	REAL diff = 1.;
	REAL conv = 0.;
	TPZVec<REAL> convdir(3,0.);
	REAL flux = 0.;
	
	material->SetParameters(diff, conv, convdir);
	material->SetInternalFlux( flux);
		
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
	

	///Inserir condicao de contorno
	TPZFMatrix val1(1,1,0.), val2(1,1,0.);
	TPZAutoPointer<TPZMaterial> BCond = material->CreateBC(mat, bc0,0, val1, val2);
	cmesh->InsertMaterialObject(BCond);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	//cmesh->AdjustBoundaryElements(); 
	//cmesh->CleanUpUnconnectedNodes();
	
	//ofstream arg("cmesh.txt");
//	cmesh->Print(arg);
	
	return cmesh;
}

void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh){
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gMesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2) gel->Divide (filhos);
		}//for i
	}//ref
}

void RefinamentoUniforme(TPZGeoMesh * gMesh, int nh, int MatId, int indexEl){
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gMesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2){
				if (gel->MaterialId()== MatId && gel-> Index()==indexEl){
					gel->Divide (filhos);
				}
			}
		}//for i
	}//ref
}

void RefinElemComp(TPZCompMesh  *cMesh, int indexEl){
	
	TPZVec<int > subindex; 
	int nel = cMesh->ElementVec().NElements(); 
	for(int el=0; el < nel; el++){
		TPZCompEl * compEl = cMesh->ElementVec()[el];
		if(!compEl) continue;
		int ind = compEl->Index();
		if(ind==indexEl){
			compEl->Divide(indexEl, subindex, 1);
		}
	}	
}

void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file)
{
	file.clear();
	int nelements = gmesh->NElements();
	
	std::stringstream node, connectivity, type;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int actualNode = -1, size = 0, nVALIDelements = 0;
	
	for(int el = 0; el < nelements; el++)
	{        
		if(gmesh->ElementVec()[el]->Type() == EPoint)//Exclude Lines and Arc3D
		{
			continue;
		}
		if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines and Arc3D
		{
			continue;
		}
		if(gmesh->ElementVec()[el]->HasSubElement())
		{
			continue;
		}
		
		int elNnodes = gmesh->ElementVec()[el]->NNodes();
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(int t = 0; t < elNnodes; t++)
		{
			for(int c = 0; c < 3; c++)
			{
				double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
				node << coord << " ";
			}            
			node << std::endl;
			
			actualNode++;
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType = -1;
		switch (gmesh->ElementVec()[el]->Type())
		{
			case (ETriangle):
			{
				elType = 5;
				break;                
			}
			case (EQuadrilateral ):
			{
				elType = 9;
				break;                
			}
			case (ETetraedro):
			{
				elType = 10;
				break;                
			}
			case (EPiramide):
			{
				elType = 14;
				break;                
			}
			case (EPrisma):
			{
				elType = 13;
				break;                
			}
			case (ECube):
			{
				elType = 12;
				break;                
			}
			default:
			{
				//ElementType NOT Found!!!
				DebugStop();
				break;    
			}
		}
		
		type << elType << std::endl;
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();
	
	file << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str();
	
	file.close();
}


void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file)
{
	TPZGeoMesh *gmesh;
	
	//    RefPatternMesh();
	//TPZGeoMesh * gmesh = refp->Mesh();
	PrintGMeshVTK(gmesh, file);
}

void GeoElMultiphysicVec(TPZManVector<TPZCompMesh *> cmeshVec, std::set<int> &geoelVec){
	
	if(cmeshVec.NElements() == 0) return;
	TPZCompMesh *cmesh = cmeshVec[0];
	TPZGeoMesh *gmesh = cmesh->Reference();
	gmesh->ResetReference();
	int isub;
	int ncm = cmeshVec.NElements();
	for (isub=0; isub<ncm; isub++) {
		cmeshVec[isub]->LoadReferences();
	}
	int ncel;
	TPZStack<TPZCompElSide> sidevec;
	for(int i = 0; i< ncm; i++){
		ncel = cmeshVec[i]->NElements();
		for (int j=0; j<ncel; j++) {
			TPZCompEl * cel = cmeshVec[i]->ElementVec()[j];
			if(cel){
				TPZGeoEl *geoel = cel->Reference();
				if (!geoel) {
					std::cout << "Geoel nulo!\n";
					DebugStop();
				}
				int ns = geoel->NSides();
				TPZGeoElSide *geoside = new TPZGeoElSide(geoel,ns-1);
				sidevec.Resize(0);
				geoside->HigherLevelCompElementList2(sidevec, 1,1);
				int nel = sidevec.NElements();
				if (nel==0){
					//std::cout << "Incluindo elemento " << geoel->Index() << std::endl;
					geoelVec.insert(geoel->Index());
				}
			}
		}
	}
		//cout<<"num Elemento geom= "<<geoelVec.size()<<endl;
	//set<int>::iterator it;
//	cout << "myset contains:";
//	for (it=geoelVec.begin() ; it != geoelVec.end(); it++ )
//		cout << " " << *it;
//	cout << endl;
	
}

void AddElements(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh)
{
	TPZGeoMesh *gmesh = MFMesh->Reference();
	gmesh->ResetReference();
	int nMFEl = MFMesh->NElements();
	int nmesh = cmeshVec.size();
	int imesh;
	for(imesh = 0; imesh<nmesh; imesh++)
	{
		cmeshVec[imesh]->LoadReferences();
		int iel, is;
		for(iel=0; iel<nMFEl; iel++)
		{
			TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *> (MFMesh->ElementVec()[iel]);
			if(!mfcel)
			{
				DebugStop();
			}
			TPZGeoEl *gel = mfcel->Reference();
			TPZStack<TPZCompElSide> celstack;
			TPZGeoElSide gelside(gel,gel->NSides()-1);
			if (gel->Reference()) {
				celstack.Push(gelside.Reference());
			}
			gelside.ConnectedCompElementList(celstack, 0, 0);
			for(is=0;is<celstack.size();is++) {
				TPZGeoElSide gelside = celstack[is].Reference();
				
				if(gelside.Element()->Dimension()!=gel->Dimension()) {
					if(is!=celstack.size()-1) {
						celstack[is] = celstack.Pop();
						is--;
					}
					else {
						celstack.Pop();
					}

				}
			}
			if(celstack.size() != 1)
			{
				DebugStop();
			}
			
			mfcel->AddElement(celstack[0].Element(), imesh);
		}
		gmesh->ResetReference();
	}
		
	
}
