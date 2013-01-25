/**
 * @file
 * @brief Contains the implementation of the Multiphysic interface methods
 * @author Agnaldo
 * @since 10/26/11.
 */

#include "TPZMultiphysicsInterfaceEl.h"
#include "pzelmat.h"
#include "pzinterpolationspace.h"
#include "pzmaterial.h"
#include "pzmultiphysicselement.h"
#include "tpzintpoints.h"
#include "pzdiscgal.h"
#include "pzmultiphysicscompel.h"
#include "pzgeoel.h"

#include "pzgraphel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "tpzgraphelt2dmapped.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel.h"



TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement() : TPZCompEl(),fLeftElSide(0), fRightElSide(0)
{
}

TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int &index,
                                                                    TPZCompElSide leftside, TPZCompElSide rightside) : TPZCompEl(mesh, ref, index),fLeftElSide(leftside), fRightElSide(rightside)
{
	
//	ref->SetReference(this);
//	ref->IncrementNumInterfaces();
//	
//	if (fLeftElSide.Side() == -1 || fRightElSide.Side() == -1){
//		PZError << "Error at " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " Side should not be -1\n";
//		DebugStop();
//	}
//	
//	this->SetLeftRightElement(fLeftElSide, fRightElSide);
//	
//	this->IncrementElConnected();
	
}

void TPZMultiphysicsInterfaceElement::IncrementElConnected(){
	const int ncon = this->NConnects();
	for(int i = 0; i < ncon; i++){
		int index = this->ConnectIndex(i);
		fMesh->ConnectVec()[index].IncrementElConnected();
	}
}

/** @brief create a copy of the given element */
TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy) : TPZCompEl(mesh,copy)
{
    TPZCompElSide left = copy.Left();
    TPZCompElSide right = copy.Right();
    if (!left || !right ) {
        DebugStop();
    }
    int leftside = left.Side();
    int rightside = right.Side();
    int leftindex = left.Element()->Index();
    int rightindex = right.Element()->Index();
    TPZCompEl *leftel = mesh.ElementVec()[leftindex];
    TPZCompEl *rightel = mesh.ElementVec()[rightindex];
    if (!leftel || !rightel) {
        DebugStop();
    }
    fLeftElSide = TPZCompElSide(leftel,leftside);
    fRightElSide = TPZCompElSide(rightel,rightside);
}

/** @brief create a copy of the given element using index mapping */
TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy, std::map<int,int> & gl2lcConMap,
                                std::map<int,int> & gl2lcElMap)
{
    TPZCompElSide left = copy.Left();
    TPZCompElSide right = copy.Right();
    if (!left || !right ) {
        DebugStop();
    }
    int leftside = left.Side();
    int rightside = right.Side();
    int leftindex = left.Element()->Index();
    int rightindex = right.Element()->Index();
    if (gl2lcElMap.find(leftindex) == gl2lcElMap.end() || gl2lcElMap.find(rightindex) == gl2lcElMap.end()) {
        DebugStop();
    }
    TPZCompEl *leftel = mesh.ElementVec()[gl2lcElMap[leftindex]];
    TPZCompEl *rightel = mesh.ElementVec()[gl2lcElMap[rightindex]];
    if (!leftel || !rightel) {
        DebugStop();
    }
    fLeftElSide = TPZCompElSide(leftel,leftside);
    fRightElSide = TPZCompElSide(rightel,rightside);
    
}




TPZMultiphysicsInterfaceElement::~TPZMultiphysicsInterfaceElement(){
}

void TPZMultiphysicsInterfaceElement::ComputeSideTransform(TPZManVector<TPZCompElSide> &Neighbor, TPZManVector<TPZTransform> &transf)
{
    TPZGeoEl *gel = Reference();
    int side = gel->NSides()-1;
    TPZGeoElSide thisside(gel,side);
    int numneigh = Neighbor.size();
    for (int in=0; in<numneigh; in++) {
        TPZGeoElSide gelside = Neighbor[in].Reference();
        if(! thisside.NeighbourExists(gelside))
        {
            DebugStop();
        }
        TPZTransform tr(thisside.Dimension());
        thisside.SideTransform3(gelside, tr);
        transf[in] = tr;
    }
}//ComputeSideTransform

/**
 * Add elements to the list of left and right elements
 */
void TPZMultiphysicsInterfaceElement::SetLeftRightElement(TPZCompElSide &leftel, TPZCompElSide &rightel)
{
    fLeftElSide = leftel;
    fRightElSide = rightel;
}

/**
 * Get left and right elements
 */
void TPZMultiphysicsInterfaceElement::GetLeftRightElement(TPZCompElSide &leftel, TPZCompElSide &rightel)
{
	leftel = fLeftElSide;
	rightel = fRightElSide;
}

/** @brief Returns the number of nodes of the element */
int TPZMultiphysicsInterfaceElement::NConnects() const
{
    return fLeftElSide.Element()->NConnects() + fRightElSide.Element()->NConnects();
}

/**
 * @brief Returns the index of the ith connectivity of the element
 * @param i connectivity index who want knows
 */
int TPZMultiphysicsInterfaceElement::ConnectIndex(int i) const
{
#ifdef DEBUG
    if (i < 0) {
        DebugStop();
    }
#endif
    int nleft = fLeftElSide.Element()->NConnects();
    if (i < nleft) {
        return fLeftElSide.Element()->ConnectIndex(i);
    }
    int nright = fRightElSide.Element()->NConnects();
    if (i < nleft+nright) {
        return fRightElSide.Element()->ConnectIndex(i-nleft);
    }
    DebugStop();
    return -1;
}



void TPZMultiphysicsInterfaceElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
	TPZDiscontinuousGalerkin  * material = dynamic_cast<TPZDiscontinuousGalerkin *> (Material());
	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
		ek.Reset();
		ef.Reset();
		return;
	}
	
	InitializeElementMatrix(ek,ef);
	
	if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
    TPZGeoEl *leftgel = leftel->Reference();
    TPZGeoEl *rightgel = rightel->Reference();
#ifdef DEBUG
    if (!leftel || !rightel) {
        DebugStop();
    }
#endif
    TPZManVector<TPZMaterialData,6> datavecleft,datavecright;
    TPZMaterialData data;
    InitMaterialData(datavecleft, leftel);
    InitMaterialData(datavecright, rightel);
    TPZManVector<TPZTransform> leftcomptr, rightcomptr;
    leftel->AffineTransform(leftcomptr);
    rightel->AffineTransform(rightcomptr);
    InitMaterialData(data);
    int intleftorder = leftel->IntegrationOrder();
    int intrightorder = rightel->IntegrationOrder();
    int integrationorder = MAX(intleftorder, intrightorder);
    TPZGeoEl *gel = Reference();
    int dimension = gel->Dimension();
    int thisside = gel->NSides()-1;
    TPZFNMatrix<9,REAL> jac(dimension,dimension),axes(dimension,3), jacInv(dimension,dimension);
    
    TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(thisside, integrationorder);
    TPZManVector<REAL,3> Point(dimension), leftPoint(leftel->Dimension()), rightPoint(rightel->Dimension());
    TPZGeoElSide neighleft(fLeftElSide.Reference()), neighright(fRightElSide.Reference());
    TPZTransform trleft(dimension),trright(dimension);
    TPZGeoElSide gelside(this->Reference(),thisside);
    // compute the transformation between neighbours
    gelside.SideTransform3(neighleft, trleft);
    gelside.SideTransform3(neighright, trright);
    TPZTransform leftloctr = leftgel->SideToSideTransform(neighleft.Side(), leftgel->NSides()-1);
    TPZTransform rightloctr = rightgel->SideToSideTransform(neighright.Side(), rightgel->NSides()-1);
    // transform from the element to the interior of the neighbours
    trleft = leftloctr.Multiply(trleft);
    trright = rightloctr.Multiply(trright);
    
    
    int nintpoints = intrule->NPoints();
    for (int ip =0; ip<nintpoints; ip++) {
        REAL weight;
        intrule->Point(ip, Point, weight);
        ComputeRequiredData(Point, data);
        weight *= fabs(data.detjac);
        trleft.Apply(Point, leftPoint);
        leftel->ComputeRequiredData(leftPoint, leftcomptr, datavecleft);
        trright.Apply(Point, rightPoint);
        rightel->ComputeRequiredData(rightPoint, rightcomptr, datavecright);
        material->ContributeInterface(data , datavecleft, datavecright, weight, ek.fMat, ef.fMat);
    }	
	
}//CalcStiff

void TPZMultiphysicsInterfaceElement::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
	const int ncon = this->NConnects();
	int numeq = 0;
	int ic;
	
	for(ic=0; ic<ncon; ic++)
	{
        TPZConnect &c = Connect(ic);
		numeq += c.NShape()*c.NState();
	}
	
    TPZMultiphysicsElement *mfcel_left = dynamic_cast<TPZMultiphysicsElement *>(fLeftElSide.Element());
    TPZMultiphysicsElement *mfcel_right = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
    if (! mfcel_left || !mfcel_right) {
        DebugStop();
    }
	//nstate=1;
    int numloadcases = 1;
    TPZMultiphysicsElement *msp  = dynamic_cast <TPZMultiphysicsElement *>(fLeftElSide.Element());
    if (!msp) {
        DebugStop();
    }
    TPZMaterial *mat = msp->Material();
    int nstate = mat->NStateVariables();
    numloadcases = mat->NumLoadCases();        
	
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,numloadcases);
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fNumStateVars = nstate;
	ef.fNumStateVars = nstate;
	
	int i;
	for(i=0; i<ncon; i++)
    {
        TPZConnect &c = Connect(i);
        int ndof = Connect(i).NShape()*c.NState();
#ifdef DEBUG
        if (c.NDof(*Mesh()) != ndof) {
            DebugStop();
        }
#endif
		ek.fBlock.Set(i,ndof);
		ef.fBlock.Set(i,ndof);
	}
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ek.fConnect)[i] = ConnectIndex(i);
		(ef.fConnect)[i] = ConnectIndex(i);
	}
	
}//void


/** @brief Initialize the material data for the neighbouring element */
void TPZMultiphysicsInterfaceElement::InitMaterialData(TPZVec<TPZMaterialData> &data, TPZMultiphysicsElement *mfcel)
{
	data.resize(mfcel->NMeshes());
	mfcel->InitMaterialData(data);
}

/** @brief initialize the material data for the geometric data */
void TPZMultiphysicsInterfaceElement::InitMaterialData(TPZMaterialData &data)
{
    TPZGeoEl *gel = Reference();
    int dim = gel->Dimension();
    int nsides = gel->NSides();
    TPZManVector<REAL> center(dim);
    gel->CenterPoint(nsides-1 , center);
    TPZGeoElSide gelside(gel,nsides-1);
    gelside.Normal(center, fLeftElSide.Element()->Reference(), fRightElSide.Element()->Reference(), data.normal);
    data.axes.Redim(dim,3);
    data.jacobian.Redim(dim,dim);
	data.jacinv.Redim(dim,dim);
	data.x.Resize(3);
    
}

/** @brief Compute the data needed to compute the stiffness matrix at the integration point */
void TPZMultiphysicsInterfaceElement::ComputeRequiredData(TPZVec<REAL> &point, TPZMaterialData &data)
{
    data.intGlobPtIndex = -1;
    TPZGeoEl *gel = Reference();
    TPZGeoElSide gelside(gel,gel->NSides()-1);
    gel->Jacobian(point, data.jacobian, data.axes, data.detjac, data.jacinv);
    //ComputeRequiredData(Point,data);
    gelside.Normal(point, fLeftElSide.Element()->Reference(), fRightElSide.Element()->Reference(), data.normal);

}

/** @brief Compute the required data from the neighbouring elements */
void TPZMultiphysicsInterfaceElement::ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform> &trvec, TPZMultiphysicsElement *Neighbour, TPZVec<TPZMaterialData> &data)
{
    DebugStop();
    //Neighbour->ComputeR
}

void TPZMultiphysicsInterfaceElement::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension)
{
	TPZGeoEl *ref = Reference();
	if (ref->Dimension() != dimension) {
		return;
	}
	
	TPZMaterial * material = Material();
	int mat = material->Id();
	int nsides = ref->NSides();
	
	if(dimension == 2 && mat > 0){
		if(nsides == 9){
			new TPZGraphElQ2dd(this,&grmesh);
			return;
		}
		if(nsides == 7){
			new TPZGraphElT2dMapped(this,&grmesh);
			return;
		}
	}//2d
	
	if(dimension == 3 && mat > 0){
		if(nsides == 27){
			new TPZGraphElQ3dd(this,&grmesh);
			return;
		}//cube
		if(nsides == 21){
			new TPZGraphElPrismMapped(this,&grmesh);
			return;
		}//prism
		if(nsides == 15){
			new TPZGraphElT3d(this,&grmesh);
			return;
		}//tetra
		if(nsides == 19){
			new TPZGraphElPyramidMapped(this,&grmesh);
			return;
		}//pyram
	}//3d
	
	if(dimension == 1 && mat > 0){
		new TPZGraphEl1dd(this,&grmesh);
	}//1d
	
}

void TPZMultiphysicsInterfaceElement::Solution(TPZVec<REAL> &qsi, int var,TPZVec<REAL> &sol)
{
	
	if(var >= 100) {
		TPZCompEl::Solution(qsi,var,sol);
		return;
	}
	
	TPZMaterial * material = this->Material();
	if(!material){
		sol.Resize(0);
		return;
	}
	
	if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
	
	TPZCompElSide LeftSide;
	TPZCompElSide RightSide;
	this->GetLeftRightElement(LeftSide, RightSide);
	
    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (LeftSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(RightSide.Element());
	
    TPZGeoEl *leftgel = leftel->Reference();
    TPZGeoEl *rightgel = rightel->Reference();	
	
    TPZManVector<TPZMaterialData,6> datavecleft,datavecright;
    TPZMaterialData data;
    InitMaterialData(datavecleft, leftel);
    InitMaterialData(datavecright, rightel);

	int nrefl = datavecleft.size();
	for(int i = 0; i<nrefl; i++)
		
	{
		datavecleft[i].SetAllRequirements(false);
        datavecleft[i].fNeedsSol = true;
		datavecleft[i].fNeedsNeighborSol = false;
		datavecleft[i].fNeedsNeighborCenter = false;
		datavecleft[i].fNeedsNormal = false;
	}
	
	int nrefr = datavecright.size();
	for(int i = 0; i<nrefr; i++)
		
	{
		datavecright[i].SetAllRequirements(false);
        datavecright[i].fNeedsSol = true;
		datavecright[i].fNeedsNeighborSol = false;
		datavecright[i].fNeedsNeighborCenter = false;
		datavecright[i].fNeedsNormal = false;
	}
	
    TPZManVector<TPZTransform> leftcomptr, rightcomptr;
    leftel->AffineTransform(leftcomptr);
    rightel->AffineTransform(rightcomptr);
    InitMaterialData(data);	
	TPZTransform lefttr;
	TPZTransform righttr;	
	
	//		Integration points in left and right elements: making transformations to neighbour elements
	this->ComputeSideTransform(LeftSide, lefttr);
	this->ComputeSideTransform(RightSide, righttr);	
	
	TPZVec<REAL> myqsi;
	myqsi.resize(qsi.size());	
	lefttr.Apply(qsi, myqsi);
	lefttr.Apply(qsi, myqsi);
	
	leftel->ComputeRequiredData(myqsi, leftcomptr, datavecleft);
	rightel->ComputeRequiredData(myqsi, rightcomptr, datavecright);
		
	material->Solution(data,datavecleft,datavecright,var, sol,LeftSide.Element(),RightSide.Element());
}

void TPZMultiphysicsInterfaceElement::ComputeSideTransform(TPZCompElSide &Neighbor, TPZTransform &transf){
	TPZGeoEl * neighel = Neighbor.Element()->Reference();
	const int dim = this->Dimension();
	TPZTransform LocalTransf(dim);
	TPZGeoElSide thisgeoside(this->Reference(), this->Reference()->NSides()-1);
	TPZGeoElSide neighgeoside(neighel, Neighbor.Side());
	thisgeoside.SideTransform3(neighgeoside, LocalTransf);
	
	TPZGeoElSide highdim(neighel, neighel->NSides()-1);
	transf = neighgeoside.SideToSideTransform(highdim).Multiply(LocalTransf);
}//ComputeSideTransform


