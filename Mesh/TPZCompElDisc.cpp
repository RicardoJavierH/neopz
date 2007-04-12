//$Id: TPZCompElDisc.cpp,v 1.88 2007-04-12 20:04:49 tiago Exp $

// -*- c++ -*-

#include "pztransfer.h"
#include "pzelmat.h"
//#include "pzelgc3d.h"
//#include "pzelgt3d.h"
//#include "pzelgpi3d.h"
//#include "pzelgpr3d.h"
//#include "pzelgt2d.h"
//#include "pzelgq2d.h"
#include "pzmatrix.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzerror.h"
#include "pzconnect.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pztempmat.h"
#include "pzmanvector.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "TPZConservationLaw.h"
#include "pzgraphel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "pzgraphel.h"
#include "pzmeshid.h"
#include <sstream>
//#include "TPZFlowCMesh.h"

#include "time.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include <math.h>
#include <stdio.h>
#include "pzmaterialdata.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcompeldisc"));
#endif


using namespace pzshape;
using namespace std;

TPZCompElDisc::TPZCompElDisc() : TPZInterpolationSpace(), fCenterPoint(3,0.)
{
  fShapefunctionType = pzshape::TPZShapeDisc::ETensorial;
}
//construtor do elemento aglomerado
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,int &index) :
		TPZInterpolationSpace(mesh,0,index), fCenterPoint(3) {
  fShapefunctionType = pzshape::TPZShapeDisc::EOrdemTotal;
}

TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy) :
                             TPZInterpolationSpace(mesh,copy), fCenterPoint(copy.fCenterPoint) {
  fShapefunctionType = copy.fShapefunctionType;
//  fReference = copy.fReference;
  TPZAutoPointer<TPZMaterial> mat = copy.Material();
}


TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,
                             const TPZCompElDisc &copy,
                             std::map<int,int> &gl2lcConMap,
                             std::map<int,int> &gl2lcElMap) : TPZInterpolationSpace(mesh,copy),
                                                              fCenterPoint(copy.fCenterPoint)
{
  fShapefunctionType = copy.fShapefunctionType;
  TPZAutoPointer<TPZMaterial> mat = copy.Material();
  gl2lcElMap[copy.fIndex] = this->fIndex;
}

TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy,int &index) :
		TPZInterpolationSpace(mesh,copy,index), fCenterPoint(copy.fCenterPoint) {
  fShapefunctionType = copy.fShapefunctionType;
//  fReference = copy.fReference;
  //criando nova malha computacional
  Reference()->SetReference(this);
  TPZAutoPointer<TPZMaterial> mat = copy.Material();
  fConstC = copy.fConstC;
  CreateMidSideConnect();
  this->SetDegree( TPZCompEl::gOrder );
  //as interfaces foram clonadas
}

//construtor do elemento descont�uo
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoEl *ref,int &index) :
		TPZInterpolationSpace(mesh,ref,index), fCenterPoint(3) {
  fShapefunctionType = pzshape::TPZShapeDisc::EOrdemTotal;
  switch(ref->Type()) {
    case EQuadrilateral:
    case ECube:
    case EPrisma:
      fShapefunctionType =
#ifdef _AUTODIFF
      pzshape::TPZShapeDisc::EOrdemTotal;
#else
      pzshape::TPZShapeDisc::ETensorial;
#endif
      break;
    default:
      break;
  }
//  fReference = ref;
  ref->SetReference(this);
  //fMesh = &mesh;
  CreateMidSideConnect();
  this->SetDegree( TPZCompEl::gOrder );
  ref->CenterPoint(ref->NSides()-1,fCenterPoint);
  TPZVec<REAL> csi(fCenterPoint);
  ref->X(csi,fCenterPoint);
  fConstC = NormalizeConst();
  //criando os elementos interface
  CreateInterfaces();
}

REAL TPZCompElDisc::NormalizeConst()
{
  TPZGeoEl *ref = Reference();
  //maior distancia entre o ponto interior e os v�tices do elemento
  int nnodes = ref->NNodes(),i;
  if(nnodes == 1) return 1.0;//elemento ponto
  REAL maxdist,dist;
  int inode = ref->NodeIndex(0);//primeiro n�do elemento
  TPZGeoNode node = ref->Mesh()->NodeVec()[inode];
  maxdist = pow(node.Coord(0)-fCenterPoint[0],2.)+pow(node.Coord(1)-fCenterPoint[1],2.);
  maxdist += pow(node.Coord(2)-fCenterPoint[2],2.);
  maxdist = sqrt(maxdist);
  for(i=1;i<nnodes;i++){
    inode = ref->NodeIndex(i);//n� sub-seguintes
    node = ref->Mesh()->NodeVec()[inode];
    dist = pow(node.Coord(0)-fCenterPoint[0],2.)+pow(node.Coord(1)-fCenterPoint[1],2.);
    dist += pow(node.Coord(2)-fCenterPoint[2],2.);
    dist = sqrt(dist);
    if(maxdist < dist) maxdist = dist;
  }
  return maxdist;
}

void TPZCompElDisc::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv,
                                 TPZFMatrix &phi, TPZFMatrix &dphix){
  TPZGeoEl * ref = this->Reference();
  if (!ref){
    PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
    return;
  }//if
  ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
  axes.Identity();//discontinuous shape does not use axes
  TPZManVector<REAL,3> XYZ(3);
  ref->X(intpoint, XYZ);
  this->ShapeX(XYZ,phi,dphix);
}

void TPZCompElDisc::ShapeX(TPZVec<REAL> &X, TPZFMatrix &phi, TPZFMatrix &dphi) {

  const int Degree = this->Degree();
  if(Degree < 0) return;

  if(Dimension()==0){
    TPZShapeDisc::Shape0D(fConstC,fCenterPoint,X,Degree,phi,dphi);
    return;
  }

  if(Dimension()==1){
     TPZShapeDisc::Shape1D(fConstC,fCenterPoint,X,Degree,phi,dphi);
     return;
  }

  if(Dimension()==2){
     TPZShapeDisc::Shape2D/*Full*/(fConstC,fCenterPoint,X,Degree,phi,dphi,fShapefunctionType);
  }

  if(Dimension()==3){
     TPZShapeDisc::Shape3D(fConstC,fCenterPoint,X,Degree,phi,dphi,fShapefunctionType);
     return;
  }

}//method

void TPZCompElDisc::Print(std::ostream &out) {

  out << "\nDiscontinous element : \n";
  if(Reference()) out << "\tGeometric reference index : " << Reference()->Index() << endl;
  out << "\tMaterial id : " << Reference()->MaterialId() << endl
      << "\tDegree of interpolation : " <<  this->Degree() << endl
      << "\tConnect index : " << fConnectIndex << endl
      << "\tNormalizing constant : " << fConstC << endl
      << "\tCenter point of the element : ";
  int size = fCenterPoint.NElements(),i;
  for(i=0;i<size-1;i++) out << fCenterPoint[i] << " , ";
  out << fCenterPoint[i] << endl;
  out << "\tDimension : " << this->Dimension() << endl;
}

int TPZCompElDisc::ConnectIndex(int side) {

  return fConnectIndex;
}

int TPZCompElDisc::NConnects(){


  return (fConnectIndex !=-1);

}

int TPZCompElDisc::CreateMidSideConnect(){
  // primeiro s� criados os elementos de volume depois os elementos BC associados aos seus lados
  // num est�io inicial o elemento BC �acoplado ao elemento ELV de volume de tal forma
  // que ambos s� vizinhos
  // o elemento BC n� pode ser dividido se o elemento ELV associado n� for dividido primeiro
  // caso o elemento ELV �dividido, ent� o elemento BC associado deveria ser dividido
  // tambem para manter a CC consistente com a malha
  // caso ELV �dividido e BC n� �ent� ELV �LowerLevelElement do elemento BC
  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    PZError << "\nTPZCompElDisc::CreateMidSideConnect Material nulo\n";
    return -1;
  }

  TPZGeoEl *ref = Reference();
  TPZStack<TPZCompElSide> list;
  int nsides = ref->NSides();
  int dimgrid = material->Dimension();
  int dim = Dimension();
  int existsconnect = 0;

  if(dimgrid == dim){
    //este �um elemento de volume
    //procura-se elemento superposto
    TPZCompElSide(this,nsides-1).EqualLevelElementList(list,0,0);
    int size = list.NElements(),i;
    for(i=0;i<size;i++){
      int dimel = list[i].Element()->Reference()->Dimension();
      if(dimel == dimgrid){
	int connectindex = list[i].Element()->ConnectIndex(0);
	list[i].Element()->SetConnectIndex(0,connectindex);
	existsconnect = 1;
	break;
      }
    }
  }

  if(dim != dimgrid/* - 1*/){ //dimgrid - 1 = interface dimension
    // o atual �um elemento BC
    fConnectIndex = -1;//=> return NshapeF() = 0
    return fConnectIndex;
  }

  if(!existsconnect){
    //o atual �um elemento de volume e
    //n� achou-se um elemento superposto
    int nvar = Material()->NStateVariables();
    int newnodeindex = Mesh()->AllocateNewConnect();
    SetConnectIndex(0,newnodeindex);
    TPZConnect &newnod = Mesh()->ConnectVec()[newnodeindex];
    int seqnum = newnod.SequenceNumber();
    const int nshape = this->NShapeF();
    Mesh()->Block().Set(seqnum,nvar*nshape);
    Mesh()->ConnectVec()[fConnectIndex].IncrementElConnected();
  }

  return fConnectIndex;
}

int TPZCompElDisc::NShapeF(){
  if(fConnectIndex == -1) return 0;
  //deve ter pelo menos um connect

  int dim = Dimension();
  return TPZShapeDisc::NShapeF(this->Degree(),dim,fShapefunctionType);
}

int TPZCompElDisc::NConnectShapeF(int inod){
#ifdef DEBUG
  if (inod != 0){
    PZError << "\nFATAL ERROR AT " << __PRETTY_FUNCTION__ 
            << " - TPZCompElDisc has only one connect and inod = " << inod << "\n";
  }
#endif
  return this->NShapeF();
}

void TPZCompElDisc::InternalPoint(TPZVec<REAL> &point){
  //ponto deformado
  point.Resize(3,0.);
  point[0] = fCenterPoint[0];
  point[1] = fCenterPoint[1];
  point[2] = fCenterPoint[2];
}


void TPZCompElDisc::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){

  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    cout << "TPZCompElDisc::CalcStiff : no material for this element\n";
    ek.Reset();
    ef.Reset();
    return;
  }
  TPZMaterialData data;
//  this->ComputeShap
  TPZGeoEl *ref = Reference();
  int ncon = NConnects();
  int dim = Dimension();
  int nstate = material->NStateVariables();
  int nshape = NShapeF();
  int numeq = nshape * nstate;

  ek.fMat.Redim(numeq,numeq);
  ef.fMat.Redim(numeq,1);
  if(ncon){
    ek.fBlock.SetNBlocks(ncon);
    ef.fBlock.SetNBlocks(ncon);
    ek.fBlock.Set(0,NShapeF()*nstate);
    ef.fBlock.Set(0,NShapeF()*nstate);
  }
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);
  for(int i=0;i<ncon;i++){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }
  if(ncon==0) return;//elemento CC no passa
  TPZFMatrix phix(nshape,1),dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL detjac,weight;
  int integ = max( 2 * Degree(), 0);
  TPZIntPoints *intrule = 0;
  intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,integ);
  if(material->HasForcingFunction())
  {
     int maxint = intrule->GetMaxOrder();
     TPZManVector<int> order(Reference()->Dimension());
     intrule->GetOrder(order);
     order.Fill(maxint);
     intrule->SetOrder(order);
  }

  int npoints = intrule->NPoints(), ip;
  TPZManVector<REAL,220> sol(nstate,0.);
  TPZFNMatrix<660> dsol(dim,nstate,0.);

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    ref->X(intpoint, x);
    weight *= fabs(detjac);
    ShapeX(x,phix,dphix);
    axes.Identity();

    if (material->NeedsSolutionToContribute()){
      this->ComputeSolution(intpoint, phix, dphix, axes, sol, dsol);
    }

    material->Contribute(x,jacinv,sol,dsol,weight,axes,phix,dphix,ek.fMat,ef.fMat);
  }
}
void TPZCompElDisc::CalcResidual(TPZElementMatrix &ef){
  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    cout << "TPZCompElDisc::CalcStiff : no material for this element\n";
    ef.Reset();
    return;
  }
  TPZGeoEl *ref = Reference();
  int ncon = NConnects();
  int dim = Dimension();
  int nstate = material->NStateVariables();
  int nshape = NShapeF();
  int numeq = nshape * nstate;

  // clean ef

  ef.fMat.Redim(numeq,1);
  if(ncon){//pode serr no m�imo ncon = 1
    ef.fBlock.SetNBlocks(ncon);
    ef.fBlock.Set(0,NShapeF()*nstate);
  }
  ef.fConnect.Resize(ncon);
  for(int i=0;i<ncon;i++){
    (ef.fConnect)[i] = ConnectIndex(i);
  }
  if(ncon==0) return;//elemento CC no passa
  TPZFMatrix phix(nshape,1),dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL detjac,weight;
  int integ = 2*Degree();
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,integ);
  int npoints = intrule->NPoints(),ip;                                              //integra fi*fj
  TPZVec<REAL> sol(nstate,0.);
  TPZFMatrix dsol(dim,nstate,0.);

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    ref->X(intpoint, x);
    weight *= fabs(detjac);
    ShapeX(x,phix,dphix);
    axes.Identity();

    this->ComputeSolution(intpoint, phix, dphix, axes, sol, dsol);

    material->Contribute(x,jacinv,sol,dsol,weight,axes,phix,dphix,ef.fMat);
  }
}

REAL TPZCompElDisc::SizeOfElement()
{
  TPZGeoEl *ref = Reference();

  int dim = ref->Dimension();
  int side = ref->NSides()-1;
  if(dim == 2) ref->SideArea(side);
  if(!dim || dim > 2){
    PZError << "TPZCompElDisc::SizeOfElement case not permited\n";
    return 0.;
  }
  if(dim == 1){
    TPZGeoNode node0 = Mesh()->Reference()->NodeVec()[ref->NodeIndex(0)];
    TPZGeoNode node1 = Mesh()->Reference()->NodeVec()[ref->NodeIndex(1)];
    TPZVec<REAL> no0(3),no1(3);
    for(int i=0;i<3;i++){
      no0[i] = node0.Coord(i);
      no1[i] = node1.Coord(i);
    }
    return ref->Distance(no0,no1);
  }
  PZError << "TPZCompElDisc::SizeOfElement this in case that it is not contemplated\n";
  return 0.;
}

void TPZCompElDisc::Divide(int index,TPZVec<int> &subindex,int interpolatesolution){

  if (Mesh()->ElementVec()[index] != this) {
    PZError << "TPZInterpolatedElement::Divide index error";
    subindex.Resize(0);
    return;
  }
  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material)
  {
    PZError << __PRETTY_FUNCTION__ << " no material\n";
    return;
  }

  TPZGeoEl *ref = Reference();
  RemoveInterfaces();

  if(0){//TESTE
    ofstream mesh("MALHADIV.out");//TESTE
    Mesh()->Reference()->Print(mesh);//TESTE
    Mesh()->Print(mesh);//TESTE
    mesh.flush();  //TESTE
    mesh.close();//TESTE
  }//TESTE

  //divide o elemento geom�rico
  int nsubs = ref->NSubElements();
  subindex.Resize(nsubs);
  TPZManVector<TPZGeoEl *> geosubs(nsubs);
  ref->Divide(geosubs);
  if(!geosubs.NElements()) {
    subindex.Resize(0);
    return;
  }

  this->Mesh()->ElementVec()[index] = NULL;
  ref->ResetReference();
  TPZCompElDisc *discel;
  int i,deg;
  deg = this->Degree();

  for (i=0;i<nsubs;i++){
    geosubs[i]->CreateCompEl(*Mesh(),subindex[i]);//aqui
    //new TPZCompElDisc(*Mesh(),geosubs[i],subindex[i]);
    discel = dynamic_cast<TPZCompElDisc *> (Mesh()->ElementVec()[subindex[i]]);
    if (!discel){
      std::stringstream mess;
      mess << __PRETTY_FUNCTION__ << " - discel is NULL ";
      LOGPZ_ERROR(logger, mess.str() );
      continue;
    }
    discel->SetDegree(deg);
  }

  if (interpolatesolution){
    Mesh()->ExpandSolution();
    for(i=0; i<nsubs; i++) {
      discel = dynamic_cast<TPZCompElDisc *> ( Mesh()->ElementVec()[subindex[i]] );
      if (!discel){
        std::stringstream mess;
        mess << __PRETTY_FUNCTION__ << " - discel is NULL ";
        LOGPZ_ERROR(logger, mess.str() );
        continue;
      }
      if(discel->Dimension() < material->Dimension()) continue;//elemento BC
      discel->InterpolateSolution(*this);
    }
  }//if interpolate

  delete this;
}

void TPZCompElDisc::InterpolateSolution(TPZCompElDisc &coarsel){
  // accumulates the transfer coefficients between the current element and the
  // coarse element into the transfer matrix, using the transformation t
  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material)
  {
    std::cout << __PRETTY_FUNCTION__ << " No material " << std::endl;
    return;
  }
  TPZGeoEl *ref = Reference();
  TPZTransform t(Dimension());
  TPZGeoEl * coarseref = coarsel.Reference();
  t = Reference()->BuildTransform2(ref->NSides()-1,coarseref,t);

  int locmatsize = NShapeF();
  int cormatsize = coarsel.NShapeF();
  int nvar = material->NStateVariables();
  int dimension = Dimension();

  TPZFMatrix loclocmat(locmatsize,locmatsize,0.);
  TPZFMatrix projectmat(locmatsize,nvar,0.);

  TPZVec<int> prevorder(dimension),order(dimension);
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,2*Degree());
//  int i;

  TPZFMatrix locphi(locmatsize,1);
  TPZFMatrix locdphi(dimension,locmatsize);	// derivative of the shape function
  // in the master domain

  TPZFMatrix corphi(cormatsize,dimension);
  TPZFMatrix cordphi(dimension,cormatsize);	// derivative of the shape function
  // in the master domain

  TPZVec<REAL> int_point(dimension),coarse_int_point(dimension);
  TPZFMatrix jacobian(dimension,dimension),jacinv(dimension,dimension);
  TPZFMatrix axes(3,3,0.);
  REAL zero = 0.;
  TPZVec<REAL> x(3,zero);
  TPZVec<REAL> u(nvar);

  int numintpoints = intrule->NPoints();
  REAL weight;
  int lin,ljn,cjn;
  TPZConnect *df;
  TPZBlock &block = Mesh()->Block();

  for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {

    intrule->Point(int_ind,int_point,weight);
    REAL jac_det = 1.;
    Reference()->Jacobian( int_point, jacobian , axes,jac_det,jacinv);
    Reference()->X(int_point, x);
    ShapeX(x,locphi,locdphi);
    axes.Identity();

    weight *= jac_det;
    t.Apply(int_point,coarse_int_point);
    TPZVec<REAL> x2(3,zero);
    coarsel.Reference()->X(coarse_int_point, x2);

#ifdef DEBUG
    if (x.NElements() == x2.NElements()){
      for(int itest = 0; itest < x.NElements(); itest++){
        if ( fabs(x[itest] - x2[itest]) > 1.e-12 ){
          std::stringstream mess;
          mess << __PRETTY_FUNCTION__ << " - x[" << itest << "] != x2[" << itest << "]";
          LOGPZ_ERROR(logger, mess.str() );
        }
      }
    }
    else{
      std::stringstream mess;
      mess << __PRETTY_FUNCTION__ << " - x.NElements() != x2.NElements()";
      LOGPZ_ERROR(logger, mess.str() );
    }
#endif

    coarsel.ShapeX(x,corphi,cordphi);
    //u eh a solucao no elemento maior
    u.Fill(0.);
    int iv = 0;

    df = &coarsel.Connect(0);
    int dfseq = df->SequenceNumber();

    int dfvar = block.Size(dfseq);
    for(ljn=0; ljn<dfvar; ljn++) {
      u[iv%nvar] += corphi(iv/nvar,0)*block(dfseq,0,ljn,0);
      iv++;
    }

    for(lin=0; lin<locmatsize; lin++) {
      for(ljn=0; ljn<locmatsize; ljn++) {
	loclocmat(lin,ljn) += weight*locphi(lin,0)*locphi(ljn,0);
      }
      for(cjn=0; cjn<nvar; cjn++) {
	projectmat(lin,cjn) += weight*locphi(lin,0)*u[cjn];
      }
    }
    jacobian.Zero();
  }
  loclocmat.SolveDirect(projectmat,ELU);
  // identify the non-zero blocks for each row
  int iv=0;

  df = &Connect(0);
  int dfseq = df->SequenceNumber();
  int dfvar = block.Size(dfseq);
#ifdef DEBUG
{
  int aux1 = NShapeF() * material->NStateVariables();
  int aux2 = material->NStateVariables();
  if ((dfvar != aux1) || (nvar != aux2) ){
  LOGPZ_ERROR(logger, __PRETTY_FUNCTION__ );
  }
}
#endif
  for(ljn=0; ljn<dfvar; ljn++) {
    block(dfseq,0,ljn,0) = projectmat(iv/nvar,iv%nvar);
    iv++;
  }
  delete intrule;
}

void TPZCompElDisc::Solution(TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol) {
  //#ifdef _AUTODIFF
  //  TPZConservationLaw2 *mat = dynamic_cast<TPZConservationLaw2 *>(fMaterial);
  //#else
  //  TPZConservationLaw *mat = dynamic_cast<TPZConservationLaw *>(fMaterial);
  //#endif

  if(var >= 100) {
    TPZCompEl::Solution(qsi,var,sol);
    return;
  }
  int nshape = NShapeF();
  int dim = Dimension();
  if(var == 99) {
    sol[0] = Degree();
    return;
  }

  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    PZError << "TPZIntEl::Solution : no Material for this element\n";
    Print(PZError);
    return;
  }

  TPZGeoEl *ref = Reference();
  int numdof = material->NStateVariables();
  TPZFNMatrix<220> phi(nshape,1);
  TPZFNMatrix<660> dphi(dim,nshape);
  TPZManVector<REAL> u(numdof);
  TPZFMatrix du(dim,numdof,0.);
  TPZFMatrix axes(3,3,0.);
  axes.Identity();
  REAL jacstore[10],jacinvstore[10];
  TPZFMatrix jacobian(dim,dim,jacstore,10);
  TPZFMatrix jacinv(dim,dim,jacinvstore,10);
  TPZManVector<REAL> x(3);

  if(var >= 0){
    ref->X(qsi,x);
  } else if(var < 0){
    //neste caso 0 ponto qsi est�no elemento deformado
    var *= -1;//recuperando var
    for(int i=0;i<3;i++) x[i] = qsi[i];
  }
  ShapeX(x,phi,dphi);
  this->ComputeSolution(qsi, phi, dphi, axes, u, du);
  material->Solution(u,du,axes,var,sol);
}

void TPZCompElDisc::Solution(TPZVec<REAL> &x, TPZVec<REAL> &uh){

  TPZCompMesh *finemesh = Mesh();
  TPZBlock &fineblock = finemesh->Block();
  int nstate = Material()->NStateVariables();
  TPZFMatrix &FineMeshSol = finemesh->Solution();
  int matsize = NShapeF(),dim = Dimension();
  TPZFMatrix phix(matsize,1,0.);
  TPZFMatrix dphix(dim,matsize,0.);
  ShapeX(x,phix,dphix);
  TPZConnect *df = &Connect(0);
  int dfseq = df->SequenceNumber();
  int dfvar = fineblock.Size(dfseq);
  int pos   = fineblock.Position(dfseq);
  int iv = 0,d;
  uh.Fill(0.);
  for(d=0; d<dfvar; d++) {
    uh[iv%nstate] += phix(iv/nstate,0)*FineMeshSol(pos+d,0);
    iv++;
  }
}



void TPZCompElDisc::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension)
{
  TPZGeoEl *ref = Reference();
  int mat = Material()->Id();
  int nsides = ref->NSides();

  if(dimension == 2 && mat > 0){
    if(nsides == 9){
      new TPZGraphElQ2dd(this,&grmesh);
      return;
    }
    if(nsides == 7){
      new TPZGraphElTd(this,&grmesh);
      return;
    }
  }
  if(dimension == 3 && mat > 0){
    new TPZGraphElQ3dd(this,&grmesh);
  }
  if(dimension == 1 && mat > 0){
    new TPZGraphEl1dd(this,&grmesh);
  }
}

int TPZCompElDisc::NSides(){

  return Reference()->NSides();
}

int TPZCompElDisc::NInterfaces(){

  int nsides = this->NSides();

  switch( nsides )
    {
    case 3: //line
      return 2;
      break;

    case 7: //triangle
      return 3;
      break;

    case 9: //square
      return 4;

    case 15: // Tetrahedra.
      return 4;
      break;

    case 19: // Prism.
      return 5;
      break;

    case 21: // Pyramid.
      return 6;
      break;

    case 27: // Hexaedra.
      return 8;
      break;

    default:
      PZError << "TPZCompElDisc::NFaces() - Unknown element shape!" << endl;
      exit (-1);
    }
}

//#include "TPZAgglomerateEl.h"
void TPZCompElDisc::AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight){

  int i,npoints;
  TPZVec<REAL> pt(3),x(3,0.0);
  TPZFMatrix jacobian(3,3),jacinv(3,3),axes(3,3);
  REAL detjac,wt;

  TPZGeoEl *subgel = Reference();
  if(!subgel) PZError << "TPZCompElDisc::AccumulateIntegrationRule data error, null geometric reference\n";
  TPZIntPoints *rule = subgel->CreateSideIntegrationRule(subgel->NSides()-1,degree);
  npoints = rule->NPoints();

  for(i=0;i<npoints;i++){

    rule->Point(i,pt,wt);
    subgel->Jacobian(pt,jacobian,axes,detjac,jacinv);
    subgel->X(pt, x);

    point.Push(x[0]);
    point.Push(x[1]);
    point.Push(x[2]);

    weight.Push(wt * fabs(detjac));
  }
  delete rule;
}


void TPZCompElDisc::CenterPoint(TPZVec<REAL> &center){

  TPZGeoEl *ref = Reference();
  if(ref || Type() == EDiscontinuous){
    ref->CenterPoint(ref->NSides()-1,center);
    return;
  } else {//aglomerado
//     TPZStack<TPZCompEl *> elvec;
//     dynamic_cast<TPZAgglomerateElement *>(this)->ListOfDiscEl(elvec);
//     TPZGeoEl *ref = elvec[0]->Reference();
//     ref->CenterPoint(ref->NSides()-1,center);
    PZError << "TPZCompElDisc::CenterPoint center points not exists!\n";
  }
}



void TPZCompElDisc::EvaluateError(  void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
				    TPZVec<REAL> &errors,TPZBlock * /*flux */) {

  int NErrors = this->Material()->NEvalErrors();
  errors.Resize(NErrors);
  errors.Fill(0.);
  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    PZError << "TPZInterpolatedElement::EvaluateError : no material for this element\n";
    Print(PZError);
    return;
  }
  int problemdimension = Mesh()->Dimension();
  if(Reference()->Dimension() < problemdimension) return;
  // Adjust the order of the integration rule
  int nsides = Reference()->NSides();
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(nsides-1,20);
  TPZGeoEl *ref = Reference();
  int dim = Dimension();

  int ndof = material->NStateVariables();
  int nflux = material->NFluxes();
  int nshape = NShapeF();
  //suficiente para ordem 5 do cubo
  TPZFNMatrix<220> phi(nshape,1);
  TPZFNMatrix<660> dphi(dim,nshape),dphix(dim,nshape);
  TPZFNMatrix<9> jacobian(dim,dim),jacinv(dim,dim);
  TPZFNMatrix<9> axes(3,3);
  TPZManVector<REAL,3> x(3);//TPZVec<REAL> x(3,0.);
  TPZManVector<REAL,6> u_exact(ndof);
  TPZFNMatrix<90> du_exact(dim,ndof);
  TPZManVector<REAL,3> intpoint(3), values(NErrors);
  values.Fill(0.0);
  REAL detjac,weight;
  TPZManVector<REAL,6> u(ndof);
  TPZFNMatrix<90> dudx(dim,ndof);
  TPZManVector<REAL,9> flux_el(nflux,0.);
  TPZAutoPointer<TPZMaterial> matp = material;

  for(int nint=0; nint<intrule->NPoints(); nint++) {

    intrule->Point(nint,intpoint,weight);
    ref->Jacobian( intpoint , jacobian, axes, detjac , jacinv);
    ref->X( intpoint , x);
    ShapeX(x,phi,dphix);

    //it is possible to have more than first derivative, as Laplacian and its derivatives.
    int nderiv = dphix.Rows();
    du_exact.Redim(nderiv, ndof);
    dudx.Redim(nderiv, ndof);

    axes.Identity();
    weight *= fabs(detjac);

    this->ComputeSolution(intpoint, phi, dphix, axes, u, dudx);

    //contribu�es dos erros
    if(fp) {
      fp(x,u_exact,du_exact);
      matp->Errors(x,u,dudx,axes,flux_el,u_exact,du_exact,values);
      for(int ier = 0; ier < NErrors; ier++)
	errors[ier] += values[ier]*weight;
    }
  }//fim for : integration rule
   //Norma sobre o elemento
  for(int ier = 0; ier < NErrors; ier++)
    errors[ier] = sqrt(errors[ier]);
  delete intrule;
}


void TPZCompElDisc::BuildTransferMatrix(TPZCompElDisc &coarsel, TPZTransfer &transfer){
  // accumulates the transfer coefficients between the current element and the
  // coarse element into the transfer matrix, using the transformation t

  int locnshape = NShapeF();
  int cornshape = coarsel.NShapeF();

  // compare interpolation orders
  // the interpolation order of this >= that interpolation order of coarse
  int locdeg = Degree(), coarsedeg = coarsel.Degree();
  if(coarsedeg > locdeg) {
    SetDegree(coarsedeg);
  }

  TPZFNMatrix<500> loclocmat(locnshape,locnshape);
  TPZFMatrix loccormat(locnshape,cornshape);
  loclocmat.Zero();
  loccormat.Zero();

  TPZGeoEl *ref = Reference();
  int integdeg = locdeg >= coarsedeg ? locdeg : coarsedeg;
  TPZIntPoints *intrule = ref->CreateSideIntegrationRule(ref->NSides()-1,2*integdeg);
  int dimension = Dimension();

  TPZFNMatrix<50> locphi(locnshape,1);
  TPZFNMatrix<150> locdphi(dimension,locnshape);
  locphi.Zero();
  locdphi.Zero();
  // derivative of the shape function
  // in the master domain

  TPZFMatrix corphi(cornshape,1);
  TPZFMatrix cordphi(dimension,cornshape);
  // derivative of the shape function
  // in the master domain

  TPZManVector<REAL> int_point(dimension);
  TPZFNMatrix<9> jacobian(dimension,dimension);
  TPZFMatrix jacinv(dimension,dimension);
  TPZFNMatrix<9> axes(3,3);
  TPZManVector<REAL> x(3);

  int_point.Fill(0.,0);
  REAL jac_det = 1.;
  ref->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
  REAL multiplier = 1./jac_det;

  int numintpoints = intrule->NPoints();
  REAL weight;
  int lin,ljn,cjn;

  for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {

    intrule->Point(int_ind,int_point,weight);
    ref->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
    ref->X(int_point, x);
    ShapeX(int_point,locphi,locdphi);
    weight *= jac_det;
    corphi.Zero();
    cordphi.Zero();
    coarsel.ShapeX(int_point,corphi,cordphi);

    for(lin=0; lin<locnshape; lin++) {
      for(ljn=0; ljn<locnshape; ljn++) {
	loclocmat(lin,ljn) += weight*locphi(lin,0)*locphi(ljn,0)*multiplier;
      }
      for(cjn=0; cjn<cornshape; cjn++) {
	loccormat(lin,cjn) += weight*locphi(lin,0)*corphi(cjn,0)*multiplier;
      }
    }
    jacobian.Zero();
  }
  loclocmat.SolveDirect(loccormat,ELDLt);


  int locblockseq = Connect(0).SequenceNumber();
  TPZStack<int> globblockvec;
  int numnonzero = 0;
  int cind = coarsel.ConnectIndex(0);
  TPZConnect &con = coarsel.Mesh()->ConnectVec()[cind];
  int corblockseq = con.SequenceNumber();
  if(locnshape == 0 || cornshape == 0)
    PZError << "TPZCompElDisc::BuilTransferMatrix error I\n";
  TPZFMatrix small(locnshape,cornshape,0.);
  loccormat.GetSub(0,0,locnshape,cornshape,small);
  REAL tol = Norm(small);
  if(tol >= 1.e-10) {
    globblockvec.Push(corblockseq);
    numnonzero++;
  }
  if(!numnonzero)
    PZError << "TPZCompElDisc::BuilTransferMatrix error II\n";
  if(transfer.HasRowDefinition(locblockseq))
    PZError << "TPZCompElDisc::BuilTransferMatrix error III\n";
  transfer.AddBlockNumbers(locblockseq,globblockvec);
  if(cornshape == 0 || locnshape == 0)
    PZError << "TPZCompElDisc::BuilTransferMatrix error IV\n";
  loccormat.GetSub(0,0,locnshape,cornshape,small);
  transfer.SetBlockMatrix(locblockseq,globblockvec[0],small);

  SetDegree(locdeg);
  delete intrule;
}

void TPZCompElDisc::AccumulateVertices(TPZStack<TPZGeoNode *> &nodes) {
  TPZGeoEl *geo = Reference();

//Code isnt place to chat
//#warning "Este metodo nao funciona para aglomerados contendo aglomerados"
  if(!geo) {
    PZError <<  "TPZCompElDisc::AccumulateVertices null reference\n";
    return;
  }
  int nvertices = geo->NNodes();
  int l;
  for(l=0;l<nvertices;l++) nodes.Push( geo->NodePtr(l) );
}

void TPZCompElDisc::SetDegree(int degree) {
  if (fConnectIndex < 0) return;
  TPZConnect &c = Mesh()->ConnectVec()[fConnectIndex];
  c.SetOrder(degree);
  int seqnum = c.SequenceNumber();
  int nvar = 1;
  TPZAutoPointer<TPZMaterial> mat = Material();
  if(mat) nvar = mat->NStateVariables();
  Mesh()->Block().Set(seqnum,NShapeF()*nvar);
}

  /**
  * returns the unique identifier for reading/writing objects to streams
  */
int TPZCompElDisc::ClassId() const
{
  return TPZCOMPELDISCID;
}

template class
    TPZRestoreClass< TPZCompElDisc, TPZCOMPELDISCID>;

  /**
  Save the element data to a stream
  */
void TPZCompElDisc::Write(TPZStream &buf, int withclassid)
{
  TPZCompEl::Write(buf,withclassid);
  WriteObjects(buf,fCenterPoint);
  buf.Write(&fConnectIndex,1);
  buf.Write(&fConstC,1);
  int matid = Material()->Id();
  buf.Write(&matid,1);
  int shapetype = fShapefunctionType;
  buf.Write(&shapetype,1);

}

  /**
  Read the element data from a stream
  */
 void TPZCompElDisc::Read(TPZStream &buf, void *context)
 {
  TPZCompEl::Read(buf,context);
  ReadObjects<3>(buf,fCenterPoint);
  buf.Read(&fConnectIndex,1);
  buf.Read(&fConstC,1);
  int matid;
  buf.Read(&matid,1);
//  fMaterial = Mesh()->FindMaterial(matid);
  int shapetype;
  buf.Read(&shapetype,1);
  fShapefunctionType = (TPZShapeDisc::MShapeType) shapetype;
 }

void TPZCompElDisc::ComputeError(int errorid, TPZVec<REAL> &error){

  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    cout << "TPZCompElDisc::ComputeError : no material for this element\n";
    return;
  }

  TPZVec<REAL> x(3,0.);
  REAL weight;
  int dim = Dimension();
  TPZVec<REAL> intpoint(dim,0.);
  int integ = max( 2 * Degree(), 0);
  TPZIntPoints *intrule = 0;
  intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,integ);
  int nstate = material->NStateVariables();
  int npoints = intrule->NPoints(), ip;

  TPZManVector<REAL,220> sol(nstate,0.);
  TPZFNMatrix<660> dsol(dim,nstate,0.);
  TPZGeoEl * ref = this->Reference();

  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;

  const int POrder = this->Degree();
  const REAL faceSize = 2. * ref->ElementRadius();
  error.Fill(0.);
  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    ref->X(intpoint, x);
    weight *= fabs(detjac);
    this->ComputeSolution(intpoint, sol, dsol, axes);
    material->ContributeErrors(x,sol, dsol,weight,error,POrder,faceSize, errorid);
  }
  delete intrule;
}

void TPZCompElDisc::Integrate(int variable, TPZVec<REAL> & value){
  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    cout << "TPZCompElDisc::ComputeError : no material for this element\n";
    return;
  }
  const int dim = this->Dimension();
  REAL weight;
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZManVector<REAL, 3> intpoint(dim,0.);
  const int varsize = material->NSolutionVariables(variable);
  TPZManVector<REAL> sol(varsize);

  value.Resize(varsize);
  value.Fill(0.);
  int integ = max( 2 * this->Degree(), 0);
  TPZIntPoints *intrule = this->Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,integ);
  int npoints = intrule->NPoints(), ip, iv;
  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    this->Reference()->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    sol.Fill(0.);
    this->Solution(intpoint, variable, sol);
    for(iv = 0; iv < varsize; iv++) value[iv] += sol[iv]*weight*fabs(detjac);
  }
  delete intrule;
}

void TPZCompElDisc::ComputeSolution(TPZVec<REAL> &qsi, TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix & axes){
  TPZGeoEl * ref = this->Reference();
  const int nshape = this->NShapeF();
  const int dim = ref->Dimension();
  TPZFMatrix phix(nshape,1),dphix(dim,nshape);

  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZManVector<REAL,3> x(3,0.);
  ref->Jacobian( qsi, jacobian, axes, detjac , jacinv);
  ref->X(qsi, x);
  this->ShapeX(x,phix,dphix);
  this->ComputeSolution(qsi, phix, dphix, axes, sol, dsol);
}//method

void TPZCompElDisc::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
                                    const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol){

  const int nstate = this->Material()->NStateVariables();
  const int ncon = this->NConnects();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();

  sol.Resize(nstate);
  sol.Fill(0.);
  dsol.Redim(dphix.Rows(), nstate);
  dsol.Zero();

  int iv = 0, d;
  for(int in=0; in<ncon; in++) {
    TPZConnect *df = &Connect(in);
    int dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    int pos = block.Position(dfseq);
    for(int jn=0; jn<dfvar; jn++) {
      sol[iv%nstate] += phi(iv/nstate,0)*MeshSol(pos+jn,0);
      for(d=0; d<dphix.Rows(); d++){
        dsol(d,iv%nstate) += dphix(d,iv/nstate)*MeshSol(pos+jn,0);
      }
      iv++;
    }
  }

}//method

void TPZCompElDisc::ComputeSolution(TPZVec<REAL> &qsi,
                                    TPZVec<REAL> &normal,
                                    TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                                    TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes){
  //TPZCompElDisc has no left/right elements. Only interface elements have it.
  leftsol.Resize(0);
  dleftsol.Resize(0,0);
  leftaxes.Zero();
  rightsol.Resize(0);
  drightsol.Resize(0,0);
  rightaxes.Zero();
  normal.Resize(0);
}//method

