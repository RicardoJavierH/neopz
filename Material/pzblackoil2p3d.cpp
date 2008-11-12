//$Id: pzblackoil2p3d.cpp,v 1.1 2008-11-12 12:45:51 fortiago Exp $ 

#include "pzblackoil2p3d.h"
#include "pzbndcond.h"

using namespace std;

#ifdef _AUTODIFF

#include "fad.h"


TPZBlackOil2P3D::EState TPZBlackOil2P3D::gState = ELastState;

/** 
 * Interpolacao linear
 */
void TPZBlackOil2P3D::Interpolate(std::map<REAL,REAL> &dados, double x, double &y, double &dy){
  double x0, xL, y0, yL;
  std::map< REAL, REAL >::iterator w;
  w = dados.lower_bound(x);

  if(dados.size() < 2){
    std::cout << "Error at " << __PRETTY_FUNCTION__ << "\n";
    std::cout.flush();
  }

  if(w == dados.end()){
    w--;
    y = w->second;

    ///derivada
    xL = w->first;
    yL = w->second;
    w--;
    x0 = w->first;
    y0 = w->second;
    dy = (yL-y0)/(xL-x0);

    return;
  }
  if (w == dados.begin()){
    y = w->second;

    ///derivada
    x0 = w->first;
    y0 = w->second;
    w++;
    xL = w->first;
    yL = w->second;
    dy = (yL-y0)/(xL-x0);
    return;
  }

  xL = w->first;
  yL = w->second;
  w--;

  x0 = w->first;
  y0 = w->second;

  y = (yL-y0)*(x-x0)/(xL-x0)+y0;
  dy = (yL-y0)/(xL-x0);

}

void TPZBlackOil2P3D::Interpolate(std::map<REAL,REAL> &dados, BFadREAL x, BFadREAL &y){
  double x0, xL, y0, yL;
  std::map< REAL, REAL >::iterator w;
  w = dados.lower_bound(x.val());

  if(dados.size() < 2){
    std::cout << "Error at " << __PRETTY_FUNCTION__ << "\n";
    std::cout.flush();
  }

  if(w == dados.end()){
    w--;
    y = w->second;

    return;
  }
  if (w == dados.begin()){
    y = w->second;

    return;
  }

  xL = w->first;
  yL = w->second;
  w--;

  x0 = w->first;
  y0 = w->second;

  y = (yL-y0)*(x-x0)/(xL-x0)+y0;

}

///Dados
/**
  * Permeabilidade relativa do oleo
  * Kro = Kro( Sw )
  */
void TPZBlackOil2P3D::Kro(double So, double &Kro, double &dKroSo){
  const int n = 8;
  double tabela[n][2] = {{0.12,1.},{0.2,0.8},{0.3,0.6},{0.4,0.45},{0.5,0.3},{0.6,0.2},{0.7,0.1},{0.82,0.}};
  std::map<REAL,REAL> dados;
  for(int i = 0; i < n; i++){
    dados[ 1.-tabela[i][0] ] = tabela[i][1];
  }
  this->Interpolate(dados,So,Kro,dKroSo);
}
void TPZBlackOil2P3D::Kro(BFadREAL So, BFadREAL &Kro){
  const int n = 8;
  double tabela[n][2] = {{0.12,1.},{0.2,0.8},{0.3,0.6},{0.4,0.45},{0.5,0.3},{0.6,0.2},{0.7,0.1},{0.82,0.}};
  std::map<REAL,REAL> dados;
  for(int i = 0; i < n; i++){
    dados[ 1.-tabela[i][0] ] = tabela[i][1];
  }
  this->Interpolate(dados,So,Kro);
}

/**
  * Permeabilidade relativa da agua
  * Krw = Krw( Sw )
  */
void TPZBlackOil2P3D::Krw(double So, double &Krw, double &dKrwSo){
  const int n = 8;
  double tabela[n][2] = {{0.12,0.},{0.2,0.1},{0.3,0.2},{0.4,0.3},{0.5,0.4},{0.6,0.55},{0.7,0.7},{0.82,1.}};
  std::map<REAL,REAL> dados;
  for(int i = 0; i < n; i++){
    dados[ 1.-tabela[i][0] ] = tabela[i][1];
  }
  this->Interpolate(dados,So,Krw,dKrwSo);
}

void TPZBlackOil2P3D::Krw(BFadREAL So, BFadREAL &Krw){
  const int n = 8;
  double tabela[n][2] = {{0.12,0.},{0.2,0.1},{0.3,0.2},{0.4,0.3},{0.5,0.4},{0.6,0.55},{0.7,0.7},{0.82,1.}};
  std::map<REAL,REAL> dados;
  for(int i = 0; i < n; i++){
    dados[ 1.-tabela[i][0] ] = tabela[i][1];
  }
  this->Interpolate(dados,So,Krw);
}

/** Bo = Bo( po )
  */
void TPZBlackOil2P3D::Bo(double po, double &Bo, double &dBoDpo){
  const int n = 10;
  double tabela[n][2] = {{14.7,1.062},{264.7,1.15},{514.7,1.207},{1014.7,1.295},{2014.7,1.435},
                         {2514.,1.5},{3014.7,1.565},{4014.7,1.695},{5014.7,1.827},{9014.7,2.352}};
  const double conv = 6894.757;
  std::map<REAL,REAL> dados;
  for(int i = 0; i < n; i++){
    dados[ tabela[i][0]*conv ] = tabela[i][1];
  }
  this->Interpolate(dados,po,Bo,dBoDpo);
}

void TPZBlackOil2P3D::Bo(BFadREAL po, BFadREAL &Bo){
  const int n = 10;
  double tabela[n][2] = {{14.7,1.062},{264.7,1.15},{514.7,1.207},{1014.7,1.295},{2014.7,1.435},
                         {2514.,1.5},{3014.7,1.565},{4014.7,1.695},{5014.7,1.827},{9014.7,2.352}};
  const double conv = 6894.757;
  std::map<REAL,REAL> dados;
  for(int i = 0; i < n; i++){
    dados[ tabela[i][0]*conv ] = tabela[i][1];
  }
  this->Interpolate(dados,po,Bo);
}

/** Viscosidade do oleo ViscOleo = ViscOleo( po )
  */
void TPZBlackOil2P3D::ViscOleo(double po, double &ViscOleo, double &dViscOleoDpo){
  const int n = 10;
  double tabela[n][2] = {{14.7,1.04},{264.7,0.975},{514.7,0.91},{1014.7,0.83},{2014.7,0.695},
                         {2514.,0.641},{3014.7,0.594},{4014.7,0.51},{5014.7,0.449},{9014.7,0.203}};
  const double conv = 6894.757;
  const double convVisc = 1e-3;
  std::map<REAL,REAL> dados;
  for(int i = 0; i < n; i++){
    dados[ tabela[i][0]*conv ] = tabela[i][1]*convVisc;
  }
  this->Interpolate(dados,po,ViscOleo,dViscOleoDpo);
}

void TPZBlackOil2P3D::ViscOleo(BFadREAL po, BFadREAL &ViscOleo){
  const int n = 10;
  double tabela[n][2] = {{14.7,1.04},{264.7,0.975},{514.7,0.91},{1014.7,0.83},{2014.7,0.695},
                         {2514.,0.641},{3014.7,0.594},{4014.7,0.51},{5014.7,0.449},{9014.7,0.203}};
  const double conv = 6894.757;
  const double convVisc = 1e-3;
  std::map<REAL,REAL> dados;
  for(int i = 0; i < n; i++){
    dados[ tabela[i][0]*conv ] = tabela[i][1]*convVisc;
  }
  this->Interpolate(dados,po,ViscOleo);
}

/** Pressao capilar
  * pc = pc( Sw )
  */
void TPZBlackOil2P3D::PressaoCapilar(double So, double &pc, double &DpcDSo){
  pc = 0.;
  DpcDSo = 0.;
}

void TPZBlackOil2P3D::PressaoCapilar(BFadREAL So, BFadREAL &pc){
  pc = 0.;
}

/** Porosidade
  * Phi = Phi( pw ) - fizemos como Phi ( po )
*/
void TPZBlackOil2P3D::Porosidade(double po, double &poros, double &dPorosDpo){
  const double comp = 3.625943e-10;
  const double pref = 101352.93;
  const double porosRef = 0.22;
  poros = porosRef*exp(comp*(po-pref));
  dPorosDpo = comp*porosRef*exp(comp*(po-pref));
}

void TPZBlackOil2P3D::Porosidade(BFadREAL po, BFadREAL &poros){
  const double comp = 3.625943e-10;
  const double pref = 101352.93;
  const double porosRef = 0.22;
  poros = porosRef*exp(comp*(po-pref));  
}

///Dados constantes

/** Densidade do oleo em condicoes padroes - kg/m3
  */
double TPZBlackOil2P3D::RhoOleoSC(){
  return 740.75782;
}

/** Densidade da agua em condicoes padroes - kg/m3
  */
double TPZBlackOil2P3D::RhoAguaSC(){
  return 996.95712;
}

/** Aceleracao da gravidade
  */
double TPZBlackOil2P3D::g(){
  return 9.81;
}

/** Bw = constante
  */
double TPZBlackOil2P3D::Bw(){
  return 1.041;
}

/** Viscosidade da agua (constante)
*/
double TPZBlackOil2P3D::ViscAgua(){
  return 0.31e-3;
}

/** Permeabilidade absoluta 
  */
void TPZBlackOil2P3D::K(TPZFMatrix &K){
  K.Resize(3,3);
  K.Zero();
  K(0,0) = 2.96077e-10;
  K(1,1) = 2.96077e-10;
  K(2,2) = 4.93462e-11;
}


///Programa

TPZBlackOil2P3D::TPZBlackOil2P3D(int id, double deltaT):TPZDiscontinuousGalerkin(id){
  this->fDeltaT = deltaT;
}

TPZBlackOil2P3D::~TPZBlackOil2P3D(){
///nothing to be done
}

TPZBlackOil2P3D::TPZBlackOil2P3D(const TPZBlackOil2P3D &cp):TPZDiscontinuousGalerkin(cp){

}

TPZAutoPointer<TPZMaterial> TPZBlackOil2P3D::NewMaterial(){
  return new TPZBlackOil2P3D(*this);
}

void TPZBlackOil2P3D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){

  ///un ou un+1
  double stateVal = 0.;
  if(gState == ELastState) stateVal = -1.;
  if(gState == ECurrentState) stateVal = +1.;

  ///pressao e saturacao
  const BFadREAL po(data.sol[0],0);
  const BFadREAL So(data.sol[1],1);
  BFadREAL pc;
  this->PressaoCapilar(So,pc);
  const BFadREAL pw = po - pc;
  const BFadREAL Sw = 1.-So;

  ///porosidade
  BFadREAL porosidade;
  this->Porosidade(po,porosidade);

  ///fator volume formacao
  BFadREAL Bo;
  this->Bo(po,Bo);
  const double Bw = this->Bw();

  ///Equacao 1
  BFadREAL VolOp1 = (porosidade*So/Bo)/this->fDeltaT;

  ///Equacao 2
  BFadREAL VolOp2 = (porosidade*Sw/Bw)/this->fDeltaT;

  ///ef = R = -1.( (f(un+1) - f(un))/dT - Fluxos )
  ef(0,0) += -1.*weight*stateVal*VolOp1.val();
  ef(1,0) += -1.*weight*stateVal*VolOp2.val();

  ///ek = -T (R)
  ek(0,0) += +1.*weight*stateVal*VolOp1.dx(0);
  ek(0,1) += +1.*weight*stateVal*VolOp1.dx(1);

  ek(1,0) += +1.*weight*stateVal*VolOp2.dx(0);
  ek(1,1) += +1.*weight*stateVal*VolOp2.dx(1);

}///method

void TPZBlackOil2P3D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc){
  cout << "Error: This method shoud not be called. " << __PRETTY_FUNCTION__ << "\n";
}///method

void TPZBlackOil2P3D::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){

  if(gState == ELastState) return;

  ///calculando distancia entre centro dos elementos
  double dist = 0;
  for(int i = 0; i < 3; i++){
    double val = data.XRightElCenter[i] - data.XLeftElCenter[i];
    dist += val*val;
  }
  dist = sqrt(dist);

  ///pressao e saturacao
  const BFadREAL poL(data.soll[0],0);
  const BFadREAL SoL(data.soll[1],1);
  const BFadREAL poR(data.solr[0],2);
  const BFadREAL SoR(data.solr[1],3);

  BFadREAL pcL,pcR;
  this->PressaoCapilar(SoL,pcL);
  this->PressaoCapilar(SoR,pcR);
  const BFadREAL pwL = poL - pcL;
  const BFadREAL pwR = poR - pcR;
  const BFadREAL SwL = 1.-SoL;
  const BFadREAL SwR = 1.-SoR;

  ///Permeabilidade
  TPZFNMatrix<9> K(3,3,0.);
  this->K(K);

  const double knormal = K(0,0)*data.normal[0]*data.normal[0] + 2.*K(0,1)*data.normal[0]*data.normal[1] +
                         K(1,1)*data.normal[1]*data.normal[1] + data.normal[2]*(2.*K(0,2)*data.normal[0] + 
                         2.*K(1,2)*data.normal[1] + K(2,2)*data.normal[2]);

  const double kgradZn = -K(0,2)*data.normal[0] - K(1,2)*data.normal[1] - K(2,2)*data.normal[2];

/// ************* Equacao 1 ******************* /

  ///Oleo
  BFadREAL BoL,BoR;
  this->Bo(poL, BoL);
  this->Bo(poR, BoR);
  const BFadREAL GammaOleoLeft  = this->g() * this->RhoOleoSC()/BoL;
  const BFadREAL GammaOleoRight = this->g() * this->RhoOleoSC()/BoR;

  ///velocidade de Darcy
  BFadREAL velocOleo = (knormal*poR-knormal*poL)/dist - (GammaOleoRight*kgradZn+GammaOleoLeft*kgradZn)/2.;

  ///Mobilidades
  BFadREAL KroL,KroR;
  this->Kro(SoL,KroL);
  this->Kro(SoR,KroR);
  BFadREAL ViscOleoLeft, ViscOleoRight;
  this->ViscOleo(poL, ViscOleoLeft);
  this->ViscOleo(poR, ViscOleoRight);
  BFadREAL LambdaOleoLeft = KroL/(ViscOleoLeft*BoL);
  BFadREAL LambdaOleoRight = KroR/(ViscOleoRight*BoR);

  ///Fluxo numerico da primeira equacao do residuo
  BFadREAL Fn1 = 0.;
  if(velocOleo < 0.){
    Fn1 = LambdaOleoLeft*velocOleo;
  }
  else{
    Fn1 = LambdaOleoRight*velocOleo;
  }

  ///ef = R = -1.( (f(un+1) - f(un))/dT - Fluxos )
  ef(0,0) += +1.*weight*( +Fn1.val() );
  ef(2,0) += +1.*weight*( -Fn1.val() );

  ///ek = -T (R)
  ek(0,0) += -1.*weight*( +Fn1.dx(0) );
  ek(0,1) += -1.*weight*( +Fn1.dx(1) );
  ek(0,2) += -1.*weight*( +Fn1.dx(2) );
  ek(0,3) += -1.*weight*( +Fn1.dx(3) );

  ek(2,0) += -1.*weight*( -Fn1.dx(0) );
  ek(2,1) += -1.*weight*( -Fn1.dx(1) );
  ek(2,2) += -1.*weight*( -Fn1.dx(2) );
  ek(2,3) += -1.*weight*( -Fn1.dx(3) );


/// ************* Equacao 2 ******************* /

  ///Agua
  const double Bw = this->Bw();
  const double GammaAgua = this->g() * this->RhoAguaSC() / Bw;

  ///velocidade de Darcy
  BFadREAL velocAgua = (knormal*pwR-knormal*pwL)/dist - (GammaAgua*kgradZn+GammaAgua*kgradZn)/2.;

  ///Mobilidades
  BFadREAL KrwL,KrwR;
  this->Krw(SoL,KrwL);
  this->Krw(SoR,KrwR);
  double ViscAgua = this->ViscAgua();
  BFadREAL LambdaAguaLeft = KrwL/(ViscAgua*Bw);
  BFadREAL LambdaAguaRight = KrwR/(ViscAgua*Bw);

  ///Fluxo numerico da segunda equacao do residuo
  BFadREAL Fn2 = 0.;
  if(velocAgua < 0.){
    Fn2 = LambdaAguaLeft*velocAgua;
  }
  else{
    Fn2 = LambdaAguaRight*velocAgua;
  }

  ///ef = R = -1.( (f(un+1) - f(un))/dT - Fluxos )
  ef(1,0) += +1.*weight*( +Fn2.val() );
  ef(3,0) += +1.*weight*( -Fn2.val() );

  ///ek = -T (R)
  ek(1,0) += -1.*weight*( +Fn2.dx(0) );
  ek(1,1) += -1.*weight*( +Fn2.dx(1) );
  ek(1,2) += -1.*weight*( +Fn2.dx(2) );
  ek(1,3) += -1.*weight*( +Fn2.dx(3) ); 

  ek(3,0) += -1.*weight*( -Fn2.dx(0) );
  ek(3,1) += -1.*weight*( -Fn2.dx(1) );
  ek(3,2) += -1.*weight*( -Fn2.dx(2) );
  ek(3,3) += -1.*weight*( -Fn2.dx(3) ); 

}///method

void TPZBlackOil2P3D::ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){

  if(gState == ELastState) return;

  if(bc.Type() == 2) return;///Parede ou simetria

  if((bc.Type() == 0) || (bc.Type() == 1)){

    ///pressao e saturacao
    const BFadREAL poL(data.soll[0],0);
    const BFadREAL SoL(data.soll[1],1);

    data.solr[0] = bc.Val2()(0,0);
    data.solr[1] = bc.Val2()(1,0);
    this->ContributeInterface(data,weight,ek,ef);

  }///Dirichlet na pressao e Dirichlet ou outflow na saturacao


/*  switch(bc.Type()) {
  case 0: /// DIRICHLET na pressao e na saturacao


    break;
  case 1: /// DIRICHLET na pressao e OUTFLOW na saturacao

    break;

  case 2: /// PAREDE / SIMETRIA

    break;

  default:
    PZError << __PRETTY_FUNCTION__ << " - Wrong boundary condition type\n";
    break;
  }*/

}///method


enum ESolutionVars { ENone = 0, EWaterPressure = 1, EOilPressure, EWaterSaturation, EOilSaturation, EDarcyVelocity };

int TPZBlackOil2P3D::VariableIndex(const std::string &name){
  if(!strcmp("WaterPressure",name.c_str()))   return  EWaterPressure;
  if(!strcmp("OilPressure",name.c_str()))     return  EOilPressure;
  if(!strcmp("WaterSaturation",name.c_str())) return  EWaterSaturation;
  if(!strcmp("OilSaturation",name.c_str()))   return  EOilSaturation;
  if(!strcmp("DarcyVelocity",name.c_str()))   return  EDarcyVelocity;
  return TPZMaterial::VariableIndex(name);
}

int TPZBlackOil2P3D::NSolutionVariables(int var){
  if(var == EWaterPressure) return 1;
  if(var == EOilPressure) return 1;
  if(var == EWaterSaturation) return 1;
  if(var == EOilSaturation) return 1;
  if(var == EDarcyVelocity) return 3;
}

void TPZBlackOil2P3D::Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
                      TPZFMatrix &axes, int var, TPZVec<REAL> &Solout){

///please implement me
}

#endif
