// -*- c++ -*-

//$Id: meshes.h,v 1.3 2006-01-18 15:22:12 tiago Exp $

class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;

#include "pzvec.h"

void SetPOrder(int p);

//Cria malha simple com solucao analitica polinomial*E^x para validar TPZCoupledTransportDarcy
TPZCompMesh * CreateSimpleMeshWithExactSolution2(int h, int p);
void ExactSol_u2(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv);
void Forcing22(TPZVec<REAL> &pt, TPZVec<REAL> &force);

//Cria malha simple com solucao analitica polinomial para validar TPZCoupledTransportDarcy
TPZCompMesh * CreateSimpleMeshWithExactSolution(int h, int p);
void ExactSol_p(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv);
void ExactSol_u(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv);
void Forcing2(TPZVec<REAL> &pt, TPZVec<REAL> &force);
void Dirichlet1(TPZVec<REAL> &pt, TPZVec<REAL> &force);

//**Para verificar material Poisson com beta variavel. Eh preciso programar beta no TPZMatPoisson3d::Contribute
TPZCompMesh * CheckBetaNonConstant(int h, int p);
void SolExata(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv);
void Force(TPZVec<REAL> &pt, TPZVec<REAL> &force);
