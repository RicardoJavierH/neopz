//
//  TPZReynoldsFlow.cpp
//  PZ
//
//  Created by Cesar Lucci on 06/02/13.
//
//

#include "TPZReynoldsFlow.h"
#include "pzbndcond.h"

TPZReynoldsFlow::TPZReynoldsFlow() : TPZMaterial()
{
    f_visc = 0.;
    f_deltaT = 0.;
    f_staticPotential = 0.;
    f_nplus1Computation = false;
}

TPZReynoldsFlow::TPZReynoldsFlow(int matId, REAL visc, REAL deltaT, REAL staticPotential) : TPZMaterial(matId)
{
    f_visc = visc;
    f_deltaT = deltaT;
    f_staticPotential = staticPotential;
    f_nplus1Computation = false;
}

TPZReynoldsFlow::TPZReynoldsFlow(const TPZReynoldsFlow &cp) : TPZMaterial(cp)
{
    f_visc = cp.f_visc;
    f_deltaT = cp.f_deltaT;
    f_staticPotential = cp.f_staticPotential;
    f_nplus1Computation = cp.f_nplus1Computation;
}


TPZReynoldsFlow::~TPZReynoldsFlow()
{

}


int TPZReynoldsFlow::Dimension()
{
    return 2;
}


int TPZReynoldsFlow::NStateVariables()
{
    return 1;
}

void TPZReynoldsFlow::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if(f_nplus1Computation == false)//estamos no passo n
    {
        REAL simmetryy = 2.;
        REAL wn = simmetryy * data.sol[1][1];//data.sol = {{p},{ux,uy,uz}} e estamos interessados em uy
        for(int i = 0; i < data.phi.Rows(); i++)
        {
            ef(i,0) += weight * (wn/f_deltaT) * data.phi(i,0);
        }
    }
    else//estamos no passo n+1
    {
        REAL simmetryy = 2.;
        REAL wnplus1 = simmetryy * data.sol[1][1];//data.sol = {{p},{ux,uy,uz}} e estamos interessados em uy
        REAL w3 = wnplus1*wnplus1*wnplus1;
        REAL carterGAMMA = 1.;//TODO : Outra cmesh ou estrutura de dados? Lembre-se que valria com o tempo e no espaco!!!
        for(int i = 0; i < data.phi.Rows(); i++)
        {
            ef(i,0) += weight * (carterGAMMA * f_staticPotential - wnplus1/f_deltaT)*data.phi(i,0);
        }
        for(int i = 0; i < data.phi.Rows(); i++)
        {
            for(int j = 0; j < data.phi.Rows(); j++)
            {
                ek(i,j) += weight * (w3/(12.*f_visc) * (data.dphix(0,j)*data.dphix(0,i) + data.dphix(1,j)*data.dphix(1,i)));
                ek(i,j) += weight * (carterGAMMA*data.phi(j,0)*data.phi(i,0));
            }
        }
    }
}


void TPZReynoldsFlow::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    if(f_nplus1Computation == false)
    {
        return;
    }
    
    if(bc.Type() != 1)
    {
        DebugStop();
    }
    
    for(int i = 0; i < data.phi.Rows(); i++)
    {
        ef(i,0) += weight * bc.Val2()(0,0) * data.phi(i,0);
    }
}


TPZMaterial * TPZReynoldsFlow::NewMaterial()
{
    return new TPZReynoldsFlow(*this);
}

int TPZReynoldsFlow::ClassId() const
{
    DebugStop();
}

void TPZReynoldsFlow::Write(TPZStream &buf, int withclassid)
{
    DebugStop();
}

void TPZReynoldsFlow::Read(TPZStream &buf, void *context)
{
    DebugStop();
}
