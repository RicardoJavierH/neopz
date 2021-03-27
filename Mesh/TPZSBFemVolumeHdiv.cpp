
//
//  TPZSBFemVolumeHdiv.cpp
//  PZ
//
//  Created by Karolinne Coelho on 25/01/2021.
//
//

#include "TPZSBFemVolumeHdiv.h"
#include "pzgeoelside.h"
#include "TPZSBFemElementGroup.h"
#include "pzintel.h"
#include "TPZMaterial.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzcmesh.h"
#include "TPZGeoLinear.h"
#include "pzmultiphysicscompel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.sbfemvolume"));
#endif

TPZSBFemVolumeHdiv::TPZSBFemVolumeHdiv(TPZMultiphysicsCompMesh & mesh, TPZGeoEl * gel, int64_t & index) : TPZInterpolationSpace(mesh, gel, index)
{
    fElementVec.Resize(7);
}

void TPZSBFemVolumeHdiv::SetElementGroupIndex(int64_t index)
{
    fElementGroupIndex = index;
    TPZCompEl *celgr = Mesh()->Element(index);
    fElementGroup = celgr;
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi.
 * @param qsi master element coordinate
 * @param sol finite element solution
 * @param dsol solution derivatives
 * @param axes axes associated with the derivative of the solution
 */
void TPZSBFemVolumeHdiv::ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol, TPZFMatrix<REAL> &axes)
{
    TPZCompMesh *cmesh = Mesh();
    sol.Resize(fCoeficients.Cols());
    dsol.Resize(fCoeficients.Cols());
    TPZGeoEl *Ref2D = Reference();
    int matid = Ref2D->MaterialId();
    TPZMaterial *mat2d = cmesh->FindMaterial(matid);

    int dim = Ref2D->Dimension();
    REAL sbfemparam = (1. - qsi[dim - 1]) / 2.;
    if (sbfemparam < 0.) {
        std::cout << "sbfemparam " << sbfemparam << std::endl;
        sbfemparam = 0.;
    }
    if (IsZero(sbfemparam)) {
        for (int i = 0; i < dim - 1; i++) {
            qsi[i] = 0.;
        }
        if (dim == 2) {
            sbfemparam = 1.e-6;
            qsi[dim - 1] = 1. - 2.e-6;
        } else {
            sbfemparam = 1.e-4;
            qsi[dim - 1] = 1. - 2.e-4;
        }
    }
    auto CSkeleton = dynamic_cast<TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear> *> (fElementVec[6]);
    TPZMaterialData data2d;
    // compute the lower dimensional shape functions
    TPZManVector<REAL, 3> qsilow(qsi);
    qsilow.Resize(dim - 1);

    TPZManVector<TPZMaterialData,2> datavec(2);
    CSkeleton->InitMaterialData(datavec);
    auto data1d = datavec[1];
    TPZGeoEl *Ref1D = CSkeleton->Reference();

    Ref1D->Jacobian(qsilow, data1d.jacobian, data1d.axes, data1d.detjac, data1d.jacinv);
    Ref2D->Jacobian(qsi, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
    // if (dim == 3) {

    //     AdjustAxes3D(data1d.axes, data2d.axes, data2d.jacobian, data2d.jacinv, data2d.detjac);
    // }
    axes = data2d.axes;

    TPZVec<TPZTransform<> > trvec;
    CSkeleton->AffineTransform(trvec);
    CSkeleton->ComputeRequiredData(qsilow, trvec, datavec);

    data1d = datavec[1];

    int nshape = data1d.phi.Rows();
    int nstate = mat2d->NStateVariables();
#ifdef PZDEBUG
    if (fPhi.Cols() != fCoeficients.Rows()) {
        DebugStop();
    }
#endif
#ifdef LOG4CXX2
    if (logger->isDebugEnabled()) {
        TPZManVector<std::complex<double> > coefcol(fCoeficients.Rows());
        for (int i = 0; i < fCoeficients.Rows(); i++) {
            coefcol[i] = fCoeficients(i, 0);
        }
        std::stringstream sout;
        sout << "coefficients " << coefcol << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    for (int s = 0; s < sol.size(); s++) {
        TPZManVector<std::complex<double>, 10> uh_xi(fPhi.Rows(), 0.), Duh_xi(fPhi.Rows(), 0.);
        int nphixi = fPhi.Rows();
        int numeig = fPhi.Cols();
        for (int c = 0; c < numeig; c++) {
            std::complex<double> xiexp;
            std::complex<double> xiexpm1;
            if (IsZero(fEigenvalues[c] + 0.5 * (dim - 2))) {
                xiexp = 1;
                xiexpm1 = 0;
            } else if (IsZero(fEigenvalues[c] + 1. + 0.5 * (dim - 2))) {
                xiexp = sbfemparam;
                xiexpm1 = 1;
            } else {
                xiexp = pow(sbfemparam, -fEigenvalues[c] - 0.5 * (dim - 2));
                xiexpm1 = pow(sbfemparam, -fEigenvalues[c] - 1. - 0.5 * (dim - 2));
            }
            for (int i = 0; i < nphixi; i++) {
                uh_xi[i] += fCoeficients(c, s) * xiexp * fPhi(i, c);
                Duh_xi[i] += -fCoeficients(c, s)*(fEigenvalues[c] + 0.5 * (dim - 2)) * xiexpm1 * fPhi(i, c);
            }
        }
#ifdef LOG4CXX2
        if (s == 0 && logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "uh_xi " << uh_xi << std::endl;
            sout << "Duh_xi " << Duh_xi << std::endl;
            data1d.phi.Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        //        std::cout << "uh_xi " << uh_xi << std::endl;
        //        std::cout << "Duh_xi " << Duh_xi << std::endl;
        sol[s].Resize(nstate);
        sol[s].Fill(0.);
        TPZFNMatrix<9, STATE> dsollow(dim - 1, nstate, 0.), dsolxieta(dim, nstate, 0.);
        TPZManVector<STATE, 3> dsolxi(nstate, 0.);
        for (int ishape = 0; ishape < nshape; ishape++) {
            for (int istate = 0; istate < nstate; istate++) {
                sol[s][istate] += data1d.phi(ishape) * uh_xi[ishape * nstate + istate].real();
                dsolxi[istate] += data1d.phi(ishape) * Duh_xi[ishape * nstate + istate].real();
                for (int d = 0; d < dim - 1; d++) {
                    dsollow(d, istate) += data1d.dphi(d, ishape) * uh_xi[ishape * nstate + istate].real();
                }
            }
        }
        for (int istate = 0; istate < nstate; istate++) {
            for (int d = 0; d < dim - 1; d++) {
                dsolxieta(d, istate) = dsollow(d, istate);
            }
            dsolxieta(dim - 1, istate) = -dsolxi[istate] / 2.;
        }
        dsol[s].Resize(dim, nstate);
        dsol[s].Zero();
        for (int istate = 0; istate < nstate; istate++) {
            for (int d1 = 0; d1 < dim; d1++) {
                for (int d2 = 0; d2 < dim; d2++) {
                    dsol[s](d1, istate) += data2d.jacinv(d2, d1) * dsolxieta(d2, istate);
                }
            }
        }
    }
}

void TPZSBFemVolumeHdiv::Solution(TPZManVector<REAL> &qsi,int var,TPZManVector<STATE> &sol)
{
    TPZGeoEl *Ref2D = this->Reference();
    int matid = Ref2D->MaterialId();
    TPZCompMesh *cmesh = Mesh();

    TPZMaterial *mat2d = cmesh->FindMaterial(matid);
    TPZMaterialData data2d;
    
    ComputeSolution(qsi, data2d.sol, data2d.dsol, data2d.axes);
    
    data2d.x.Resize(3, 0.);
    Reference()->X(qsi, data2d.x);
    mat2d->Solution(data2d, var, sol);
}

TPZCompEl * CreateSBFemMultiphysicsCompEl(TPZMultiphysicsCompMesh &mesh, TPZGeoEl *gel, int64_t &index)
{
    new TPZSBFemVolumeHdiv(mesh, gel, index);    
}


void TPZSBFemVolumeHdiv::Shape(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphidxi)
{
    TPZCompMesh *cmesh = Mesh();
    TPZCompEl *celgroup = cmesh->Element(fElementGroupIndex);
    auto elgr = dynamic_cast<TPZSBFemMultiphysicsElGroup *> (celgroup);
    TPZFMatrix<std::complex<double> > &CoefficientLoc = elgr->PhiInverse();
#ifdef LOG4CXX2
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        CoefficientLoc.Print("Coefficients = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZGeoEl *Ref2D = Reference();
    int matid = Ref2D->MaterialId();
    TPZMaterial *mat2d = cmesh->FindMaterial(matid);
    int dim = Ref2D->Dimension();
    int nstate = mat2d->NStateVariables();

    phi.Redim(CoefficientLoc.Cols() * nstate, 1);
    dphidxi.Redim(dim*nstate, CoefficientLoc.Cols());

    REAL sbfemparam = (1. - qsi[dim - 1]) / 2.;
    if (sbfemparam < 0.) {
        std::cout << "sbfemparam " << sbfemparam << std::endl;
        sbfemparam = 0.;
    }
    if (IsZero(sbfemparam)) {
        for (int i = 0; i < dim - 1; i++) {
            qsi[i] = 0.;
        }
        if (dim == 2) {
            sbfemparam = 1.e-6;
            qsi[dim - 1] = 1. - 2.e-6;
        } else {
            sbfemparam = 1.e-4;
            qsi[dim - 1] = 1. - 2.e-4;
        }
    }
    auto CSkeleton = dynamic_cast<TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear> *> (fElementVec[6]);
    TPZMaterialData data2d;
    // compute the lower dimensional shape functions
    TPZManVector<REAL, 3> qsilow(qsi);
    qsilow.Resize(dim - 1);

    TPZManVector<TPZMaterialData,2> datavec(2);
    CSkeleton->InitMaterialData(datavec);

    auto data1d = datavec[1];
    TPZGeoEl *Ref1D = CSkeleton->Reference();

    Ref1D->Jacobian(qsilow, data1d.jacobian, data1d.axes, data1d.detjac, data1d.jacinv);
    Ref2D->Jacobian(qsi, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
    if (dim == 3) {

        // AdjustAxes3D(data1d.axes, data2d.axes, data2d.jacobian, data2d.jacinv, data2d.detjac);
    }
    TPZVec<TPZTransform<> > trvec;
    CSkeleton->AffineTransform(trvec);
    CSkeleton->ComputeRequiredData(qsilow, trvec, datavec);

    data1d = datavec[1];

    int nshape = data1d.phi.Rows();
#ifdef PZDEBUG
    if (fPhi.Cols() != fCoeficients.Rows()) {
        DebugStop();
    }
#endif
    phi.Zero();
#ifdef LOG4CXX2
    if (logger->isDebugEnabled()) {
        int eq = 1;
        TPZManVector<std::complex<double> > coefcol(CoefficientLoc.Rows());
        for (int i = 0; i < CoefficientLoc.Rows(); i++) {
            coefcol[i] = CoefficientLoc(i, eq);
        }
        std::stringstream sout;
        sout << "coefficients " << coefcol << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    for (int s = 0; s < CoefficientLoc.Cols(); s++) {
        TPZManVector<std::complex<double>, 10> uh_xi(fPhi.Rows(), 0.), Duh_xi(fPhi.Rows(), 0.);
        int nphixi = fPhi.Rows();
        int numeig = fPhi.Cols();
        for (int c = 0; c < numeig; c++) {
            std::complex<double> xiexp;
            std::complex<double> xiexpm1;
            if (IsZero(fEigenvalues[c] + 0.5 * (dim - 2))) {
                xiexp = 1;
                xiexpm1 = 0;
            } else if (IsZero(fEigenvalues[c] + 1. + 0.5 * (dim - 2))) {
                xiexp = sbfemparam;
                xiexpm1 = 1;
            } else {
                xiexp = pow(sbfemparam, -fEigenvalues[c] - 0.5 * (dim - 2));
                xiexpm1 = pow(sbfemparam, -fEigenvalues[c] - 1. - 0.5 * (dim - 2));
            }
            for (int i = 0; i < nphixi; i++) {
                uh_xi[i] += CoefficientLoc(c, s) * xiexp * fPhi(i, c);
                Duh_xi[i] += -CoefficientLoc(c, s)*(fEigenvalues[c] + 0.5 * (dim - 2)) * xiexpm1 * fPhi(i, c);
            }
        }
#ifdef LOG4CXX2
        if (s == 1 && logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "uh_xi " << uh_xi << std::endl;
            sout << "Duh_xi " << Duh_xi << std::endl;
            data1d.phi.Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        //        std::cout << "uh_xi " << uh_xi << std::endl;
        //        std::cout << "Duh_xi " << Duh_xi << std::endl;
        TPZFNMatrix<9, STATE> dsollow(dim - 1, nstate, 0.), dsolxieta(dim, nstate, 0.);
        TPZManVector<STATE, 3> dsolxi(nstate, 0.);
        for (int ishape = 0; ishape < nshape; ishape++) {
            for (int istate = 0; istate < nstate; istate++) {
                phi(s * nstate + istate, 0) += data1d.phi(ishape) * uh_xi[ishape * nstate + istate].real();
                //                sol[s][istate] += data1d.phi(ishape)*uh_xi[ishape*nstate+istate].real();
                dsolxi[istate] += data1d.phi(ishape) * Duh_xi[ishape * nstate + istate].real();
                for (int d = 0; d < dim - 1; d++) {
                    dsollow(d, istate) += data1d.dphi(d, ishape) * uh_xi[ishape * nstate + istate].real();
                }
            }
        }
        for (int istate = 0; istate < nstate; istate++) {
            for (int d = 0; d < dim - 1; d++) {
                dsolxieta(d, istate) = dsollow(d, istate);
            }
            dsolxieta(dim - 1, istate) = -dsolxi[istate] / 2.;
        }
        for (int istate = 0; istate < nstate; istate++) {
            for (int d1 = 0; d1 < dim; d1++) {
                for (int d2 = 0; d2 < dim; d2++) {
                    //                    dsol[s](d1,istate) += data2d.jacinv(d2,d1)*dsolxieta(d2,istate);
                    dphidxi(istate * nstate + d1, s) += data2d.jacinv(d2, d1) * dsolxieta(d2, istate);
                }
            }
        }
    }
}