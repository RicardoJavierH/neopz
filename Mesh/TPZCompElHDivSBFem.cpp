/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDivSBFem methods.
 */


#include "TPZCompElHDivSBFem.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzlog.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"
#include "pzelchdiv.h"
#include "pzvec_extras.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDivSBFem"));
#endif

// Initialize with the geometry of the SBFemVolume
template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem(TPZCompMesh &mesh, TPZGeoEl *gel, TPZGeoElSide &gelside, int64_t &index) :
TPZCompElHDivCollapsed<TSHAPE>(mesh, gel, index), fGelVolSide(gelside)
{
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZCompElHDivCollapsed<TSHAPE>(mesh, gel, index), fGelVolSide(0)
{
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem(TPZCompMesh &mesh, const TPZCompElHDivSBFem<TSHAPE> &copy) :
TPZCompElHDivCollapsed<TSHAPE>(mesh,copy)
{
    fGelVolSide = copy.fGelVolSide;
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem() : TPZCompElHDivCollapsed<TSHAPE>()
{
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::~TPZCompElHDivSBFem()
{
}

// Set the GeoElSide for the volumetric SBFEM element
template<class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::SetGelVolumeSide(TPZGeoElSide &gelside)
{
    fGelVolSide = gelside;
}

// Get the GeoElSide for the volumetric SBFEM element
template<class TSHAPE>
TPZGeoElSide & TPZCompElHDivSBFem<TSHAPE>::GetGelVolumeSide()
{
    return fGelVolSide;
}

template<class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi)
{
    // Compute data for the 1d functions first
    bool needssol = data.fNeedsSol;
    data.fNeedsSol = false;
    TPZCompElHDivCollapsed<TSHAPE>::ComputeRequiredData(data, qsi);
    data.fNeedsSol = needssol;

    auto nshape2d = data.phi.Rows();
    auto nshape1d = 0;
    for (int i = 0; i < 3; i++)
    {
        nshape1d += this->NConnectShapeF(i, this->fPreferredOrder);
    }

    HDivCollapsedDirections(data, nshape1d);

    data.ComputeFunctionDivergence();
    
    int64_t nshape = (nshape2d - nshape1d) / 2;
    // Adjusting divergence values
    for (int64_t i = 0; i < nshape1d; i++)
    {
        data.divphi(i) = 2*data.dphi(0,i);
        data.phi(i) *= 2;
    }
    for (int64_t i = 0; i < nshape; i++)
    {
        data.divphi(i+nshape1d) = 0;
        data.divphi(i+nshape1d+nshape) = 2*data.phi(i+nshape1d);
        data.phi(i+nshape1d+nshape) = 0;
    }
    for (auto i = 0; i < nshape1d; i++)
    {
        // data.divphi(i) *= 2*data.detjac;
        // data.phi(i) *= 2*(data.detjac);
    }
    for (auto i = 0; i < nshape; i++)
    {
        // data.divphi(i) *= data.detjac;
        // data.phi(i+nshape1d) *= 1/data.detjac;
    }
    
    if (data.fNeedsSol)
    {
        std::cout << "Compute Solution not implemented yet \n";
        DebugStop();
        // ComputeSolution(qsi, data);
    }
}

// This function will compute the directions for the HDiv collapsed based on the information of neighbourhood
template<class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::HDivCollapsedDirections(TPZMaterialData &data, int64_t nshape1d)
{
    // Computing the deformed directions for the 2d functions using the information of the neighbourhood
    // Inspired in TPZSBFemVolume::ComputeKMatrices
    auto detjac1d = data.detjac;

    // The Reference element will be the skeleton
    TPZGeoEl *Ref1D = this->Reference();
    int dim1 = TSHAPE::Dimension;

    // The Volume element will be the one passed as an argument by the TPZSBFemVolumeHDiv - GEO DO EL 1D
    // The SBFemVolumeHdiv will compute the Contribute for the SBFemHdiv element
    // Before that, the material data must be computed (It will call ComputeRequiredData from TPZCompElHDivSBFem)
    TPZGeoEl * gelvolume = fGelVolSide.Element();
    int dim2 = gelvolume->Dimension();

    // Find the higher side of the skeleton element
    TPZGeoElSide SkeletonSide(Ref1D, Ref1D->NSides() - 1);

    // Create a transformation between these sides
    TPZTransform<REAL> tr(dim2, dim1);
    tr = SkeletonSide.NeighbourSideTransform(fGelVolSide);
    TPZTransform<REAL> t2 = gelvolume->SideToSideTransform(fGelVolSide.Side(), gelvolume->NSides() - 1);
    tr = t2.Multiply(tr);

    // Applying the transformation
    TPZManVector<REAL, 3> xi(dim1), xiquad(dim2), xivol(dim2);
    xi = data.xParametric;
    tr.Apply(xi, xiquad);
    xivol = xiquad;
    xivol[dim2 - 1] = -0.5;

    // Computing the data for the volumetric element
    gelvolume->X(xivol, data.x);
    gelvolume->Jacobian(xiquad, data.jacobian, data.axes, data.detjac, data.jacinv);

    auto axes = data.axes;
    if (dim2 == 3) {
        AdjustAxes3D(axes, data.axes, data.jacobian, data.jacinv, data.detjac);
    }

    // Adjusting derivatives
    ExtendShapeFunctions(data);

    auto ndir = data.fDeformedDirections.Cols()-1;

    if (dim2 == 3)
    {
        std::cout << "Function not verified for 3D cases yet\n";
        DebugStop();
    }
    
    TPZFNMatrix<9,REAL> grad(dim2,dim2,0);
    gelvolume->GradX(xivol,grad);

    // data.detjac = sqrt(2*data.detjac);
    
    for (auto j = 0; j < 3; j++)
    {
        for (auto i = 0; i < dim2; i++)
        {
            data.fDeformedDirections(i,j) = data.jacobian(i,0)/(detjac1d);
        }
        data.fDeformedDirections(2,j) = 0.;
    }
    for (auto j = 3; j < data.fDeformedDirections.Cols(); j++)
    {
        for (auto i = 0; i < dim2; i++)
        {
            // for (auto k = 0; k < dim2; k++)
            // {
                data.fDeformedDirections(i,j) = data.jacobian(i,1)/(detjac1d/2);
            // }
        }
        data.fDeformedDirections(2,j) = 0.;
    }
}

template <class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::AdjustAxes3D(const TPZFMatrix<REAL> &axes2D, TPZFMatrix<REAL> &axes3D, TPZFMatrix<REAL> &jac3D, TPZFMatrix<REAL> &jacinv3D, REAL detjac)
{
    TPZManVector<REAL, 3> ax1(3), ax2(3), ax3(3);
    for (int i = 0; i < 3; i++) {
        ax1[i] = axes2D.g(0, i);
        ax2[i] = axes2D.g(1, i);
        Cross(ax1, ax2, ax3);
    }
    for (int i = 0; i < 3; i++) {
        axes3D(0, i) = ax1[i];
        axes3D(1, i) = ax2[i];
        axes3D(2, i) = ax3[i];
        if (detjac < 0.) {
            axes3D(2, i) *= -1.;
        }
    }
    TPZFNMatrix<9, REAL> jacnew(3, 3), axest(3, 3), jacinv(3, 3);
    axes3D.Transpose(&axest);
    axes3D.Multiply(jac3D, jacnew);
    jacinv3D.Multiply(axest, jacinv);
    jac3D = jacnew;
    jacinv3D = jacinv;
#ifdef PZDEBUG
    // check whether the axes are orthogonal and whether the jacobian is still the inverse of jacinv
    {
        TPZFNMatrix<9, REAL> ident1(3, 3, 0.), ident2(3, 3, 0.), identity(3, 3);
        identity.Identity();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    ident1(i, j) += axes3D(i, k) * axes3D(j, k);
                    ident2(i, j) += jac3D(i, k) * jacinv3D(k, j);
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (fabs(ident1(i, j) - identity(i, j)) > 1.e-6) {
                    DebugStop();
                }
                if (fabs(ident2(i, j) - identity(i, j)) > 1.e-6) {
                    DebugStop();
                }
            }
        }
    }
#endif
}

template <class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::ExtendShapeFunctions(TPZMaterialData &data)
{
    int dim = fGelVolSide.Element()->Dimension();

    auto nshape = data.phi.Rows();
    auto nshape1d = 0;
    for (int i = 0; i < 3; i++)
    {
        nshape1d += this->NConnectShapeF(i, this->fPreferredOrder);
    }
    
    int64_t nshape2d = (nshape - nshape1d) / 2;
    for (int ish = 0; ish < nshape2d; ish++)
    {
        data.phi(ish + nshape1d, 0) = data.phi(ish, 0);
        data.phi(ish + nshape1d+nshape2d, 0) = data.phi(ish, 0);
        for (int d = 0; d < dim - 1; d++)
        {
            data.dphi(d, ish + nshape1d) = 0.;
            data.dphi(d, ish + nshape1d+nshape2d) = data.phi(ish);
        }
        data.dphi(dim-1, ish+nshape1d) = - data.phi(ish) / 2.;
        data.dphi(dim-1, ish+nshape1d+nshape2d) = 0.;
    }
    TPZInterpolationSpace::Convert2Axes(data.dphi, data.jacinv, data.dphix);
}

#include "pzshapetriang.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"

using namespace pzshape;

template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeQuad>>;

template class TPZCompElHDivSBFem<TPZShapeTriang>;
template class TPZCompElHDivSBFem<TPZShapeLinear>;
template class TPZCompElHDivSBFem<TPZShapeQuad>;