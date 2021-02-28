//
//  TPZSBFemVolumeHdiv.hpp
//  PZ
//
//  Created by Karolinne Coelho on 25/01/2021.
//
//

#pragma once

#include <stdio.h>
#include "pzcompel.h"
#include "pzelmat.h"
#include "TPZSBFemVolume.h"
#include "TPZCompElHDivSBFem.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSBFemMultiphysicsElGroup.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzgeoquad.h"

#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzelchdiv.h"

using namespace std;
using namespace pzshape;

class TPZSBFemVolumeHdiv : public TPZInterpolationSpace
{
    /// index of element group
    int64_t fElementGroupIndex = -1;

    int fSkeleton = -1;

    TPZCompEl *fElementGroup = 0;

	/** @brief List of pointers to computational elements */
	TPZManVector<TPZCompEl* ,7> fElementVec;
    	
	/** @brief Indexes of the connects of the element */
	TPZVec<int64_t> fConnectIndexes;
    
    /// pointer to the integration rule
    TPZIntPoints *fIntRule = 0;
    
    /// vector of local indices of multipliers in the group
    TPZManVector<int64_t> fLocalIndices;
    
    /// Section of the phi vector associated with this volume element
    TPZFNMatrix<30,std::complex<double> > fPhi;
    
    /// Eigenvlues associated with the internal shape functions
    TPZManVector<std::complex<double> > fEigenvalues;
    
    /// Multiplier coeficients associated with the solution
    TPZFNMatrix<30,std::complex<double> > fCoeficients;

public:
    
    TPZSBFemVolumeHdiv(TPZMultiphysicsCompMesh & mesh, TPZGeoEl * gel, int64_t & index);
    
    virtual ~TPZSBFemVolumeHdiv()
    {
        // Reference()->ResetReference();
    }

    void AddElement(TPZCompEl * cel, int localindex)
    {
        fElementVec[localindex] = cel;
        auto ncon = fConnectIndexes.size();
        auto nconcel =cel->NConnects();
        fConnectIndexes.Resize(ncon+nconcel);
        for (auto i = 0; i < nconcel; i++)
        {
            fConnectIndexes[i+ncon] = cel->ConnectIndex(i);
        }
    }

    void SetSkeleton(int64_t skeleton)
    {
        fSkeleton = skeleton;
    }
    
    int64_t SkeletonIndex()
    {
        return fSkeleton;
    }

    void SetPressureIds(int64_t leftpressure, int64_t rightpressure)
    {
        DebugStop();
    }

    void SetFluxIds(int64_t leftflux, int64_t rightflux)
    {
        DebugStop();
    }

    /** @brief Method for creating a copy of the element */
    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const
    {
        DebugStop();
        return 0;
    }

    virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                                    std::map<int64_t,int64_t> & gl2lcConMap,
                                    std::map<int64_t,int64_t> & gl2lcElMap) const
    {
        DebugStop();
        return 0;
    }

    void SetElementGroupIndex(int64_t index);

    /** @brief Returns the number of nodes of the element */
    virtual int NConnects() const
    {
        return fConnectIndexes.size();   
    }

    virtual int64_t ConnectIndex(int i) const
    {
        if (i > fConnectIndexes.size())
        {
            DebugStop();
        }        
        return fConnectIndexes[i];
    }

    virtual int Dimension() const
    {
        TPZGeoEl *reference = Reference();
        return reference->Dimension();
    }

    virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef) override
    {
        DebugStop();
    }

    virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) override
    {
        DebugStop();
    }

    virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override
    {
        DebugStop();
    }

    virtual void SetConnectIndex(int inode, int64_t index)
    {
        DebugStop();
    }

    virtual TPZIntPoints & GetIntegrationRule() const
    {
        if(!fIntRule) DebugStop();
		return *fIntRule;
    }

    virtual TPZIntPoints & GetIntegrationRule()
    {
        if(!fIntRule) InitializeIntegrationRule();
        return *fIntRule;
    }

    void InitializeIntegrationRule()
    {
        if (fIntRule) {
            DebugStop();
        }
        int nsides = Reference()->NSides();
        fIntRule = Reference()->CreateSideIntegrationRule(nsides-1, 1);
    }

    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const
    {
        DebugStop();
    }

    virtual int NSideConnects(int iside) const
    {
        DebugStop();
        return 0;
    }

    virtual int SideConnectLocId(int icon,int is) const
    {
        DebugStop();
        return 0;
    }

    virtual int NShapeF() const
    {
        DebugStop();
        return 0;
    }

    virtual int NConnectShapeF(int icon, int order) const
    {
        DebugStop();
        return 0;
    }

    virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphidxi)
    {
        DebugStop();
    }

    virtual void SetPreferredOrder ( int order )
    {
        DebugStop();
    }

    virtual void PRefine ( int order )
    {
        DebugStop();
    }

    virtual TPZCompEl *Element(int elindex)
    {
        return fElementVec[elindex];
    }

    virtual TPZManVector<TPZCompEl * ,7> & ElementVec()
    {
        return fElementVec;
    }

    virtual void SetConnectIndexes(TPZVec<int64_t> &indexes)
    {
        fConnectIndexes = indexes;
    }

    // virtual void AffineTransform(TPZVec<TPZTransform<> > &trVec) const
    // {
    //     DebugStop();
    // }

    // virtual void InitMaterialData(TPZVec<TPZMaterialData > &dataVec, TPZVec<int64_t> *indices = 0)
    // {
    //     DebugStop();
    // }

    // virtual void PolynomialOrder(TPZVec<int> &order) const
    // {
    //     DebugStop();
    // }


    // virtual void Print(std::ostream &out = std::cout) const
    // {
    //     out << "Printing " << __PRETTY_FUNCTION__ << std::endl;
    //     TPZCompEl::Print(out);
    //     out << "Group Element Index " << fElementGroupIndex << std::endl;
    //     out << "Skeleton Element Index " << fSkeleton << std::endl;
    //     out << "Local Indices " << fLocalIndices << std::endl;
    //     fCoeficients.Print("Coef =",out,EMathematicaInput);
    //     fPhi.Print("Phi = ",out,EMathematicaInput);
    //     if (fCoeficients.Rows())
    //     {
    //         TPZManVector<std::complex<double>,5> prod(fPhi.Rows(),0.);
    //         for (int i=0; i<fPhi.Rows(); i++) {
    //             for (int j=0; j<fPhi.Cols(); j++) {
    //                 prod[i] += fPhi.GetVal(i,j)*fCoeficients.GetVal(j,0);
    //             }
    //         }
    //         out << "Values at border " << prod << std::endl;
    //     }
    //     for (int i=0; i<NConnects(); i++) {
    //         Connect(fConnectIndexes[i]).Print(*Mesh(),out);
    //     }
    // }
};

TPZCompEl * CreateSBFemMultiphysicsCompEl(TPZMultiphysicsCompMesh &mesh, TPZGeoEl *gel, int64_t &index);