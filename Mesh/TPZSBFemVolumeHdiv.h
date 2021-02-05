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

#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzelchdiv.h"

using namespace std;
using namespace pzshape;

class TPZSBFemVolumeHdiv : public TPZMultiphysicsElement
{
    /// index of element group
    int64_t fElementGroupIndex = -1;

    int fSkeleton;

    TPZCompEl *fElementGroup = 0;

	/** @brief List of pointers to computational elements */
	TPZManVector<TPZCompElSide ,5> fElementVec;
    	
	/** @brief Indexes of the connects of the element */
	TPZVec<int64_t> fConnectIndexes;
    
    /// pointer to the integration rule
    TPZIntPoints *fIntRule = 0;
    
    /// vector of local indices of multipliers in the group
    TPZManVector<int64_t> fLocalIndices;

public:
    
    TPZSBFemVolumeHdiv(TPZMultiphysicsCompMesh & mesh, TPZGeoEl * gel, int64_t & index);
    
    virtual ~TPZSBFemVolumeHdiv()
    {
        // Reference()->ResetReference();
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
        if (fElementGroup == 0) {
            return 0;
        }
        return fElementGroup->NConnects();
    }

    virtual int64_t ConnectIndex(int i) const
    {
        if (fElementGroup == 0) {
            DebugStop();
        }
        return fElementGroup->ConnectIndex(i);
    }

    virtual int Dimension() const
    {
        TPZGeoEl *reference = Reference();
        return reference->Dimension();
    }

    virtual void AddElement(const TPZCompElSide &cel, int64_t mesh)
    {
        if (fElementVec.size() <= mesh)
        {
            fElementVec.resize(mesh+1);
            fActiveApproxSpace.Resize(mesh+1, 1);
        }
        fElementVec[mesh] = cel;
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

    /** @brief add an element to the datastructure */
    // NEED TO CHECK IT
    virtual void AddElement(TPZCompEl *cel, int64_t meshindex) override
    {
		if (fElementVec.size() <= meshindex) 
		{
			fElementVec.resize(meshindex+1);
            fActiveApproxSpace.Resize(meshindex+1, 1);
		}
        if (cel)
        {
            TPZGeoEl *gel = cel->Reference();
            TPZCompElSide celside(cel,gel->NSides()-1);
            fElementVec[meshindex] = celside;
        }
        else
        {
            fElementVec[meshindex] = TPZCompElSide();
        }
    }

    virtual TPZCompEl *Element(int64_t elindex)
    {
        DebugStop();
        return 0;
    }

    virtual TPZManVector<TPZCompElSide,5> & ElementVec()
    {
        return fElementVec;
    }

    virtual TPZCompEl *ReferredElement(int64_t mesh)
    {
        DebugStop();
        return 0;
    }

    virtual int64_t NMeshes() override
    {
        return fElementVec.size();
    }

    virtual void SetConnectIndexes(TPZVec<int64_t> &indexes)
    {
        fConnectIndexes = indexes;
    }

    virtual void AffineTransform(TPZVec<TPZTransform<> > &trVec) const
    {
        DebugStop();
    }

    virtual void InitMaterialData(TPZVec<TPZMaterialData > &dataVec, TPZVec<int64_t> *indices = 0)
    {
        DebugStop();
    }

    virtual void PolynomialOrder(TPZVec<int> &order) const
    {
        DebugStop();
    }
};

TPZCompEl * CreateSBFemMultiphysicsCompEl(TPZMultiphysicsCompMesh &mesh, TPZGeoEl *gel, int64_t &index);