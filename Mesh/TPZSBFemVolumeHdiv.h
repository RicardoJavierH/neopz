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
#include "TPZSBFemMultiphysicsElGroup.h"

#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzelchdiv.h"

using namespace std;
using namespace pzshape;

class TPZSBFemVolumeHdiv : public TPZInterpolationSpace
{
    int fSkeleton;

    int fLeftfluxMatId, fRightfluxMatId;

    TPZCompEl *fElementGroup = 0;

    TPZIntPoints *fIntRule = 0;

public:
    
    TPZSBFemVolumeHdiv(TPZCompMesh &mesh, TPZGeoEl * gel, TPZGeoEl * gel1d, int64_t &index);
    
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

    void SetElementGroupIndex(int64_t index)
    {
        DebugStop();
    }

    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const
    {
        // till I remember how this works
        DebugStop();
        return 0;
    }

    virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                                    std::map<int64_t,int64_t> & gl2lcConMap,
                                    std::map<int64_t,int64_t> & gl2lcElMap) const
    {
        // till I remember how this works
        DebugStop();
        return 0;
    }

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

    virtual void PRefine (int order)
    {
        TPZCompEl *cel = Mesh()->Element(fSkeleton);
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> (cel);
        intel->PRefine(order);
    }
    
    void SetPreferredOrder(int order)
    {
        fPreferredOrder = order;
        TPZCompEl *cel = Mesh()->Element(fSkeleton);
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> (cel);
        intel->SetPreferredOrder(order);
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

    void Shape(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphidxi)
    {
        DebugStop();
    }

    virtual int NConnectShapeF(int icon, int order) const
    {
        TPZConnect &c = Connect(icon);
        return c.NShape();
    }

    virtual int NShapeF() const
    {
        int nc = NConnects();
        int nshape = 0;
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = Connect(ic);
            nshape += c.NShape();
        }
        return nshape;
    }

    virtual int SideConnectLocId(int icon,int is) const
    {
        DebugStop();
		return 0;
    }

    virtual int NSideConnects(int iside) const
    {
        DebugStop();
		return 0;
    }

    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const
    {
        DebugStop();
    }

    virtual void SetConnectIndex(int inode, int64_t index)
    {
        DebugStop();
    }
};

TPZCompEl * CreateSBFemHdivCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, TPZGeoEl * gel1d, int64_t &index);

// #include "pzshapetriang.h"
// #include "pzshapepoint.h"
// #include "pzshapelinear.h"
// #include "pzshapequad.h"

// using namespace pzshape;

// template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeLinear>>;
// template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeTriang>>;
// template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeQuad>>;

// template class TPZCompElHDivSBFem<TPZShapeTriang>;
// template class TPZCompElHDivSBFem<TPZShapeLinear>;
// template class TPZCompElHDivSBFem<TPZShapeQuad>;