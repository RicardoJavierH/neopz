//
//  TPZBuildSBFemHdiv.h
//  PZ
//
//  Created by Karolinne Coelho on 18/01/21.
//
//

#pragma once

#include "TPZBuildSBFem.h"
#include "TPZMultiphysicsCompMesh.h"

using namespace std;

class TPZBuildSBFemHdiv : public TPZBuildSBFem
{
    int fDifpressure, fInterface, fExternalleftflux, fInternal, fExternalrightflux;

    set<int> fCondensedMatids;

    map<int64_t,TPZCompEl *> fGeltocel;

    TPZManVector<int64_t> fElGroupIndexes;
    
public:
    
    /// simple constructor
    TPZBuildSBFemHdiv(TPZAutoPointer<TPZGeoMesh> &gmesh, int skeletonmatid, std::map<int,int> &matidtranslation) : TPZBuildSBFem(gmesh,skeletonmatid, matidtranslation)
    {
        fDifpressure = fSkeletonMatId + 1;
        fInterface = fDifpressure +1;
        fExternalleftflux = fInterface +1;
        fInternal = fExternalleftflux+1;
        fExternalrightflux = fInternal+1;
        fCondensedMatids = {fExternalleftflux, fInternal, fExternalrightflux};
    }

    int GetSideSkeletonEl(TPZGeoEl * gel);

    int GetSideCollapsedEl(TPZGeoEl * gel);

    void BuildMultiphysicsCompMesh(TPZMultiphysicsCompMesh &cmesh);

    void CreateExternalElements(TPZGeoMesh * gmesh, set<int> & matidtarget);

    void CreateCollapsedGeoEls(TPZCompMesh & cmeshpressure, set<int> & matidstarget);

    void CreateCompElPressure(TPZCompMesh & cmeshpressure, set<int> & matids1d);

    void CreateCompElFlux(TPZCompMesh &cmeshflux, set<int> & matidtarget);

    void CreateSBFEMMultiphysicsMesh(TPZMultiphysicsCompMesh & cmeshm);

    void AddInterfaceElements(TPZMultiphysicsCompMesh & cmeshm, set<int> & matids1d);
    
    void CreateSBFEMMultiphysicsVol(TPZMultiphysicsCompMesh & cmeshm, set<int> &matid1d, set<int> & matidtarget);

    void CreateSBFEMMultiphysicsElGroups(TPZMultiphysicsCompMesh & cmeshm, set<int> & matidtarget);

    void GroupandCondense(TPZMultiphysicsCompMesh & cmeshm);

// NOT READY YET
    void BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh);
};