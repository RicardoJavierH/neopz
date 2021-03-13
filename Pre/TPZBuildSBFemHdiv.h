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
    int fDifPressure, fInternal, fLeftFlux, fRightFlux, fInterface;

    set<int> fCondensedMatids;

    map<int64_t,TPZCompEl *> fGeltocel;

    TPZManVector<int64_t> fElGroupIndexes;
    
public:
    
    /// simple constructor
    TPZBuildSBFemHdiv(TPZGeoMesh * gmesh, int skeletonmatid, std::map<int,int> &matidtranslation) : TPZBuildSBFem(gmesh,skeletonmatid, matidtranslation)
    {
        // fSkeletonMatId represents the external average pressure;
        fInterface = fSkeletonMatId+1;

        // fRightFlux represents the external right flux (next to fSkeletonMatId)
        fRightFlux = fInterface+1;

        // fInternal will gather internal flux and pressures (next to fRightFlux)
        fInternal = fRightFlux+1;
        
        // fLeftFlux represents the external left flux (next to fInternal)
        fLeftFlux = fInternal +1;

        // fDifPressure represents the differential of the pressure (next to fLeftFlux)
        fDifPressure = fLeftFlux + 1;
        
        // Order of the elements, from left to right:
        // fDifPressure -> fInterface -> fLeftFlux -> fInternal -> fRightFlux -> fInterface -> fSkeletonMatId
        fCondensedMatids = {fInterface, fLeftFlux, fInternal, fRightFlux};
    }

    int GetSideSkeletonEl(TPZGeoEl * gel);

    int GetSideCollapsedEl(TPZGeoEl * gel);

    void BuildMultiphysicsCompMesh(TPZMultiphysicsCompMesh &cmesh);

    void CreateExternalElements(TPZGeoMesh * gmesh, set<int> & matidtarget);

    void CreateCollapsedGeoEls(TPZCompMesh & cmeshpressure, set<int> & matidstarget, set<int> & matids1d);

    void CreateCompElPressure(TPZCompMesh &cmeshpressure, set<int> & matids1d);

    void CreateCompElFlux(TPZCompMesh &cmeshflux, set<int> & matidtarget, set<int> & matid1d);

    void CreateSBFEMMultiphysicsMesh(TPZMultiphysicsCompMesh & cmeshm);

    void AddInterfaceElements(TPZMultiphysicsCompMesh & cmeshm, set<int> & matids1d);
    
    void CreateSBFEMMultiphysicsVol(TPZMultiphysicsCompMesh & cmeshm, set<int> &matid1d, set<int> & matidtarget);

    void CreateSBFEMMultiphysicsElGroups(TPZMultiphysicsCompMesh & cmeshm, set<int> & matidtarget);

    void GroupandCondense(TPZMultiphysicsCompMesh & cmeshm);

    void StandardConfigurationHdiv();

    void AddInternalElements();

// NOT READY YET
    void BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh);
};