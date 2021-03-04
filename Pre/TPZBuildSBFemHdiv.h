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
    int fDifPressure, fInterface, fLeftFlux, fRightFlux, fAverPressure, fBC;

    set<int> fCondensedMatids;

    map<int64_t,TPZCompEl *> fGeltocel;

    TPZManVector<int64_t> fElGroupIndexes;
    
public:
    
    /// simple constructor
    TPZBuildSBFemHdiv(TPZGeoMesh * gmesh, int skeletonmatid, std::map<int,int> &matidtranslation) : TPZBuildSBFem(gmesh,skeletonmatid, matidtranslation)
    {
        fDifPressure = fSkeletonMatId + 1;

        fInterface = fDifPressure +1;
        
        fLeftFlux = fInterface +1;
        fRightFlux = fLeftFlux+1;
        
        fAverPressure = fRightFlux+1;
        
        fCondensedMatids = {fLeftFlux, fInterface, fSkeletonMatId, fRightFlux};
    }

    int GetSideSkeletonEl(TPZGeoEl * gel);

    int GetSideCollapsedEl(TPZGeoEl * gel);

    void BuildMultiphysicsCompMesh(TPZMultiphysicsCompMesh &cmesh);

    void CreateExternalElements(TPZGeoMesh * gmesh, set<int> & matidtarget);

    void CreateCollapsedGeoEls(TPZCompMesh & cmeshpressure, set<int> & matidstarget);

    void CreateCompElPressure(TPZCompMesh & cmeshpressure, set<int> & matids1d);

    void CreateCompElFlux(TPZCompMesh &cmeshflux, set<int> & matidtarget, set<int> & matid1d);

    void CreateSBFEMMultiphysicsMesh(TPZMultiphysicsCompMesh & cmeshm);

    void AddInterfaceElements(TPZMultiphysicsCompMesh & cmeshm, set<int> & matids1d);
    
    void CreateSBFEMMultiphysicsVol(TPZMultiphysicsCompMesh & cmeshm, set<int> &matid1d, set<int> & matidtarget);

    void CreateSBFEMMultiphysicsElGroups(TPZMultiphysicsCompMesh & cmeshm, set<int> & matidtarget);

    void GroupandCondense(TPZMultiphysicsCompMesh & cmeshm);

// NOT READY YET
    void BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh);
};