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

    set<int> fCondensedMatids = {fDifpressure, fExternalleftflux, fInternal, fExternalrightflux};

    map<int64_t,TPZCompEl *> fGeltocel;

    TPZManVector<int64_t> fElGroupIndexes;
    
public:
    
    /// simple constructor
    TPZBuildSBFemHdiv(TPZAutoPointer<TPZGeoMesh> &gmesh, int skeletonmatid, std::map<int,int> &matidtranslation) : TPZBuildSBFem(gmesh,skeletonmatid, matidtranslation)
    {
    }

    int GetSideSkeletonEl(TPZGeoEl * gel);

    int GetSideCollapsedEl(TPZGeoEl * gel);

    void BuildMultiphysicsCompMesh(TPZMultiphysicsCompMesh &cmesh);

    void CreateExternalElements(TPZGeoMesh * gmesh, set<int> & matidtarget);

    void CreateCollapsedGeoEls(TPZCompMesh & cmeshpressure, set<int> & matidstarget);

    void CreateCompElPressure(TPZCompMesh & cmeshpressure);

    void CreateCompElFlux(TPZCompMesh &cmeshflux, set<int> & matidtarget);

    void CreateSBFEMMultiphysicsMesh(TPZMultiphysicsCompMesh & cmeshm);

    void AddInterfaceElements(TPZMultiphysicsCompMesh & cmeshm);
    
    void CreateSBFEMMultiphysicsVol(TPZMultiphysicsCompMesh & cmeshm, set<int> & matidtarget);

// ONGOING
    // void AdjustExternalPressureConnectivity(TPZMultiphysicsCompMesh & cmeshm);

// NOT READY YET
    void BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh);

    void GroupandCondense();
};