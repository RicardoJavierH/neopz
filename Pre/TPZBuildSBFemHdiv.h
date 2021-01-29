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
    //Order of MatIds: fLeftpressure, fRightpressure, fLeftflux, fRightflux;
    // TPZManVector<int, 4> fMatIdsHdiv;
    int fLeftfluxMatId, fRightfluxMatId;

    int fLeftpressureMatId, fRightpressureMatId;

    int fInterfaceMatId;

    map<int64_t,TPZCompEl *> fGeltocel;
    
public:
    
    /// simple constructor
    TPZBuildSBFemHdiv(TPZAutoPointer<TPZGeoMesh> &gmesh, int skeletonmatid, std::map<int,int> &matidtranslation) : TPZBuildSBFem(gmesh,skeletonmatid, matidtranslation)
    {
        fLeftpressureMatId = skeletonmatid+1;
        fRightpressureMatId = fLeftpressureMatId+1;
        fLeftfluxMatId = fRightpressureMatId+1;
        fRightfluxMatId = fLeftfluxMatId+1;
        fInterfaceMatId = fRightfluxMatId+1;
    }

    void BuildMultiphysicsCompMesh(TPZCompMesh &cmesh);

    void CreateExternalElements(TPZGeoMesh * gmesh);

    void CreateVolumetricElementsHdiv(TPZCompMesh &cmesh);

    void CreateSBFemVolumeHdiv(TPZCompMesh & cmesh, set<int> & matidstarget);

    void CreateSBFemDiscontinuousElements(TPZCompMesh &cmeshpressure);

    void UpdateMultiphysicsMesh(TPZManVector<TPZCompMesh*, 2> & cmeshvec, TPZMultiphysicsCompMesh & cmeshm);

    void CreateSBFemMultiphysicsElGroups(TPZMultiphysicsCompMesh & cmeshm);

// NOT READY YET
    void BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh);

    void AdjustFluxConnectivities(TPZCompMesh & cmesh);

    void AdjustExternalPressureConnectivity(TPZCompMesh & cmesh);

    void CreateSBFemInterfaceElementGroups(TPZCompMesh & cmeshm);

    void GroupandCondense();
};