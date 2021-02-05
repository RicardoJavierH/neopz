//
//  TPZBuildSBFemHdiv.cpp
//  PZ
//
//  Created by Karolinne Coelho on 18/01/21.
//
//

#include "TPZBuildSBFemHdiv.h"

#include "TPZCompElHDivSBFem.h"
#include "TPZNullMaterial.h"
#include "TPZLagrangeMultiplier.h"
#include "pzgeoelbc.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzbndcond.h"
#include "TPZVTKGeoMesh.h"

#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"

#include "TPZSBFemMultiphysicsElGroup.h"
#include "TPZSBFemVolumeHdiv.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzbuildsbfem"));
#endif

// This function will create both TPZSBFemMultiphysicsVol and TPZSBFemMultiphysicsElGroup.
// Sequence of the method:
// 1. Define geometry with the hybrid dim-1 elements;
// 2. Update atomic meshes (cmeshflux and cmeshpressure);
// 3. Update Multiphysics mesh;
// 4. Create interface elements;
// 5. Create SBFemMultiphysicsVol objects based on the multiphysics elements;
// 6. Update the multiphysics mesh to contemplate only the SBFem elements (ignore FE els);
// 7. Define SBFemMultiphysicsElGroups as a group of SBFemMultiphysicsVol;
// 8. Condense the DOFs of the SBFemMultiphysicsElGroups.
void TPZBuildSBFemHdiv::BuildMultiphysicsCompMesh(TPZMultiphysicsCompMesh & cmeshm)
{
    TPZManVector<TPZCompMesh*, 2> cmeshvec = cmeshm.MeshVector();

    // The flux cmesh must have:
    // 1. Material ID defining the map EGroup -> EMatVol only!
    // The Skeleton will be defined in the pressure and multiphysics mesh only, since it represents the external pressure.
    auto cmeshflux = cmeshvec[0];

    // The pressure cmesh must have:
    // 1. Material ID defining the map EGroup -> EMatVol.
    // 2. The Material ID related to the Skeleton. It will be a Neumann BC.
    // 3. The boundary conditions.
    // All these materials must be defined in the cmeshm too.
    auto cmeshpressure = cmeshvec[1];

    // ********** CREATING THE GEOMETRY
    // Skeleton Elements + Collapsed Elements + External elements:

    // Before calling this function the Skeleton elements has been already created
    // So I just need to create the collapsed and external elements

    // Creating dim-1 elements CompEls - Skeleton and BCs
    int dim = cmeshm.Dimension();
    auto fGMesh = cmeshm.Reference();

    set<int> matids;
    for (auto gel : fGMesh->ElementVec())
    {
        if (!gel)
        {
            continue;
        }
        if (gel->Dimension() < dim) {
            matids.insert(gel->MaterialId());
        }
    }
    cmeshpressure->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshpressure->AutoBuild(matids);

    // Creating volumetric collapsed elements
    set<int> matidstarget; // matid of the collapsed element (output parameter)
    // The updated mesh will be fGMesh
    // But the comp mesh used must be cmeshpressure because the connectivity of the multiphysics mesh
    // hasn't been created yet.
    CreateCollapsedGeoEls(*cmeshpressure, matidstarget);
    cmeshm.SetReference(fGMesh);

    // Creating dim-1 elements: External flux and external pressure.
    CreateExternalElements(fGMesh, matidstarget);

#ifdef PZDEBUG
    ofstream outvtk("GeometryHybridSBFEM.vtk");
    TPZVTKGeoMesh vtk;
    vtk.PrintGMeshVTK(fGMesh, outvtk, true);
    ofstream gout("GeometryHybridSBFEM.txt");
    fGMesh->Print(gout);
#endif

    // ********** CREATING THE COMPUTATIONAL MESH
    // At this point the code created all geometric elements needed for ther SBFEM simulation and it's stored in fGMesh
    // Here it's created the SBFemVolumeHdiv elements -> CompElHdivElement -> Adjust connectivities
    // But firstly the atomic meshes must be properly created
    CreateCompElPressure(*cmeshpressure);
    CreateCompElFlux(*cmeshflux, matidstarget);

    // Then, creating the multiphysics mesh
    CreateSBFEMMultiphysicsMesh(cmeshm);
    TPZManVector<int> active(2,1);
    cmeshm.BuildMultiphysicsSpace(active, cmeshvec);
    cmeshm.LoadReferences();
    cmeshm.CleanUpUnconnectedNodes();

    AddInterfaceElements(cmeshm);

    CreateSBFEMMultiphysicsVol(cmeshm, matidstarget);

//     // Based on the multiphysics mesh, create SBFemElementGroups
//     // The SBFemMultiphysicsElGroups = Group of collapsed flux elements (SBFemVolumes) + pressure compels
//     // Then, condense the flux into the pressures.
//     // CreateSBFemFluxElGroups(*cmeshm);
//     // CreateSBFemPressureElGroups(*cmeshm);

//     // Ajusting the connectivity of these elements
//     // AdjustExternalPressureConnectivity(*cmeshm);

// #ifdef PZDEBUG
//     ofstream mout("cmeshmultiphysics.txt");
//     cmeshm->Print(mout);
// #endif
}

// Creates geometric collapsed elements
void TPZBuildSBFemHdiv::CreateCollapsedGeoEls(TPZCompMesh & cmeshpressure, set<int> & matidstarget)
{
    // all computational elements for internal DOFs have been loaded for flux and pressure meshes
    set<int> matids;
    for (auto it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++)
    {
        int64_t mat = it->second;
        if (cmeshpressure.FindMaterial(mat))
        {
            matids.insert(it->first);
            matidstarget.insert(it->second);
        }
    }

    // fGMesh = gmesh for the multiphysics mesh
    // Creating geometric collapsed elements for the multiphysics mesh based on the pressure mesh
    auto dim = fGMesh->Dimension();
    auto gmeshpressure = cmeshpressure.Reference();

    for (auto gelpressure : gmeshpressure->ElementVec())
    {
        if (!gelpressure || gelpressure->HasSubElement() || gelpressure->Reference())
        {
            continue;
        }
        // we create SBFemVolume elements by partitioning the volume elements
        auto el = gelpressure->Index();
        if (gelpressure->Dimension() != dim || fElementPartition[el] == -1)
        {
            continue;
        }
        // I am searching elements with dim dimension and with FE matids
        // to transform it in collapsed elements with matidtargets
        if (matids.find(gelpressure->MaterialId()) == matids.end())
        {
            continue;
        }
        int nsides = gelpressure->NSides();
        for (int is = 0; is<nsides; is++)
        {
            if (gelpressure->SideDimension(is) != dim-1)
            {
                continue;
            }
            TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gelpressure,is);
            bool onlyinterpolated = true;
            bool removeduplicates = true;
            // we identify all computational elements connected to this geometric element side
            gelside.EqualorHigherCompElementList2(celstack, onlyinterpolated, removeduplicates);

            // we create a volume element based on all smaller elements linked to this side
            for (auto icelstack : celstack)
            {
                // Creating Duffy elements:
                int64_t index;
                TPZGeoElSide subgelside = icelstack.Reference();
                if (subgelside.Dimension() != dim-1)
                {
                    continue;
                }
                auto nnodes = subgelside.NSideNodes();
                TPZManVector<int64_t,8> Nodes(nnodes*2,-1);
                for (int in=0; in<nnodes; in++)
                {
                    Nodes[in] = subgelside.SideNodeIndex(in); // Boundary nodes
                }
                int elpartition = fElementPartition[el];
                for (int in=nnodes; in < 2*nnodes; in++)
                {
                    Nodes[in] = fPartitionCenterNode[elpartition]; // Scaling center nodes
                }
                auto matid = fMatIdTranslation[gelpressure->MaterialId()];

                // Then, the collapsed volumetric elements will be created in the multiphysics gmesh
                if (subgelside.IsLinearMapping())
                {
                    switch(nnodes)
                    {
                        case 2:
                            fGMesh->CreateGeoElement(EQuadrilateral, Nodes, matid, index);
                            break;
                        case 4:
                            fGMesh->CreateGeoElement(ECube, Nodes, matid, index);
                            DebugStop();
                            break;
                        case 3:
                            fGMesh->CreateGeoElement(EPrisma, Nodes, matid, index);
                            DebugStop();
                            break;
                        default:
                            std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                            DebugStop();
                    }
                    
                }
                else
                {
                    int64_t elementid = fGMesh->NElements()+1;
                    switch(nnodes)
                    {
                        case 2:
                            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (Nodes, matid, *fGMesh,index);
                            break;
                        case 4:
                            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > (Nodes, matid, *fGMesh,index);
                            break;
                        case 3:
                            fGMesh->CreateGeoElement(EPrisma, Nodes, matid, index);
                            break;
                        default:
                            std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                            DebugStop();
                    }
                }
                if (index >= fElementPartition.size())
                {
                    fElementPartition.resize(index+1);
                }
                fElementPartition[index] = elpartition;
            }
        }
    }
    fGMesh->BuildConnectivity();

#ifdef PZDEBUG
    {
        ofstream gout("gmeshwithvol.txt");
        fGMesh->Print(gout);
    }
#endif
}

// The order is: fDifPressure - fInterface - fExternalLeftflux - fInternal - fExternalRightflux - fInterface - fSkeleton
// fSkeleton is the external pressure.
void TPZBuildSBFemHdiv::CreateExternalElements(TPZGeoMesh * gmesh, set<int> & matidtarget)
{
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel)
        {
            continue;
        }
        // matidtarget contains the material id of the collapsed elements
        auto it = matidtarget.find(gel->MaterialId());
        if (it == matidtarget.end())
        {
            continue;
        }

        // getting the side in which the 1d elements will be constructed as neighbours
        auto iside = GetSideCollapsedEl(gel);
        
        TPZGeoElBC(gel, iside, fDifpressure); // "External" pressure - physically the derivative of the pressure
        TPZGeoElBC(gel, iside, fInterface); // TPZMultiphysicsInterface
        TPZGeoElBC(gel, iside, fExternalleftflux); // HdivBound
        TPZGeoElBC(gel, iside, fInternal); // Hdiv + Internal pressure
        TPZGeoElBC(gel, iside, fExternalrightflux); // HdivBound
        TPZGeoElBC(gel, iside, fInterface); // TPZMultiphysicsInterface
        // fSkeleton will be the external pressures
    }
}

// 
void TPZBuildSBFemHdiv::CreateCompElPressure(TPZCompMesh &cmeshpressure)
{
    auto dim = cmeshpressure.Dimension()-1; // materials with dim-1 dimensional elements
    auto nstate = cmeshpressure.MaterialVec().begin()->second->NStateVariables();

    // Pressure mesh will be composed of:
    // fDifpressure + fInternal + fSkeleton (external pressure)
    auto matinternal = new TPZNullMaterial(fInternal, dim, nstate);
    cmeshpressure.InsertMaterialObject(matinternal);

    auto matleft = new TPZNullMaterial(fDifpressure, dim, nstate);
    cmeshpressure.InsertMaterialObject(matleft);

    set<int> matids = {fInternal, fDifpressure};
    cmeshpressure.SetAllCreateFunctionsContinuous();
    cmeshpressure.ApproxSpace().CreateDisconnectedElements(true);
    cmeshpressure.AutoBuild(matids); // BCs and fSkeleton elements have already built

    for(auto newnod : cmeshpressure.ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }
}

void TPZBuildSBFemHdiv::CreateCompElFlux(TPZCompMesh &cmeshflux, set<int> & matidtarget)
{
    auto dim = cmeshflux.Dimension()-1; // materials with dim-1 dimensional elements
    auto nstate = cmeshflux.MaterialVec().begin()->second->NStateVariables();

    cmeshflux.SetReference(fGMesh);

    auto mat = new TPZNullMaterial(fInternal, dim, nstate);
    cmeshflux.InsertMaterialObject(mat);
    
    auto matboundleft = new TPZNullMaterial(fExternalleftflux, dim, nstate);
    cmeshflux.InsertMaterialObject(matboundleft);
    
    auto matboundright = new TPZNullMaterial(fExternalrightflux, dim, nstate);
    cmeshflux.InsertMaterialObject(matboundright);

    map<int64_t,TPZCompEl *> geltocel;
    for (auto gel1d : fGMesh->ElementVec())
    {
        if (!gel1d) continue;
        if (gel1d->MaterialId() != fInternal) continue;
        int64_t index;
        auto side = GetSideSkeletonEl(gel1d);
        TPZGeoElSide gelside(gel1d, side);
        auto gelcollapsedside = gelside.HasNeighbour(matidtarget); // TEST IT
        if (!gelcollapsedside)
        {
            DebugStop();
        }
        auto celhdivc = new TPZCompElHDivSBFem<pzshape::TPZShapeLinear>(cmeshflux, gel1d, gelcollapsedside, index);
        geltocel[gel1d->Index()] = celhdivc;
    }
    cmeshflux.ExpandSolution();
    cmeshflux.Reference()->ResetReference();
    
    // Now setting the geoelement for the HdivBound
    for (auto gel : fGMesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->MaterialId() != fExternalleftflux) continue;

        TPZGeoElSide gelside(gel,2);
        auto intfluxside = gelside.Neighbour();
        auto cel  = geltocel[intfluxside.Element()->Index()];
        TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * celhdivc = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * >(cel);
        if (!celhdivc)
        {
            DebugStop();
        }

        int64_t index;
        auto hdivboundleft = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmeshflux,gel,index);
        hdivboundleft->SetConnectIndex(0,celhdivc->ConnectIndex(3));
        gel->ResetReference();

        auto rightfluxside = gelside.Neighbour();
        auto gel1dright = rightfluxside.Element();
        auto hdivboundright = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmeshflux,gel1dright,index);
        hdivboundright->SetConnectIndex(0,celhdivc->ConnectIndex(4));
        gel1dright->ResetReference();
    }
    
    cmeshflux.LoadReferences();
    cmeshflux.Reference()->ResetReference();
    cmeshflux.CleanUpUnconnectedNodes();
    cmeshflux.ExpandSolution();
}

void TPZBuildSBFemHdiv::CreateSBFEMMultiphysicsMesh(TPZMultiphysicsCompMesh & cmeshm)
{
    auto dim = cmeshm.Dimension();
    auto nstate = cmeshm.MaterialVec().begin()->second->NStateVariables();

    for (auto matid : fCondensedMatids)    
    {
        auto mat = new TPZNullMaterial(matid, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }
    {
        auto mat = new TPZLagrangeMultiplier(fInterface, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }
}

void TPZBuildSBFemHdiv::AddInterfaceElements(TPZMultiphysicsCompMesh & cmeshm)
{
    for (auto gel : fGMesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != fDifpressure)
        {
            continue;
        }

        auto nsides = gel->NSides();
        TPZGeoElSide gelsidedpr(gel,nsides-1); // gelside for the differential of the pressure
        auto gelsideint = gelsidedpr.Neighbour(); // the neighbour must be the interface
        if (!gelsideint || gelsideint.Element()->MaterialId() != fInterface) //
        {
            DebugStop();
        }
        auto gelsidefll = gelsideint.Neighbour(); // the neighbour of the interface must be the left flux
        if (!gelsidefll || gelsidefll.Element()->MaterialId() != fExternalleftflux)
        {
            DebugStop();
        }
        TPZCompElSide celsidedpr = gelsidedpr.Reference();
        TPZCompElSide celsidefll = gelsidefll.Reference();
        if (!celsidedpr || !celsidefll) {
            DebugStop();
        }
        int64_t index;
        TPZMultiphysicsInterfaceElement *intfl = new TPZMultiphysicsInterfaceElement(cmeshm, gelsideint.Element(),index,celsidefll,celsidedpr);


        auto gelsideinternal = gelsidefll.Neighbour(); // Next neighbour will be internal
        auto gelsideflr = gelsideinternal.Neighbour(); // Next neighbour will be flux right
        if (!gelsideflr || gelsideflr.Element()->MaterialId() != fExternalrightflux)
        {
            DebugStop();
        }
        gelsideint = gelsideflr.Neighbour(); // the neighbour of the right flux must be the interface again
        if (!gelsideint || gelsideint.Element()->MaterialId() != fInterface) //
        {
            DebugStop();
        }
        auto gelsidepr = gelsideint.Neighbour(); // Next neighbour will be the external pressure
        if (!gelsidepr || gelsidepr.Element()->MaterialId() != fSkeletonMatId)
        {
            DebugStop();
        }
        TPZCompElSide celsidepr = gelsidepr.Reference();
        TPZCompElSide celsideflr = gelsideflr.Reference();
        if (!celsidepr || !celsideflr) {
            DebugStop();
        }
        TPZMultiphysicsInterfaceElement *intfr = new TPZMultiphysicsInterfaceElement(cmeshm, gelsideint.Element(),index,celsideflr,celsidepr);

        
    }
}

void TPZBuildSBFemHdiv::CreateSBFEMMultiphysicsVol(TPZMultiphysicsCompMesh & cmeshm, set<int> & matidtarget)
{
    for (auto gel : fGMesh->ElementVec())
    {
        // I am searching for a multiphysics element in which the geometric element is a collapsed element
        // It only creates the volumetric multiphysics element, working as an AutoBuild()
        if (!gel)
        {
            continue;
        }
        auto it = matidtarget.find(gel->MaterialId());
        if (it == matidtarget.end())
        {
            continue;
        }
        int64_t index;
        auto cel = CreateSBFemMultiphysicsCompEl(cmeshm, gel, index);
    }
    cmeshm.AddElements(); // NEED TO TEST IT
    cmeshm.AddConnects();
}

// void TPZBuildSBFemHdiv::AdjustExternalPressureConnectivity(TPZMultiphysicsCompMesh & cmeshm)
// {
//     for (auto cel : cmeshm.ElementVec())
//     {
//         if (!cel)
//         {
//             DebugStop();
//         }
//         auto celgr = dynamic_cast<TPZSBFemMultiphysicsElGroup * >(cel);
//         if (!celgr)
//         {
//             DebugStop();
//         }
        
//         // auto newid = -1;
//         // auto nsides = 4;
//         // auto dim = cel->Reference()->Dimension();

//         // TPZStack<int64_t> internalprcon;

//         // auto ncon = cel->NConnects() - nsides * dim;
//         // perm.resize(ncon);

//         // for (auto con : cmeshm->ConnectVec())
//         // {
//         //     if (con.SequenceNumber() == -1)
//         //     {
//         //         continue;
//         //     }
//         //     perm[con.SequenceNumber()] = con.SequenceNumber();
//         // }

//         // int64_t nf = cmeshf->NConnects() - 2*nsides;
//         // auto id = nf+3*nsides;

//         // for (int is = 0; is < nsides; is++)
//         // {
//         //     for (int ic = 0; ic < 3; ic++)
//         //     {
//         //         auto pos = nf + 3*nsides + is*6 + ic;
//         //         perm[pos] = id;
//         //         id++;
//         //     }
//         // }
//         // for (int is = 0; is < nsides; is++)
//         // {
//         //     for (int ic = 0; ic < 3; ic++)
//         //     {
//         //         auto pos = nf + 3*nsides + is*6 + ic + 3;
//         //         perm[pos] = id;
//         //         id++;
//         //     }
//         // }
//     }
// }

void TPZBuildSBFemHdiv::GroupandCondense()
{
    DebugStop();
}


void TPZBuildSBFemHdiv::BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh)
{
    DebugStop();
}


// ******** UTILITY FUNCTIONS

int TPZBuildSBFemHdiv::GetSideSkeletonEl(TPZGeoEl * gel)
{
    auto side = -1;
    switch (gel->Type())
    {
    case EOned:
        side = 3;
        break;
    case ETriangle:
        side = 6;
        break;
    case EQuadrilateral:
        side = 8;
    default:
        break;
    }
    return side;
}

int TPZBuildSBFemHdiv::GetSideCollapsedEl(TPZGeoEl * gel)
{
    auto side = -1;
    switch (gel->Type())
    {
    case EQuadrilateral:
        side = 4;
        break;
    case EPrisma:
        side = 15;
        break;
    case ECube:
        side = 20;
        break;
    default:
        break;
    }
    return side;
}