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

#include "pzcondensedcompel.h"
#include "pzgeoelbc.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzbuildsbfem"));
#endif

// This function will create both TPZSBFemMultiphysicsVol and TPZSBFemMultiphysicsElGroup.
// Sequence of the method:
// 1. Define geometry with the external dim-1 elements and collapsed element;
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
    // 1. Material ID defined in the map EGroup -> EMatVol
    auto cmeshflux = cmeshvec[0];    

    // The pressure cmesh must have:
    // 1. Material ID defining the map EGroup -> EMatVol.
    // 2. The Material ID related to the Skeleton. It will be a Neumann BC - external average pressure.
    // 3. The boundary conditions.
    // All these materials must be defined in the cmeshm too.
    auto cmeshpressure = cmeshvec[1];

    // ********** CREATING THE GEOMETRY
    // Skeleton Elements + Collapsed Elements + External elements:

    // Before calling this function the Skeleton elements has been already created
    // So I just need to create the collapsed and external elements

    // Creating dim-1 elements CompEls - Skeleton
    int dim = cmeshm.Dimension();
    fGMesh->SetName("gmesh with collapsed els");
    cmeshpressure->SetName("cmesh pressure");
    
    fGMesh->SetReference(cmeshpressure);
    cmeshpressure->SetReference(fGMesh);

    set<int> matids1d;
    for (auto gel : fGMesh->ElementVec())
    {
        if (!gel)
        {
            continue;
        }
        if (gel->Dimension() == fGMesh->Dimension()-1) {
            matids1d.insert(gel->MaterialId());
            if(gel->MaterialId() != fSkeletonMatId)
            {
                fBC = gel->MaterialId();
                fCondensedMatids.insert(gel->MaterialId());
            }
        }
    }
    cmeshpressure->ApproxSpace().CreateDisconnectedElements(true);
    cmeshpressure->AutoBuild(matids1d);
#ifdef PZDEBUG
    if(1)
    {
        ofstream sout("cmeshpressure0.txt");
        cmeshpressure->Print(sout);
        ofstream gout("gmeshwithoutcollapsedels.txt");
        fGMesh->Print(gout);
    }
#endif

    // Creating volumetric collapsed elements
    set<int> matidstarget; // matid of the collapsed element (output parameter)
    // The updated mesh will be fGMesh
    // But the comp mesh used must be cmeshpressure because the connectivity of the multiphysics mesh
    // hasn't been created yet.
    CreateCollapsedGeoEls(*cmeshpressure, matidstarget, matids1d);
#ifdef PZDEBUG
    if(1)
    {
        ofstream gout("gmeshwithcollapsedels.txt");
        fGMesh->Print(gout);
    }
#endif
    
    cmeshm.SetReference(fGMesh);

    // Creating dim-1 elements: External flux and external pressure.
    CreateExternalElements(fGMesh, matidstarget);

    // ********** CREATING THE COMPUTATIONAL MESH
    // At this point the code created all geometric elements needed for ther SBFEM simulation and it's stored in fGMesh
    // Here it's created the SBFemVolumeHdiv elements -> CompElHdivElement -> Adjust connectivities
    // But firstly the atomic meshes must be properly created
    CreateCompElPressure(*cmeshpressure, matids1d);
    CreateCompElFlux(*cmeshflux, matidstarget, matids1d);

    // Then, creating the multiphysics mesh
    // Next method will setup the multiphysics mesh
    // Deleting the Hdiv elements existing in the mesh.
    // The elements must be specified as TPZCompElHDivSBFem elements
    CreateSBFEMMultiphysicsMesh(cmeshm);
    cmeshm.SetName("multiphysicssbfem");
    TPZManVector<int> active(2,1);
    cmeshm.BuildMultiphysicsSpace(active, cmeshvec);
    cmeshm.LoadReferences();
    cmeshm.CleanUpUnconnectedNodes();
    {
        ofstream sout("cmeshmultiphysics.txt");
        cmeshm.Print(sout);
    }


    // Adding the interface elements
    AddInterfaceElements(cmeshm, matidstarget);
    {
        ofstream sout("cmeshmultiphysics.txt");
        cmeshm.Print(sout);
    }

    // Creating the SBFEMVolumeHdiv elements
    CreateSBFEMMultiphysicsVol(cmeshm, matids1d, matidstarget);

    // Based on the multiphysics mesh, create SBFemElementGroups
    // The SBFemMultiphysicsElGroups = Group of collapsed flux elements (SBFemVolumes) + pressure compels
    // Then, condense the flux into the pressures.
    CreateSBFEMMultiphysicsElGroups(cmeshm, matidstarget);
    cmeshm.ExpandSolution();
    cmeshm.ComputeNodElCon();
    cmeshm.CleanUpUnconnectedNodes();

#ifdef PZDEBUG
    {
        ofstream mout("cmeshmultiphysics.txt");
        cmeshm.Print(mout);
    }
#endif
    GroupandCondense(cmeshm);
    cmeshm.ExpandSolution();
    cmeshm.ComputeNodElCon();
    cmeshm.CleanUpUnconnectedNodes();

#ifdef PZDEBUG
    {
        ofstream mout("cmeshmultiphysicscondensed.txt");
        cmeshm.Print(mout);
    }
#endif

    // auto ElementVec = cmeshm.ElementVec();
    // auto nel = cmeshm.NElements();
    // for (int64_t el = 0; el<nel; el++) {
    //     TPZCompEl *cel = cmeshm.Element(el);
    //     auto sbfemgr = dynamic_cast<TPZSBFemMultiphysicsElGroup * >(cel);
    //     if (sbfemgr)
    //     {
    //         continue;
    //     }
    //     auto sbfemvol = dynamic_cast<TPZSBFemVolumeHdiv * >(cel);
    //     if (sbfemvol)
    //     {
    //         continue;
    //     }
    //     auto sbfemcondensed = dynamic_cast<TPZCondensedCompEl * >(cel);
    //     if (sbfemcondensed)
    //     {
    //         continue;
    //     }
        
    //     if(cel)
    //     {
    //         delete cel;
    //         ElementVec[el] = 0;
    //     }
    // }
#ifdef PZDEBUG
    ofstream sout("cmeshmultiphysics.txt");
    cmeshm.Print(sout);
    ofstream fout("cmeshflux.txt");
    cmeshflux->Print(fout);
    ofstream pout("cmeshpressure.txt");
    cmeshpressure->Print(pout);
#endif
}

// Creates geometric collapsed elements
void TPZBuildSBFemHdiv::CreateCollapsedGeoEls(TPZCompMesh & cmeshpressure, set<int> & matidstarget, set<int> & matids1d)
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

    // Creating geometric collapsed elements for the multiphysics mesh based on the pressure mesh
    auto gmeshpressure = cmeshpressure.Reference();
    auto dim = gmeshpressure->Dimension();

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
            // TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gelpressure,is);

            auto subgelside = gelside.Neighbour();
            auto it = matids1d.find(subgelside.Element()->MaterialId());
            while (subgelside != gelside)
            {
                if (it != matids1d.end()) break;
                subgelside = subgelside.Neighbour();
                it = matids1d.find(subgelside.Element()->MaterialId());
            }
            
            // if(it != matids1d.end())
            {

            // bool onlyinterpolated = false;
            // bool removeduplicates = true;
            // we identify all computational elements connected to this geometric element side
            // gelside.EqualorHigherCompElementList2(celstack, onlyinterpolated, removeduplicates);

            // we create a volume element based on all smaller elements linked to this side
            // for (auto icelstack : celstack)
            // {
                // Creating Duffy elements:
                int64_t index;
                // TPZGeoElSide subgelside = icelstack.Reference();
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

// The order is: fDifPressure - fFluxRight - fFluxLeft - fAverPressure - fSkeleton
// fSkeleton is the internal flux and pressures.
void TPZBuildSBFemHdiv::CreateExternalElements(TPZGeoMesh * gmesh, set<int> & matidtarget)
{
    auto nelpartitions = fElementPartition.size();
    fElementPartition.Resize(nelpartitions*6);
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

        auto idcollapsed = fElementPartition[gel->Index()];

        // getting the side in which the 1d elements will be constructed as neighbours
        auto iside = GetSideCollapsedEl(gel);
        
        TPZGeoElBC(gel, iside, fAverPressure);
        fElementPartition[gmesh->NElements()-1] = idcollapsed;

        TPZGeoElBC(gel, iside, fRightFlux);
        fElementPartition[gmesh->NElements()-1] = idcollapsed;

        TPZGeoElBC(gel, iside, fLeftFlux);
        fElementPartition[gmesh->NElements()-1] = idcollapsed;

        TPZGeoElBC(gel, iside, fDifPressure);
        fElementPartition[gmesh->NElements()-1] = idcollapsed;
    }
}

// 
void TPZBuildSBFemHdiv::CreateCompElPressure(TPZCompMesh &cmeshpressure, set<int> & matids1d)
{
    auto dim = cmeshpressure.Dimension()-1; // materials with dim-1 dimensional elements
    auto nstate = cmeshpressure.MaterialVec().begin()->second->NStateVariables();
    cmeshpressure.SetDimModel(dim);
    cmeshpressure.SetAllCreateFunctionsContinuous();
    cmeshpressure.ApproxSpace().CreateDisconnectedElements(true);

    // Pressure mesh will be composed of:
    // Internal DOFs - Materials already created!

    // Dif pressure
    auto matleft = new TPZNullMaterial(fDifPressure, dim, nstate);
    cmeshpressure.InsertMaterialObject(matleft);

    // Average pressure
    auto matright = new TPZNullMaterial(fAverPressure, dim, nstate);
    cmeshpressure.InsertMaterialObject(matright);

    // ORDERING ELEMENTS

    // 1st. Internal connects
    // BCs and fSkeleton elements have already built

    // 2nd. Dif. pressure
    set<int> matids = {fDifPressure};
    cmeshpressure.AutoBuild(matids);
    
    // 3rd. Average pressure
    matids = {fAverPressure};
    cmeshpressure.AutoBuild(matids);

    for(auto newnod : cmeshpressure.ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }

#ifdef PZDEBUG
    ofstream pout("cmeshpressure.txt");
    cmeshpressure.Print(pout);
#endif

}

void TPZBuildSBFemHdiv::CreateCompElFlux(TPZCompMesh &cmeshflux, set<int> & matidtarget, set<int> & matid1d)
{
    // Now creating the material for the external elements:
    auto dim = cmeshflux.Dimension()-1; // materials with dim-1 dimensional elements
    auto nstate = cmeshflux.MaterialVec().begin()->second->NStateVariables();
    
    auto matboundleft = new TPZNullMaterial(fLeftFlux, dim, nstate);
    cmeshflux.InsertMaterialObject(matboundleft);
    
    auto matboundright = new TPZNullMaterial(fRightFlux, dim, nstate);
    cmeshflux.InsertMaterialObject(matboundright);

    map<int64_t,TPZCompEl *> geltocel;

    fGMesh->ResetReference();
    for (auto gel1d : fGMesh->ElementVec())
    {
        if (!gel1d) continue;
        auto it = matid1d.find(gel1d->MaterialId());
        if (it == matid1d.end()) continue;

        auto side = GetSideSkeletonEl(gel1d);
        TPZGeoElSide gelside1d(gel1d, side);

        auto gelsidecollapsed = gelside1d.Neighbour();
        auto matidcollapsed = gelsidecollapsed.Element()->MaterialId();
        it = matidtarget.find(matidcollapsed);
        while (it == matidtarget.end())
        {
            gelsidecollapsed = gelsidecollapsed.Neighbour();
            matidcollapsed = gelsidecollapsed.Element()->MaterialId();
            it = matidtarget.find(matidcollapsed);
            if(gelsidecollapsed.Element() == gelside1d.Element())
            {
                continue;
            }
        }

        int64_t index;
        auto celhdivc = new TPZCompElHDivSBFem<pzshape::TPZShapeLinear>(cmeshflux, gel1d, gelsidecollapsed, index);
        geltocel[gel1d->Index()] = celhdivc;
    }
        cmeshflux.ExpandSolution();
        cmeshflux.Reference()->ResetReference();

    // Now setting the geoelement for the HdivBound
    for (auto gelleft : fGMesh->ElementVec())
    {
        if (!gelleft || gelleft->MaterialId() != fLeftFlux) continue;

        TPZGeoElSide gelsideleft(gelleft, gelleft->NSides()-1);

        auto gelsideright = gelsideleft.Neighbour();
#ifdef PZDEBUG
        if(gelsideright.Element()->MaterialId() != fRightFlux)
        {
            DebugStop();
        }
#endif

        auto gelsideskeleton = gelsideright.Neighbour();
        auto it = matid1d.find(gelsideskeleton.Element()->MaterialId());
        while(it == matid1d.end())
        {
            gelsideskeleton = gelsideskeleton.Neighbour();
            it = matid1d.find(gelsideskeleton.Element()->MaterialId());
        }

        auto cel  = geltocel[gelsideskeleton.Element()->Index()];
        TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * celhdivc = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * >(cel);
#ifdef PZDEBUG
        if (!celhdivc)
        {
            DebugStop();
        }
#endif

        int64_t index;

        auto hdivboundleft = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmeshflux,gelsideleft.Element(),index);
        hdivboundleft->SetConnectIndex(0,celhdivc->ConnectIndex(3));
        gelsideleft.Element()->ResetReference();

        auto hdivboundright = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmeshflux,gelsideright.Element(),index);
        hdivboundright->SetConnectIndex(0,celhdivc->ConnectIndex(4));
        gelsideright.Element()->ResetReference();
    }

    cmeshflux.LoadReferences();
    cmeshflux.Reference()->ResetReference();
    cmeshflux.CleanUpUnconnectedNodes();
    cmeshflux.ExpandSolution();

#ifdef PZDEBUG
    ofstream sout("cmeshflux0.txt");
    cmeshflux.Print(sout);
#endif
}

void TPZBuildSBFemHdiv::CreateSBFEMMultiphysicsMesh(TPZMultiphysicsCompMesh & cmeshm)
{
    // The materials for ESkeleton and Emat1 were already created before.
    auto dim = cmeshm.Dimension();
    auto nstate = cmeshm.MaterialVec().begin()->second->NStateVariables();
    
    {
        auto mat = new TPZNullMaterial(fDifPressure, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }

    {
        auto mat = new TPZNullMaterial(fAverPressure, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }

    {
        auto mat = new TPZNullMaterial(fLeftFlux, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }

    {
        auto mat = new TPZNullMaterial(fRightFlux, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }

    {
        auto mat = new TPZLagrangeMultiplier(fInterface, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }
    cmeshm.SetAllCreateFunctionsMultiphysicElem();
}

void TPZBuildSBFemHdiv::AddInterfaceElements(TPZMultiphysicsCompMesh & cmeshm, set<int> & matidtarget)
{
    for (auto geldifpr : fGMesh->ElementVec())
    {
        if (!geldifpr || geldifpr->MaterialId() != fDifPressure) continue;

        auto side = GetSideSkeletonEl(geldifpr);
        TPZGeoElSide gelsidedifpr(geldifpr, side);

        auto gelsideleft = gelsidedifpr.Neighbour(); // fLeftFlux
#ifdef PZDEBUG
        if(gelsideleft.Element()->MaterialId() != fLeftFlux)
        {
            DebugStop();
        }
#endif
        int64_t index;

        TPZCompElSide celsidedifpr = gelsidedifpr.Reference();
        TPZCompElSide celsideleft = gelsideleft.Reference();
        if (!celsidedifpr || !celsideleft)
        {
            DebugStop();
        }
        if (fElementPartition[gelsidedifpr.Element()->Index()] != fElementPartition[gelsideleft.Element()->Index()])
        {
            DebugStop();
        }
        TPZGeoElBC gelbc(gelsidedifpr, fInterface);
        TPZMultiphysicsInterfaceElement *intl = new TPZMultiphysicsInterfaceElement(cmeshm, gelbc.CreatedElement(),index,celsidedifpr,celsideleft);
        fElementPartition[gelbc.CreatedElement()->Index()] = fElementPartition[gelsidedifpr.Element()->Index()];
    }

    for (auto gelright : fGMesh->ElementVec())
    {
        if (!gelright || gelright->MaterialId() != fRightFlux) continue;

        // 1st neighbour is fDifPressure, 2nd is fLeftflux, 3rd is fRightflux
        // 4th is fAverPressure, 5th fSkeleton (InternalFlux)
        auto side = GetSideSkeletonEl(gelright);
        TPZGeoElSide gelsideright(gelright, side);

        auto gelsideaverpr = gelsideright.Neighbour(); // fLeftFlux
#ifdef PZDEBUG
        if(gelsideaverpr.Element()->MaterialId() != fAverPressure)
        {
            DebugStop();
        }
#endif
        int64_t index;

        TPZCompElSide celsideaverpr = gelsideaverpr.Reference();
        TPZCompElSide celsideright = gelsideright.Reference();
        if (!celsideaverpr || !celsideright)
        {
            DebugStop();
        }
        if (fElementPartition[gelsideaverpr.Element()->Index()] != fElementPartition[gelsideright.Element()->Index()])
        {
            DebugStop();
        }
        TPZGeoElBC gelbc(gelsideaverpr, fInterface);
        TPZMultiphysicsInterfaceElement *intr = new TPZMultiphysicsInterfaceElement(cmeshm, gelbc.CreatedElement(),index,celsideright,celsideaverpr);
        fElementPartition[gelbc.CreatedElement()->Index()] = fElementPartition[gelsideaverpr.Element()->Index()];
    }
}

// Here SBFEM Multiphysics Volume elements are created.
// The Multiphysics element must be fSkeleton
// But I still need to pass the Volumetric Collapsed element so the TPZSBFemVolumeHdiv class will be able to perform integrations.
void TPZBuildSBFemHdiv::CreateSBFEMMultiphysicsVol(TPZMultiphysicsCompMesh & cmeshm, set<int> &matids1d, set<int> & matidtarget)
{

#ifdef PZDEBUG
    ofstream outvtk("GeometryHybridSBFEM.vtk");
    TPZVTKGeoMesh vtk;
    vtk.PrintGMeshVTK(fGMesh, outvtk, true);
    ofstream gout("GeometryHybridSBFEM.txt");
    fGMesh->Print(gout);
#endif

    // This for creates the volumetric multiphysics element, working as an AutoBuild()
    for (auto gelcollapsed : fGMesh->ElementVec())
    {
        // Now I'm searching for the collapsed element, in which the matid is in matidtarget.
        if (!gelcollapsed)
        {
            continue;
        }
        auto it = matidtarget.find(gelcollapsed->MaterialId());
        if (it == matidtarget.end())
        {
            continue;
        }

        auto idvol = fElementPartition[gelcollapsed->Index()];
        if (idvol == -1)
        {
            DebugStop();
        }

        int64_t index;
        auto cel = CreateSBFemMultiphysicsCompEl(cmeshm, gelcollapsed, index);

        // Getting the TPZSBFemVolumeHdiv compel
        auto celsbfem = cmeshm.Element(index);
        auto sbfem = dynamic_cast<TPZSBFemVolumeHdiv * >(celsbfem);
        if(!sbfem)
        {
            DebugStop();
        }

        // Adding multiphysics compels:
        // The neighbour will be in the surface identified as side. So I get the gelside and search for its neighbour.
        // I'll include all computational elements related to geo collapsed element
        auto side = GetSideCollapsedEl(gelcollapsed);
        TPZGeoElSide gelside(gelcollapsed, side);

        // 1st element will be fDifPressure
        auto neigh = gelside.Neighbour();
        auto gel1d = neigh.Element();
        if(!gel1d || !(gel1d->Reference()))
        {
            DebugStop();
        }
        while (gel1d->MaterialId() != fDifPressure || fElementPartition[gel1d->Index()] != fElementPartition[idvol])
        {
            neigh = neigh.Neighbour();
            gel1d = neigh.Element();
        }
        auto cel1d = gel1d->Reference();
        sbfem->AddElement(cel1d, 0);
        int count = 1;

        while (gelside != neigh)
        {
            neigh = neigh.Neighbour();
            gel1d = neigh.Element();
            if(gel1d->Reference() && gel1d->Dimension() == fGMesh->Dimension()-1)
            {
                if(fElementPartition[gel1d->Index()] == idvol && count != 6)
                {
                    auto cel1d = gel1d->Reference();
                    sbfem->AddElement(cel1d, count);
                    count++;
                    if(count == 7) break;
                }
                auto it = matids1d.find(gel1d->MaterialId());
                if(fElementPartition[gel1d->Index()] == -1 && it != matids1d.end() && count == 6)
                {
                    auto cel1d = gel1d->Reference();
                    sbfem->AddElement(cel1d, count);
                    count++;
                    if(count == 7) break;
                }
                if(gel1d->MaterialId() == fSkeletonMatId && count == 6)
                {
                    auto cel1d = gel1d->Reference();
                    sbfem->AddElement(cel1d, count);
                    count++;
                    if(count == 7) break;
                }
            }
        }
        sbfem->Print(cout);
    }
}

void TPZBuildSBFemHdiv::GroupandCondense(TPZMultiphysicsCompMesh & cmeshm)
{
    for (auto cel : cmeshm.ElementVec())
    {
        if (!cel)
        {
            continue;
        }
        auto sbfemgr = dynamic_cast<TPZSBFemMultiphysicsElGroup *>(cel);
        if (!sbfemgr)
        {
            continue;
        }
        sbfemgr->GroupandCondense(fCondensedMatids);
        cmeshm.ComputeNodElCon();
        cmeshm.CleanUpUnconnectedNodes();
    }
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
        side = 2;
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

void TPZBuildSBFemHdiv::CreateSBFEMMultiphysicsElGroups(TPZMultiphysicsCompMesh & cmesh, set<int> & matidtarget)
{
    int64_t numgroups = fPartitionCenterNode.size();
    int64_t groupelementindices(numgroups);
    
    TPZManVector<int64_t> elementgroupindices(numgroups); 
    
    for (int64_t el=0; el<numgroups; el++)
    {
        int64_t index;
        new TPZSBFemMultiphysicsElGroup(cmesh,index);
        elementgroupindices[el] = index;
    }
    int dim = cmesh.Dimension();
    for (auto cel : cmesh.ElementVec())
    {
        if (!cel) continue;
        
        auto sbfem = dynamic_cast<TPZSBFemVolumeHdiv *>(cel);
        if (sbfem)
        {
            TPZGeoEl *gel = sbfem->Reference();
            int64_t gelindex = gel->Index();

            int side = GetSideCollapsedEl(gel);
            TPZGeoElSide gelside(gel,side);
            int geldim = gel->Dimension();
            int nsidenodes = gel->NSideNodes(side);
            TPZManVector<int64_t,8> cornernodes(nsidenodes);
            for (int node = 0; node<nsidenodes; node++)
            {
                cornernodes[node] = gel->SideNodeIndex(side, node);
            }
            
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if(neighbour.Element()->Dimension() == geldim-1 && neighbour.Element()->Reference())
                {
                    int nsidenodesneighbour = neighbour.Element()->NCornerNodes();
                    if (nsidenodesneighbour == nsidenodes)
                    {
                        TPZManVector<int64_t,8> neighnodes(nsidenodesneighbour);
                        for (int i=0; i<nsidenodesneighbour; i++) {
                            neighnodes[i] = neighbour.Element()->NodeIndex(i);
                        }
                        int numequal = 0;
                        for (int i=0; i<nsidenodesneighbour; i++) {
                            if (neighnodes[i] == cornernodes[i]) {
                                numequal++;
                            }
                        }
                        if (numequal == nsidenodesneighbour) {
                            break;
                        }
                    }
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                // we are not handling open sides (yet)
                DebugStop();
            }
            int64_t skelindex = neighbour.Element()->Reference()->Index();
            sbfem->SetSkeleton(skelindex);
            
            int64_t gelgroup = fElementPartition[gelindex];
            if (gelgroup == -1) {
                DebugStop();
            }
            int64_t celgroupindex = elementgroupindices[gelgroup];
            auto celgr = cmesh.Element(celgroupindex);
            auto sbfemgr = dynamic_cast<TPZSBFemMultiphysicsElGroup *>(celgr);
            if (!sbfemgr) {
                DebugStop();
            }
            sbfemgr->AddElement(sbfem);
        }
    }
    
    for (auto index : elementgroupindices)
    {
        auto cel = cmesh.Element(index);
        auto sbfemgroup = dynamic_cast<TPZSBFemMultiphysicsElGroup *>(cel);
        if (!sbfemgroup)
        {
            DebugStop();
        }
        const TPZVec<TPZCompEl *> &subgr = sbfemgroup->GetElGroup();
        int64_t nsub = subgr.NElements();
        for (auto cels : subgr)
        {
            TPZSBFemVolumeHdiv *femvol = dynamic_cast<TPZSBFemVolumeHdiv *>(cels);
            if (!femvol) {
                DebugStop();
            }
            femvol->SetElementGroupIndex(index);
        }
    }
}