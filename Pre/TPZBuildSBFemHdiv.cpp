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
    auto fGMesh = cmeshm.Reference();

    set<int> matids1d;
    for (auto gel : fGMesh->ElementVec())
    {
        if (!gel)
        {
            continue;
        }
        if (gel->Dimension() == fGMesh->Dimension()-1) {
            matids1d.insert(gel->MaterialId());
        }
    }
    cmeshpressure->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshpressure->AutoBuild(matids1d);

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
    // Next method will setup the multiphysics mesh
    CreateSBFEMMultiphysicsMesh(cmeshm);
    cmeshm.SetName("multiphysicssbfem");
    TPZManVector<int> active(2,1);
    cmeshm.BuildMultiphysicsSpace(active, cmeshvec);
    cmeshm.LoadReferences();
    cmeshm.CleanUpUnconnectedNodes();

    // Adding the interface elements
    AddInterfaceElements(cmeshm, matids1d);
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

    auto ElementVec = cmeshm.ElementVec();
    auto nel = cmeshm.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmeshm.Element(el);
        auto sbfemgr = dynamic_cast<TPZSBFemMultiphysicsElGroup * >(cel);
        if (sbfemgr)
        {
            continue;
        }
        auto sbfemvol = dynamic_cast<TPZSBFemVolumeHdiv * >(cel);
        if (sbfemvol)
        {
            continue;
        }
        auto sbfemcondensed = dynamic_cast<TPZCondensedCompEl * >(cel);
        if (sbfemcondensed)
        {
            continue;
        }
        
        if(cel)
        {
            delete cel;
            ElementVec[el] = 0;
        }
    }
    cmeshm.InitializeBlock();
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
            bool onlyinterpolated = false;
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

    cmeshpressure.SetAllCreateFunctionsContinuous();
    cmeshpressure.ApproxSpace().CreateDisconnectedElements(true);
    set<int> matids = {fDifpressure};
    cmeshpressure.AutoBuild(matids); // BCs and fSkeleton elements have already built
    matids = {fInternal};
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

    auto mat = new TPZNullMaterial(fInternal, dim, nstate);
    cmeshflux.InsertMaterialObject(mat);
    
    auto matboundleft = new TPZNullMaterial(fExternalleftflux, dim, nstate);
    cmeshflux.InsertMaterialObject(matboundleft);
    
    auto matboundright = new TPZNullMaterial(fExternalrightflux, dim, nstate);
    cmeshflux.InsertMaterialObject(matboundright);

    map<int64_t,TPZCompEl *> geltocel;
    for (auto gelcollapsed : fGMesh->ElementVec())
    {
        if (!gelcollapsed) continue;

        int midcollapsed = gelcollapsed->MaterialId();
        auto it = matidtarget.find(midcollapsed);
        if (it == matidtarget.end())
        {
            continue;
        }

        // the 3rd neighbour is fInternal
        auto side = GetSideCollapsedEl(gelcollapsed);
        TPZGeoElSide gelcollapsedside(gelcollapsed, side);

        auto gel1dside = gelcollapsedside.Neighbour(); // fInterface
        gel1dside = gel1dside.Neighbour(); // fExtrightflux
        gel1dside = gel1dside.Neighbour(); // fInternal

        auto gel1d = gel1dside.Element();

        if (gel1d->MaterialId() != fInternal)
        {
            DebugStop();
        }

        int64_t index;
        auto celhdivc = new TPZCompElHDivSBFem<pzshape::TPZShapeLinear>(cmeshflux, gel1d, gelcollapsedside, index);
        geltocel[gel1d->Index()] = celhdivc;
        cmeshflux.ExpandSolution();
        cmeshflux.Reference()->ResetReference();
    }
    
    // Now setting the geoelement for the HdivBound
    // Order of the neighbours:
    // fExtrightflux -> fInternal -> fExtleftflux
    for (auto gel : fGMesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->MaterialId() != fExternalrightflux) continue;

        TPZGeoElSide gelside(gel,2);
        auto intfluxside = gelside.Neighbour();
        int mid = intfluxside.Element()->MaterialId();
        while (mid != fInternal)
        {
            DebugStop();
        }
        
        auto cel  = geltocel[intfluxside.Element()->Index()];
        TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * celhdivc = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * >(cel);
        if (!celhdivc)
        {
            DebugStop();
        }

        int64_t index;
        auto hdivboundright = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmeshflux,gel,index);
        hdivboundright->SetConnectIndex(0,celhdivc->ConnectIndex(3));
        gel->ResetReference();

        auto leftfluxside = intfluxside.Neighbour();
        auto gel1dleft = leftfluxside.Element();
        auto hdivboundleft = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmeshflux,gel1dleft,index);
        hdivboundleft->SetConnectIndex(0,celhdivc->ConnectIndex(4));
        gel1dleft->ResetReference();
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

    auto mat = new TPZNullMaterial(fDifpressure, dim, nstate);
    cmeshm.InsertMaterialObject(mat);

    for (auto matid : fCondensedMatids)    
    {
        auto mat = new TPZNullMaterial(matid, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }
    {
        auto mat = new TPZLagrangeMultiplier(fInterface, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }
    cmeshm.SetAllCreateFunctionsMultiphysicElem();
}

void TPZBuildSBFemHdiv::AddInterfaceElements(TPZMultiphysicsCompMesh & cmeshm, set<int> & matids1d)
{
    for (auto gel : fGMesh->ElementVec())
    {
        if (!gel)
        {
            continue;
        }
        auto matidskeleton = gel->MaterialId();
        auto it = matids1d.find(matidskeleton);
        if (it == matids1d.end())
        {
            continue;
        }
        auto side = GetSideSkeletonEl(gel);
        TPZGeoElSide gelside(gel,side); // gelside for the external pressure

        auto gelsideint = gelside.Neighbour();
        auto matidinterface = gelsideint.Element()->MaterialId();
        while (matidinterface != fInterface)
        {
            gelsideint = gelsideint.Neighbour();
            matidinterface = gelsideint.Element()->MaterialId();
        }

        auto gelsiderightf = gelsideint.Neighbour();
        auto matidrightflux = gelsiderightf.Element()->MaterialId();
        if (matidrightflux != fExternalrightflux)
        {
            DebugStop();
        }
        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celsiderightf = gelsiderightf.Reference();
        if (!celside || !celsiderightf) {
            DebugStop();
        }
        int64_t index;
        TPZMultiphysicsInterfaceElement *intl = new TPZMultiphysicsInterfaceElement(cmeshm, gelsideint.Element(),index,celsiderightf,celside);

        auto gelsideinternal = gelsiderightf.Neighbour(); // Next neighbour will be internal
        auto gelsideleftf = gelsideinternal.Neighbour(); // Next neighbour will be flux left
        if (!gelsideleftf || gelsideleftf.Element()->MaterialId() != fExternalleftflux)
        {
            DebugStop();
        }
        gelsideint = gelsideleftf.Neighbour(); // the neighbour of the right flux must be the interface again
        if (!gelsideint || gelsideint.Element()->MaterialId() != fInterface)
        {
            DebugStop();
        }
        auto gelsidepr = gelsideint.Neighbour(); // Next neighbour will be the internal pressure (dif. of pressure)
        if (!gelsidepr || gelsidepr.Element()->MaterialId() != fDifpressure)
        {
            DebugStop();
        }
        TPZCompElSide celsidepr = gelsidepr.Reference();
        TPZCompElSide celsideleftf = gelsideleftf.Reference();
        if (!celsidepr || !celsideleftf) {
            DebugStop();
        }
        TPZMultiphysicsInterfaceElement *intr = new TPZMultiphysicsInterfaceElement(cmeshm, gelsideint.Element(),index,celsideleftf,celsidepr);
        
    }
}

// Here SBFEM Multiphysics Volume elements are created.
// The Multiphysics element must be fSkeleton
// But I still need to pass the Volumetric Collapsed element so the TPZSBFemVolumeHdiv class will be able to perform integrations.
// Remember the order of the material IDs:
// (fSkeleton -> fElgroup -> fCollapsedMatId) -> fInterface -> fExtrightflux -> fInternal ->...
void TPZBuildSBFemHdiv::CreateSBFEMMultiphysicsVol(TPZMultiphysicsCompMesh & cmeshm, set<int> &matids, set<int> & matidtarget)
{
    // This for creates the volumetric multiphysics element, working as an AutoBuild()
    for (auto gel : fGMesh->ElementVec())
    {
        // Now I'm searching for the collapsed element, in which the matid is in matidtarget.
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

        // Adding multiphysics compels:
        // The neighbour will be in the surface identified as side. So I get the gelside and search for its neighbour.
        // I'll include all computational elements related to geo collapsed element
        auto side = GetSideCollapsedEl(gel);
        TPZGeoElSide gelside(gel, side);

        // 1st neighbour - Interface;
        // 2nd neighbour - fExtrightflux;
        // 3rd neighbour - fInternal;
        // 4th neighbour - fExtleftflux;
        // 5th neighbour - fInterface;
        // 6th neighbour - fDifPressure;
        // 7th neighbour - Skeleton.
        auto neigh = gelside.Neighbour();
        auto gel1d = neigh.Element();
        if(!gel1d || !(gel1d->Reference()))
        {
            DebugStop();
        }
        auto cel1d = gel1d->Reference();
        auto celsbfem = cmeshm.Element(index);
        auto sbfem = dynamic_cast<TPZSBFemVolumeHdiv * >(celsbfem);
        if(sbfem)
        {
            sbfem->AddElement(cel1d, 0);
        } else
        {
            DebugStop();
        }   
        for (auto i = 1; i < 7; i++)
        {
            neigh = neigh.Neighbour();
            gel1d = neigh.Element();
            if(!gel1d || !(gel1d->Reference()))
            {
                DebugStop();
            }
            cel1d = gel1d->Reference();
            celsbfem = cmeshm.Element(index);
            auto sbfem = dynamic_cast<TPZSBFemVolumeHdiv * >(celsbfem);
            if(sbfem)
            {
                sbfem->AddElement(cel1d, i);
            } else
            {
                DebugStop();
            }   
        }
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