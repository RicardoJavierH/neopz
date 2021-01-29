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

void TPZBuildSBFemHdiv::BuildMultiphysicsCompMesh(TPZCompMesh &cmesh)
{
    // Getting the multiphysics mesh
    auto cmeshm = dynamic_cast<TPZMultiphysicsCompMesh * >(&cmesh);
    if (!cmeshm)
    {
        DebugStop();
    }
    
    TPZManVector<TPZCompMesh*, 2> cmeshvec = cmeshm->MeshVector();
    auto cmeshflux = cmeshvec[0];
    auto cmeshpressure = cmeshvec[1];

    // Creating the geometry for the the SBFEM simulation
    // Skeleton Elements + Collapsed Elements + External elements:

    // Creating dim-1 elements CompEls
    int dim = cmeshflux->Dimension();
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
    cmeshflux->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmeshflux->AutoBuild(matids);

    auto gmeshflux = cmeshflux->Reference();
    CreateExternalElements(gmeshflux);

    // Creating volumetric and external elements
    CreateVolumetricElementsHdiv(*cmeshflux); // Here it's created the SBFemVolumeHdiv elements

    fGMesh = gmeshflux;

#ifdef PZDEBUG
    ofstream outvtk("GeometrySBFEM.vtk");
    TPZVTKGeoMesh vtk;
    vtk.PrintGMeshVTK(fGMesh, outvtk, true);
    ofstream gout("buildsbfemgmesh.txt");
    cmeshflux->Reference()->Print(gout);
#endif

    // Create discontinuous comp elements for the external pressures
    CreateSBFemDiscontinuousElements(*cmeshpressure);

#ifdef PZDEBUG
    ofstream fout("cmeshflux.txt");
    cmeshflux->Print(fout);
    ofstream pout("cmeshpressure.txt");
    cmeshpressure->Print(pout);
#endif

    UpdateMultiphysicsMesh(cmeshvec, *cmeshm);

    CreateSBFemMultiphysicsElGroups(*cmeshm);

#ifdef PZDEBUG
    ofstream mout("cmeshmultiphysics.txt");
    cmeshm->Print(mout);
#endif

    // Creating Element Groups
    CreateSBFemInterfaceElementGroups(*cmeshm);

    // Ajusting the connectivity of these elements
    // AdjustExternalPressureConnectivity();

    // Group and Condense
    // GroupandCondense();
}

void TPZBuildSBFemHdiv::CreateExternalElements(TPZGeoMesh * gmesh)
{
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel)
        {
            continue;
        }
        if (gel->MaterialId() != fSkeletonMatId)
        {
            continue;
        }
        auto iside = -1;
        switch (gel->Type())
        {
        case EOned:
            iside = 3;
            break;
        case ETriangle:
            iside = 6;
            break;
        case EQuadrilateral:
            iside = 8;
        default:
            break;
        }
        TPZGeoElBC(gel, iside, fLeftpressureMatId);
        TPZGeoElBC(gel, iside, fRightpressureMatId);
        TPZGeoElBC(gel, iside, fLeftfluxMatId);
        TPZGeoElBC(gel, iside, fRightfluxMatId);
    }
}

void TPZBuildSBFemHdiv::CreateSBFemDiscontinuousElements(TPZCompMesh &cmeshpressure)
{
    auto dim = cmeshpressure.Dimension()-1;
    auto nstate = cmeshpressure.MaterialVec().begin()->second->NStateVariables();

    auto matleft = new TPZNullMaterial(fLeftpressureMatId, dim, nstate);
    cmeshpressure.InsertMaterialObject(matleft);

    auto matright = new TPZNullMaterial(fRightpressureMatId, dim, nstate);
    cmeshpressure.InsertMaterialObject(matright);

    set<int> matids = {fLeftpressureMatId, fRightpressureMatId};
    cmeshpressure.SetAllCreateFunctionsContinuous();
    cmeshpressure.ApproxSpace().CreateDisconnectedElements(true);
    cmeshpressure.AutoBuild(matids);

    for(auto newnod : cmeshpressure.ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }
}

void TPZBuildSBFemHdiv::UpdateMultiphysicsMesh(TPZManVector<TPZCompMesh*, 2> & cmeshvec, TPZMultiphysicsCompMesh & cmeshm)
{
    auto dim = fGMesh->Dimension();
    auto nstate = cmeshm.MaterialVec().begin()->second->NStateVariables();
    cmeshm.ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    set<int> matids = {fLeftfluxMatId, fRightpressureMatId, fLeftpressureMatId, fRightpressureMatId};
    for(auto mId : matids)
    {
        auto mat = new TPZNullMaterial(mId,dim,nstate);
        cmeshm.InsertMaterialObject(mat);
    }
    {
        auto mat = new TPZLagrangeMultiplier(fInterfaceMatId, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }

    TPZManVector<int> active(2,1);
    cmeshm.BuildMultiphysicsSpace(active, cmeshvec);
    cmeshm.LoadReferences();
    cmeshm.CleanUpUnconnectedNodes();
}

void TPZBuildSBFemHdiv::CreateVolumetricElementsHdiv(TPZCompMesh &cmesh)
{
    TPZGeoMesh *gmesh = cmesh.Reference();

    // all computational elements have been loaded
    set<int> matids, matidstarget;
    for (auto it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++)
    {
        int64_t mat = it->second;
        if (cmesh.FindMaterial(mat))
        {
            matids.insert(it->first);
            matidstarget.insert(it->second);
        }
    }

    // Creating geometric collapsed elements
    auto dim = gmesh->Dimension();
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->HasSubElement() || gel->Reference())
        {
            continue;
        }
        // we create SBFemVolume elements by partitioning the volume elements
        auto el = gel->Index();
        if (gel->Dimension() != dim || fElementPartition[el] == -1)
        {
            continue;
        }
        if (matids.find(gel->MaterialId()) == matids.end())
        {
            continue;
        }
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++)
        {
            if (gel->SideDimension(is) != dim-1)
            {
                continue;
            }
            TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gel,is);
            bool onlyinterpolated = false; // TEST HERE!! - I guess it should be false because I am using Hdiv approx space
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
                auto matid = fMatIdTranslation[gel->MaterialId()];
                if (subgelside.IsLinearMapping())
                {
                    switch(nnodes)
                    {
                        case 2:
                            gmesh->CreateGeoElement(EQuadrilateral, Nodes, matid, index);
                            break;
                        case 4:
                            gmesh->CreateGeoElement(ECube, Nodes, matid, index);
                            DebugStop();
                            break;
                        case 3:
                            gmesh->CreateGeoElement(EPrisma, Nodes, matid, index);
                            DebugStop();
                            break;
                        default:
                            std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                            DebugStop();
                    }
                    
                }
                else
                {
                    int64_t elementid = gmesh->NElements()+1;
                    switch(nnodes)
                    {
                        case 2:
                            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (Nodes, matid, *gmesh,index);
                            break;
                        case 4:
                            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > (Nodes, matid, *gmesh,index);
                            break;
                        case 3:
                            gmesh->CreateGeoElement(EPrisma, Nodes, matid, index);
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
    gmesh->BuildConnectivity();

#ifdef PZDEBUG
    {
        ofstream gout("gmeshwithvol.txt");
        gmesh->Print(gout);
        ofstream cout("cmeshhdiv.txt");
        cmesh.Print(cout);
    }
#endif

    // Creating SBFemVolumeHdiv elements
    gmesh->ResetReference();
    CreateSBFemVolumeHdiv(cmesh, matidstarget);
    cmesh.ExpandSolution();

    AdjustFluxConnectivities(cmesh);

#ifdef PZDEBUG
    {
        ofstream cout("cmeshhdiv.txt");
        cmesh.Print(cout);
    }
#endif
}


void TPZBuildSBFemHdiv::CreateSBFemVolumeHdiv(TPZCompMesh & cmesh, set<int> & matidstarget)
{
    auto gmesh = cmesh.Reference();
    auto dim = gmesh->Dimension();
    auto porder = cmesh.GetDefaultOrder();

    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->Dimension() != gmesh->Dimension())
        {
            continue;
        }
        auto it = matidstarget.find(gel->MaterialId());
        if (it == matidstarget.end())
        {
            continue;
        }
        // Get Side of the volumetric element with the dim-1 element
        auto nnodes = gel->NNodes();
        auto iside = -1;
        switch (nnodes)
        {
        case 4:
            iside = 4;
            break;
        case 6:
            DebugStop();
            break;
        case 8:
            DebugStop();
            break;
        default:
            DebugStop();
            break;
        }
        
        int64_t index;
        TPZGeoElSide gelside(gel,iside);
        auto neigh  = gelside.Neighbour();
        while (neigh != gelside)
        {
            auto gelneigh = neigh.Element();
            auto matid = gelneigh->MaterialId();
            if (gelneigh->Dimension() == dim-1 && matid != fLeftfluxMatId && matid != fRightfluxMatId)
            {
                break;
            }
            neigh = neigh.Neighbour();
        }
        if (neigh == gelside)
        {
            DebugStop();
        }
        auto gel1d = neigh.Element();

        auto celhdivc = CreateSBFemHdivCompEl(cmesh, gel, gel1d, index);
        fGeltocel[gel1d->Index()] = celhdivc;
    }
}



void TPZBuildSBFemHdiv::CreateSBFemMultiphysicsElGroups(TPZMultiphysicsCompMesh & cmesh)
{
    auto numgroups = fPartitionCenterNode.size();
    auto groupelementindices(numgroups);
    
    TPZManVector<int64_t> elementgroupindices(numgroups); 
    
    for (int64_t el=0; el<numgroups; el++)
    {
        int64_t index;
        new TPZSBFemMultiphysicsElGroup(cmesh,index);
        elementgroupindices[el] = index;
    }
    
    auto dim = cmesh.Dimension();
    for (auto cel : cmesh.ElementVec())
    {
        if (!cel || !(cel->Reference()) )
        {
            continue;
        }
        auto gel = cel->Reference();
        auto sbfem = dynamic_cast<TPZSBFemVolumeHdiv * >(cel);
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
        if (sbfem)
        {
            auto gel = sbfem->Reference();
            auto gelindex = gel->Index();

            TPZGeoElSide gelside(gel,side);
            auto geldim = gel->Dimension();
            auto nsidenodes = gel->NSideNodes(side);

            TPZManVector<int64_t,8> cornernodes(nsidenodes);
            for (int node = 0; node<nsidenodes; node++)
            {
                cornernodes[node] = gel->SideNodeIndex(side, node);
            }
            
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside)
            {
                if(neighbour.Element()->Dimension() == geldim-1 && neighbour.Element()->Reference())
                {
                    int nsidenodesneighbour = neighbour.Element()->NCornerNodes();
                    if (nsidenodesneighbour == nsidenodes)
                    {
                        TPZManVector<int64_t,8> neighnodes(nsidenodesneighbour);
                        for (int i=0; i<nsidenodesneighbour; i++)
                        {
                            neighnodes[i] = neighbour.Element()->NodeIndex(i);
                        }
                        int numequal = 0;
                        for (int i=0; i<nsidenodesneighbour; i++)
                        {
                            if (neighnodes[i] == cornernodes[i]) {
                                numequal++;
                            }
                        }
                        if (numequal == nsidenodesneighbour)
                        {
                            break;
                        }
                    }
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside)
            {
                // we are not handling open sides (yet)
                DebugStop();
            }
            int64_t skelindex = neighbour.Element()->Reference()->Index();
            sbfem->SetSkeleton(skelindex);
            
            int64_t gelgroup = fElementPartition[gelindex];
            if (gelgroup == -1)
            {
                DebugStop();
            }
            int64_t celgroupindex = elementgroupindices[gelgroup];
            TPZCompEl *celgr = cmesh.Element(celgroupindex);
            TPZSBFemMultiphysicsElGroup *sbfemgr = dynamic_cast<TPZSBFemMultiphysicsElGroup *>(celgr);
            if (!sbfemgr)
            {
                DebugStop();
            }
            sbfemgr->AddElement(sbfem);
            sbfem->SetElementGroupIndex(celgroupindex);
        }
    }
    
    // for (int64_t el=0; el<numgroups; el++) {
    //     int64_t index;
        
    //     index = elementgroupindices[el];
    //     TPZCompEl *cel = cmesh.Element(index);
    //     TPZSBFemElementGroup *sbfemgroup = dynamic_cast<TPZSBFemElementGroup *>(cel);
    //     if (!sbfemgroup) {
    //         DebugStop();
    //     }
    //     const TPZVec<TPZCompEl *> &subgr = sbfemgroup->GetElGroup();
    //     int64_t nsub = subgr.NElements();
    //     for (int64_t is=0; is<nsub; is++) {
    //         TPZCompEl *cel = subgr[is];
    //         TPZSBFemVolume *femvol = dynamic_cast<TPZSBFemVolume *>(cel);
    //         if (!femvol) {
    //             DebugStop();
    //         }
    //         femvol->SetElementGroupIndex(index);
    //     }
    // }
}

// NOT READY YET

void TPZBuildSBFemHdiv::AdjustFluxConnectivities(TPZCompMesh & cmesh)
{
    auto gmesh = cmesh.Reference();
    gmesh->ResetReference();

    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->MaterialId() != fLeftfluxMatId) continue;

        auto nnodes1d = gel->NNodes();
        auto iside = -1;
        switch (nnodes1d)
        {
        case 2:
            iside = 2;
            break;
        
        default:
            break;
        }

        TPZGeoElSide gelside(gel,iside);
        auto intfluxside = gelside.HasNeighbour(fSkeletonMatId);
        if (!intfluxside)
        {
            continue;
        }
        auto cel = fGeltocel[intfluxside.Element()->Index()];
        auto celhdivc = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * >(cel);
        if (!celhdivc)
        {
            DebugStop();
        }

        int64_t index;
        auto hdivboundleft = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmesh,gel,index);
        hdivboundleft->SetConnectIndex(0,celhdivc->ConnectIndex(3));
        gel->ResetReference();

        auto rightfluxside = gelside.HasNeighbour(fRightfluxMatId);
        auto gel1dright = rightfluxside.Element();
        auto hdivboundright = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmesh,gel1dright,index);
        hdivboundright->SetConnectIndex(0,celhdivc->ConnectIndex(4));
        gel1dright->ResetReference();
    }
    
    cmesh.LoadReferences();
    cmesh.Reference()->ResetReference();
    cmesh.CleanUpUnconnectedNodes();
    cmesh.ExpandSolution();
}

void TPZBuildSBFemHdiv::CreateSBFemInterfaceElementGroups(TPZCompMesh & cmeshm)
{
    auto gmesh = cmeshm.Reference();
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != fLeftfluxMatId)
        {
            continue;
        }
        auto nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        auto gelsidepr = gelside.HasNeighbour(fLeftpressureMatId);
        if (!gelsidepr)
        {
            DebugStop();
        }

        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celneigh = gelsidepr.Reference();
        if (!celside || !celneigh) {
            DebugStop();
        }
        TPZGeoElBC gelbc(gelside, fInterfaceMatId);
        int64_t index;
        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(cmeshm,gelbc.CreatedElement(),index,celneigh,celside);
    }

    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != fRightfluxMatId)
        {
            continue;
        }
        auto nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        auto gelsidepr = gelside.HasNeighbour(fRightpressureMatId);
        if (!gelsidepr)
        {
            DebugStop();
        }

        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celneigh = gelsidepr.Reference();
        if (!celside || !celneigh) {
            DebugStop();
        }
        TPZGeoElBC gelbc(gelside, fInterfaceMatId);
        int64_t index;
        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(cmeshm,gelbc.CreatedElement(),index,celneigh,celside);
    }
}

void TPZBuildSBFemHdiv::GroupandCondense()
{
    DebugStop();
}

void TPZBuildSBFemHdiv::BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh)
{
    DebugStop();
}