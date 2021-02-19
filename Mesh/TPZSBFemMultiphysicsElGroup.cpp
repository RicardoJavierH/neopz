//
//  TPZSBFemElementGroup.cpp
//  PZ
//
//  Created by Karolinne Coelho on 22/01/2021.
//
//
#include "TPZSBFemMultiphysicsElGroup.h"
#include "TPZMultiphysicsCompMesh.h"

void TPZSBFemMultiphysicsElGroup::AddElement(TPZCompEl *cel)
{
    std::set<int64_t> connects;
    for (auto ic : fConnectIndexes)
    {
        connects.insert(ic);
    }
    
    auto celvol = dynamic_cast<TPZSBFemVolumeHdiv *>(cel);
    TPZCompEl *celskeleton = Mesh()->Element(celvol->SkeletonIndex());

    auto nc = celskeleton->NConnects();
    for (int ic=0; ic<nc; ic++)
    {
        connects.insert(celskeleton->ConnectIndex(ic));
    }

    nc = connects.size();
    if (nc != fConnectIndexes.size())
    {
        fConnectIndexes.Resize(nc, 0);
        auto it = connects.begin();
        for (int ic = 0; it != connects.end(); it++,ic++)
        {
            fConnectIndexes[ic] = *it;
        }
    }
    TPZElementGroup::AddElement(cel);
}

void TPZSBFemMultiphysicsElGroup::Print(std::ostream &out) const
    {
        out << __PRETTY_FUNCTION__ << std::endl;
        TPZElementGroup::Print(out);

        out << "Element indexes of the volume elements ";
        for (auto cel : fElGroup)
        {
            out << cel->Index() << " ";
        }
        out << std::endl;
        out << "Indices of the associated computational skeleton elements\n";
        for (auto cel : fElGroup)
        {
            auto vol = dynamic_cast<TPZSBFemVolumeHdiv *>(cel);
            if(!vol) DebugStop();
            out << vol->SkeletonIndex() << " ";
        }
        out << std::endl;
        out << "Connect indexes of the contained elements\n";
        for (auto cel : fElGroup)
        {
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++)
            {
                out << cel->ConnectIndex(ic) << " ";
            }
            out << std::endl;
        }
        out << "End of " << __PRETTY_FUNCTION__ << std::endl;
    }


// void TPZSBFemMultiphysicsElGroup::GroupandCondense(TPZMultiphysicsCompMesh * cmesh)
// {
//     for (auto cel : cmesh->ElementVec())
//     {
//         if (!cel)
//         {
//             continue;
//         }
//         std::cout << "Element index " << cel->Index() << " ";
//         if (cel->Reference())
//         {
//             std::cout << "matid " << cel->Reference()->MaterialId();
//             if (cel->Reference()->MaterialId() == Eleftpressure || cel->Reference()->MaterialId() == Erightpressure)
//             {
//                 std::cout << " not added\n";
//                 continue;
//             }
//         }
//         if(cel == elgr)
//         {
//             std::cout << " group not added\n";
//             continue;
//         }
//         elgr->AddElement(cel);

//     }
//     cmesh->ComputeNodElCon();
    
//     {
//         std::ofstream out("cmesh.txt");
//         cmesh->Print(out);
//     }
//     bool keepmatrix = false;
//     auto cond = new TPZCondensedCompEl(elgr, keepmatrix);
//     cmesh->ExpandSolution();
// }