//
//  TPZSBFemElementGroup.cpp
//  PZ
//
//  Created by Karolinne Coelho on 22/01/2021.
//
//
#include "TPZSBFemMultiphysicsElGroup.h"
#include "TPZMultiphysicsCompMesh.h"

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