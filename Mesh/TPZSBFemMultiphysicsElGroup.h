//
//  TPZSBFemElementGroup.hpp
//  PZ
//
//  Created by Karolinne Coelho on 22/01/21.
//
//

#pragma once

#include <stdio.h>
#include "TPZSBFemVolumeHdiv.h"
#include "TPZSBFemElementGroup.h"

using namespace std;

class TPZSBFemMultiphysicsElGroup : public TPZSBFemElementGroup
{
    
private:
    
    TPZElementGroup * fCondensedEls;
    
public:
    
    TPZSBFemMultiphysicsElGroup() : TPZSBFemElementGroup()
    {
        
    }
    
    /// constructor
    TPZSBFemMultiphysicsElGroup(TPZCompMesh &mesh, int64_t &index) : TPZSBFemElementGroup(mesh,index)
    {
        
    }

    void AddElement(TPZCompEl *cel);

    virtual void Print(std::ostream &out) const;

    void GroupandCondense(set<int> & matidscondensed);


    void ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2, TPZElementMatrix &M0);
    
    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);


    void ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2);
};