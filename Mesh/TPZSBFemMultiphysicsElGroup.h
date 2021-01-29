//
//  TPZSBFemElementGroup.hpp
//  PZ
//
//  Created by Karolinne Coelho on 22/01/21.
//
//

#pragma once

#include <stdio.h>

#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"


class TPZSBFemMultiphysicsElGroup : public TPZSBFemElementGroup
{
    
private:
    

    
public:
    
    TPZSBFemMultiphysicsElGroup() : TPZSBFemElementGroup()
    {
        
    }
    
    /// constructor
    TPZSBFemMultiphysicsElGroup(TPZCompMesh &mesh, int64_t &index) : TPZSBFemElementGroup(mesh,index)
    {
        
    }

    // void GroupandCondense(TPZMultiphysicsCompMesh * cmesh);
};