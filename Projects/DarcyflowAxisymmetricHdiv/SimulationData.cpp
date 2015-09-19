/*
 *  problemdata.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "SimulationData.h"

SimulationData::SimulationData()
{
    /** @brief Delta t - s */
    fDeltaT = 0.0;
    
    /** @brief Max Time - s */
    fMaxTime = 0.0;    
    
    /** @brief Time - s */
    fTime = 0.0;
    
    /** @brief Maximum number of newton iterations */
    fMaxiterations = 0;
    
    /** @brief Maximum number of newton iterations */
    fFixedJacobianIterations = 0;
    
    /** @brief DeltaX tolerance for newton iterations */
    ftoleranceDeltaX = 0.0;
    
    /** @brief Residual tolerance for newton iterations */
    ftoleranceResiual = 0.0;
    
    /** @brief Number of uniform mesh refinement */
    fHref = 0;
    
    /** @brief Number of uniform mesh refinement in postprocessing */
    fHrefpost = 0;
    
    /** @brief Approximation Order for velocity */
    fqorder = 1;
    
    /** @brief Approximation Order for pressure */
    fporder = 1;
    
    /** @brief Approximation Order for saturations */
    fsorder = 0;
    
    /** @brief Used a direct Solver */
    fIsDirect = false;
    
    /** @brief Optimize band width */
    fOptband = false;
    
    /** @brief Use Conjugated Gradient method */
    fIsCG = false;
    
    /** @brief Broyden iterations */
    fIsBroyden = false;
    
    /** @brief Define the use of static condensation */
    fCondenseElements = false;
    
    /** @brief Define the use dimensionless formulation */
    fIsDimensionless = false;
    
    /** @brief Define the use of a linear approxiamtion of S using a gradient reconstruction procedure */
    fUseGR = false;
    
    /** @brief void material being used for GR */
    fMatL2 = 1000;
    
    /** @brief State: n or n+1 temporal state */
    fnStep = true;
    
    /** @brief Definition of the Top bc */
    fTopBCini.Resize(4,0.0);
    
    /** @brief Definition of the Bottom bc */
    fBottomBCini.Resize(4,0.0);
    
    /** @brief Definition of the Right bc */
    fRightBCini.Resize(4,0.0);
    
    /** @brief Definition of the Left bc */
    fLeftBCini.Resize(4,0.0);
  
    /** @brief Definition of the Top bc */
    fTopBC.Resize(4,0.0);
    
    /** @brief Definition of the Bottom bc */
    fBottomBC.Resize(4,0.0);
    
    /** @brief Definition of the Right bc */
    fRightBC.Resize(4,0.0);
    
    /** @brief Definition of the Left bc */
    fLeftBC.Resize(4,0.0);
  
    /** @brief Gravity */
    fGravity.Resize(2,1);
//    fGravity.Zero();
    
    
}

SimulationData::~SimulationData()
{
    
}