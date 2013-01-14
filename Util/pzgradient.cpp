/*
 *  pzgradient.cpp
 *  PZ
 *
 *  Created by Agnaldo on 14/01/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <string>

#include "pzgradient.h"

TPZGradient::TPZGradient(): TPZFunction<STATE>(){
    
    fCenter.Resize(3,0.);
    fGradient.Redim(3,1);
    fUc=0.;
}

TPZGradient::TPZGradient(const TPZGradient &cp): TPZFunction<STATE>(cp), fCenter(cp.fCenter), fGradient(cp.fGradient){
    
    fUc = cp.fUc;
}

