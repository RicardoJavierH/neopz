//
//  TPZGeoMeshBluider.h
//  pz
//
//  Created by Omar Durán on 2/12/19.
//

#ifndef TPZGeoMeshBluider_h
#define TPZGeoMeshBluider_h

#include <stdio.h>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"


#include "tpzpoint.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "tpzcube.h"
#include "pzgeopyramid.h"


#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticcube.h"
#include "tpzquadratictetra.h"
#include "tpzquadraticprism.h"
#include "tpzgeoblend.h"

#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include <tpzarc3d.h>

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoElement.h"

class TPZGeoMeshBluider {
    
public:
    
    static void InsertNodes(TPZGeoMesh * gmesh, std::vector<int> & node_identifiers, std::vector<double> & coord);
    
    static void InsertElement(TPZGeoMesh * gmesh, int & physical_identifier, int & el_type, int & el_identifier, std::vector<int> & node_identifiers);
    
    static int GetNumberofNodes(int & el_type);
    
    static void PrintGeometry(TPZGeoMesh * gmesh, std::string  & name);
    
};

#endif /* TPZGeoMeshBluider_h */
