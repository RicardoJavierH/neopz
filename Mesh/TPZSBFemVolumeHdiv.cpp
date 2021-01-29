
//
//  TPZSBFemVolumeHdiv.cpp
//  PZ
//
//  Created by Karolinne Coelho on 25/01/2021.
//
//

#include "TPZSBFemVolumeHdiv.h"
#include "pzgeoelside.h"
#include "TPZSBFemElementGroup.h"
#include "pzintel.h"
#include "TPZMaterial.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzcmesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.sbfemvolume"));
#endif

TPZSBFemVolumeHdiv::TPZSBFemVolumeHdiv(TPZCompMesh &mesh, TPZGeoEl * gel, TPZGeoEl * gel1d, int64_t &index)
{
    switch (gel1d->Type())
    {
    case EOned:
    {
        auto iside = 4;
        TPZGeoElSide gelside(gel,iside);
        new TPZCompElHDivSBFem<TPZShapeLinear>(mesh, gel1d, gelside, index);
    }
        break;
    case ETriangle:
    {
        auto iside = 15;
        TPZGeoElSide gelside(gel,iside);
        new TPZCompElHDivSBFem<TPZShapeTriang>(mesh, gel1d, gelside, index);
    }
        break;
    case EQuadrilateral:
    {
        auto iside = 20;
        TPZGeoElSide gelside(gel,iside);
        new TPZCompElHDivSBFem<TPZShapeQuad>(mesh, gel1d, gelside, index);
    }
        break;
    default:
        DebugStop();
        break;
    }
}

TPZCompEl * CreateSBFemHdivCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, TPZGeoEl * gel1d, int64_t &index)
{
    new TPZSBFemVolumeHdiv(mesh, gel, gel1d, index);    
}