
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

TPZSBFemVolumeHdiv::TPZSBFemVolumeHdiv(TPZMultiphysicsCompMesh & mesh, TPZGeoEl * gel, int64_t & index) : TPZInterpolationSpace(mesh, gel, index)
{
    fElementVec.Resize(7);
}

// NEED TO CHECK IT
void TPZSBFemVolumeHdiv::SetElementGroupIndex(int64_t index)
{
    fElementGroupIndex = index;
    std::map<int64_t, int> globtolocal;
    TPZCompEl *celgr = Mesh()->Element(index);
    fElementGroup = celgr;
    int nc = celgr->NConnects();
    TPZManVector<int, 10> firsteq(nc + 1, 0);
    for (int ic = 0; ic < nc; ic++)
    {
        globtolocal[celgr->ConnectIndex(ic)] = ic;
        TPZConnect &c = celgr->Connect(ic);
        firsteq[ic + 1] = firsteq[ic] + c.NShape() * c.NState();
    }
    int neq = 0;
    TPZCompEl *celskeleton = Mesh()->Element(fSkeleton);
    nc = celskeleton->NConnects();
    for (int ic = 0; ic < nc; ic++)
    {
        TPZConnect &c = celskeleton->Connect(ic);
        neq += c.NShape() * c.NState();
    }
    fLocalIndices.Resize(neq);
    int count = 0;
    for (int ic = 0; ic < nc; ic++)
    {
        int64_t cindex = celskeleton->ConnectIndex(ic);
#ifdef PZDEBUG
        if (globtolocal.find(cindex) == globtolocal.end()) {
            DebugStop();
        }
#endif
        TPZConnect &c = celskeleton->Connect(ic);
        int neq = c.NShape() * c.NState();
        int locfirst = firsteq[globtolocal[cindex]];
        for (int eq = 0; eq < neq; eq++)
        {
            fLocalIndices[count++] = locfirst + eq;
        }
    }
#ifdef PZDEBUG
    if (count != neq) DebugStop();
#endif
}

TPZCompEl * CreateSBFemMultiphysicsCompEl(TPZMultiphysicsCompMesh &mesh, TPZGeoEl *gel, int64_t &index)
{
    new TPZSBFemVolumeHdiv(mesh, gel, index);    
}