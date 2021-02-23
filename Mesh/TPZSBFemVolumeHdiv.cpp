
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
#include "TPZGeoLinear.h"
#include "pzmultiphysicscompel.h"

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

TPZElementMatrix TPZSBFemVolumeHdiv::ComputeEKlocal()
{
    // ekloc and efloc will be initialized with the connectivity related to the multiphysics element
    // the InitializeElementMatrix is performed inside CalcStiff from TPZMultiphysicsCompEl
    TPZElementMatrix ekloc(Mesh(), TPZElementMatrix::EK);
    TPZElementMatrix efloc(Mesh(), TPZElementMatrix::EF);

    // ekvol and efvol will be initialized with the connectivity of the TPZSBFemVolumeHdiv
    // which means that this ekvol is composed by the contributions of ekloc computed for each element in fElementVec.
    TPZElementMatrix ekvol(Mesh(), TPZElementMatrix::EK);
    TPZElementMatrix efvol(Mesh(), TPZElementMatrix::EF);
    InitializeElementMatrix(ekvol, efvol);

    map<int64_t,int64_t> locindex;
    auto ncon = fConnectIndexes.size();
    for (int64_t ic=0; ic<ncon ; ic++)
    {
        locindex[fConnectIndexes[ic]] = ic;
    }
    
    for (int iel=0; iel<7; iel++)
    {
        if(iel == 0 || iel == 1 || iel == 3 || iel == 4) continue;
        auto cel = fElementVec[iel];
#ifdef PZDEBUG
        // continue if it's a interface element
        if(!cel || !(cel->Reference())) continue;
#endif
        switch (cel->Reference()->Type())
        {
        case EOned:
        {
            auto celm = dynamic_cast<TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear> *>(cel);
            celm->CalcStiff(ekloc, efloc);
        }    
            break;
        default:
            break;
        }

        // Contributing ekloc into ekvol
        int nelcon = ekloc.NConnects();
        for (int ic=0; ic<nelcon; ic++)
        {
            int iblsize = ekloc.fBlock.Size(ic);
            int icindex = ekloc.fConnect[ic];
            int ibldest = locindex[icindex];
            for (int jc = 0; jc<nelcon; jc++)
            {
                int jblsize = ekloc.fBlock.Size(jc);
                int jcindex = ekloc.fConnect[jc];
                int jbldest = locindex[jcindex];
                for (int idf = 0; idf<iblsize; idf++)
                {
                    for (int jdf=0; jdf<jblsize; jdf++)
                    {
                        ekvol.fBlock(ibldest,jbldest,idf,jdf) += ekloc.fBlock(ic,jc,idf,jdf);
                    }
                }
            }
        }
    }
    return ekvol;
}

// Not sure how if this method is right
void TPZSBFemVolumeHdiv::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
    const int ncon = this->NConnects();
    int numeq = 0;
    int ic;
    
    for(ic=0; ic<ncon; ic++)
    {
        int64_t neqThisConn = Connect(ic).NDof(*Mesh());
        numeq += neqThisConn;
    }
    
    int64_t nref = this->fElementVec.size();
    int nstate = 0;
    int numloadcases = 1;
    
    ek.fMat.Redim(numeq,numeq);
    ef.fMat.Redim(numeq,numloadcases);
    ek.fBlock.SetNBlocks(ncon);
    ef.fBlock.SetNBlocks(ncon);
    
    int i;
    for(i=0; i<ncon; i++){
        unsigned int ndof = Connect(i).NDof(*Mesh());
#ifdef PZDEBUG
        TPZConnect &c = Connect(i);
        if (c.NShape()*c.NState() != ndof) {
            DebugStop();
        }
#endif
        ek.fBlock.Set(i,ndof);
        ef.fBlock.Set(i,ndof);
    }
    ek.fConnect.Resize(ncon);
    ef.fConnect.Resize(ncon);
    for(i=0; i<ncon; i++){
        (ek.fConnect)[i] = ConnectIndex(i);
        (ef.fConnect)[i] = ConnectIndex(i);
    }
    ek.fOneRestraints = GetShapeRestraints();
    ef.fOneRestraints = GetShapeRestraints();
}

TPZCompEl * CreateSBFemMultiphysicsCompEl(TPZMultiphysicsCompMesh &mesh, TPZGeoEl *gel, int64_t &index)
{
    new TPZSBFemVolumeHdiv(mesh, gel, index);    
}