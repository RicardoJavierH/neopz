//
//  TPZSBFemElementGroup.cpp
//  PZ
//
//  Created by Karolinne Coelho on 22/01/2021.
//
//
#include "TPZSBFemMultiphysicsElGroup.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzmaterialdata.h"

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


void TPZSBFemMultiphysicsElGroup::GroupandCondense(set<int> & condensedmatid)
{
#ifdef PZDEBUG
    if (fElGroup.size() == 0 || condensedmatid.size() != 3)
    {
        DebugStop();
    }
#endif
    int64_t index;
    fCondensedEls = new TPZElementGroup(*Mesh(), index);

    for (auto celvol : fElGroup)
    {
        if (!celvol)
        {
            continue;
        }
        auto sbfemvol = dynamic_cast<TPZSBFemVolumeHdiv * >(celvol);

#ifdef PZDEBUG
        if (!sbfemvol)
        {
            DebugStop();
        }
#endif
        for (auto cel : sbfemvol->ElementVec())
        {
            if (!cel)
            {
                DebugStop();
            }
            auto matid = cel->Reference()->MaterialId();
            auto it = condensedmatid.find(matid);
            if (it == condensedmatid.end())
            {
                continue;
            }
            fCondensedEls->AddElement(cel);
        }
    }
    // comp el pressure -> set node el con

    // Mesh()->ComputeNodElCon();

    bool keepmatrix = false;
    fCondEl = new TPZCondensedCompEl(fCondensedEls, keepmatrix);

    // verificacao: imprimir em um log, pegar connectindexes do cond
    // pegar ultima versao no refactorgeom
}

void TPZSBFemMultiphysicsElGroup::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    InitializeElementMatrix(ek, ef);

    TPZElementMatrix E0, E1, E2;
    ComputeMatrices(E0, E1, E2);

}

void TPZSBFemMultiphysicsElGroup::ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2)
{
    TPZElementMatrix ek(Mesh(),TPZElementMatrix::EK);
    TPZElementMatrix ef(Mesh(),TPZElementMatrix::EF);
    fCondEl->CalcStiff(ek,ef);
    
    auto n = ek.fMat.Rows()/2;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            E0.fMat(i,j) = ek.fMat(i,j)*1/4;
            E1.fMat(i,j) = ek.fMat(i+n,j)*1/2;
            E2.fMat(i,j) = ek.fMat(i+n,j+n);
        }
    }

}