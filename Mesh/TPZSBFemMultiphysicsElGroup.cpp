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
    if (fElGroup.size() == 0 || condensedmatid.size() == 0)
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
    auto ncon = fCondensedEls->NConnects();
    for (auto icon = 0; icon < ncon; icon++)
    {
        auto &c = fCondensedEls->Connect(icon);
        c.SetCondensed(true);
    }
    for (auto cel : fElGroup)
    {
        auto sbfem = dynamic_cast <TPZSBFemVolumeHdiv*>(cel);
        for (auto sbfemel : sbfem->ElementVec())
        {
            auto matid = sbfemel->Reference()->MaterialId();
            auto it = condensedmatid.find(matid);
            if(it == condensedmatid.end())
            {
                auto ncon = sbfemel->NConnects();
                for (auto icon = 0; icon < ncon; icon++)
                {
                    auto &c = sbfemel->Connect(icon);
                    c.SetCondensed(false);
                }
            }
        }
        
    }
    
    
    // comp el pressure -> set node el con
    fCondensedEls->Print(std::cout);
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

    int n = E0.fMat.Rows();
    auto dim = Mesh()->Dimension();
    
    TPZFMatrix<STATE> E0Inv(E0.fMat);

    TPZVec<int> pivot(E0Inv.Rows(),0);
    int nwork = 4*n*n + 2*n;
    TPZVec<STATE> work(2*nwork,0.);
    int info=0;
#ifdef STATEdouble
    dgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info);
#endif
#ifdef STATEfloat
    sgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info);
#endif
    if (info != 0) {
        DebugStop();
    }
#ifdef STATEdouble
    dgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info);
#endif
#ifdef STATEfloat
    sgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info);
#endif
    if (info != 0) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> globmat(2*n,2*n,0.);
    
#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n);
#else
    cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmat(i,j+n) = -E0Inv(i,j);
        }
    }
#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., &E1.fMat(0,0), n, &globmat(0,0), 2*n, 0., &globmat(n,0), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., &E1.fMat(0,0), n, &globmat(0,0), 2*n, 0., &globmat(n,0), 2*n);
#else
    cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmat(i+n,j) -= E2.fMat(i,j);
        }
    }

#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1., &E1.fMat(0,0), n, &E0Inv(0,0), n, 0., &globmat(n,n), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1., &E1.fMat(0,0), n, &E0Inv(0,0), n, 0., &globmat(n,n), 2*n);
#else
    cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif

    for (int i=0; i<n; i++) {
        globmat(i,i) -= (dim-2)*0.5;
        globmat(i+n,i+n) += (dim-2)*0.5;
    }
    
    TPZFMatrix<STATE> globmatkeep(globmat);
    TPZFMatrix<complex<double> > eigenVectors;
    TPZManVector<complex<double> > eigenvalues;
    globmatkeep.SolveEigenProblem(eigenvalues, eigenVectors);

    cout << eigenvalues << "\n";


}

void TPZSBFemMultiphysicsElGroup::ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2)
{
    TPZElementMatrix ek(Mesh(),TPZElementMatrix::EK);
    TPZElementMatrix ef(Mesh(),TPZElementMatrix::EF);
    fCondEl->CalcStiff(ek,ef);
    
    auto n = ek.fMat.Rows()/2;
    E0.fMat.Resize(n,n);
    E1.fMat.Resize(n,n);
    E2.fMat.Resize(n,n);
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