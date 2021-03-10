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
    
    bool keepmatrix = false;
    fCondEl = new TPZCondensedCompEl(fCondensedEls, keepmatrix);    
}

void TPZSBFemMultiphysicsElGroup::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    TPZElementMatrix E0, E1, E2;
    ComputeMatrices(E0, E1, E2);
    
    InitializeElementMatrix(ek, ef);

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
    
    
    static pthread_mutex_t mutex_serial = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&mutex_serial);

    TPZFMatrix<STATE> globmatkeep(globmat);
    TPZFMatrix<complex<double> > eigenVectors;
    TPZManVector<complex<double> > eigenvalues;
    globmatkeep.SolveEigenProblem(eigenvalues, eigenVectors);
    cout << eigenvalues << "\n";
    
    pthread_mutex_unlock(&mutex_serial);
        
    TPZFNMatrix<200,std::complex<double> > QVectors(n,n,0.);
    fPhi.Resize(n, n);
    TPZManVector<std::complex<double> > eigvalsel(n,0);
    TPZFMatrix<std::complex<double> > eigvecsel(2*n,n,0.),eigvalmat(1,n,0.);
    int count = 0;
    for (int i=0; i<2*n; i++) {
        if (eigenvalues[i].real() < -1.e-6) {
            double maxvaleigenvec = 0;
            for (int j=0; j<n; j++) {
                QVectors(j,count) = eigenVectors(j+n,i);
                eigvecsel(j,count) = eigenVectors(j,i);
                eigvecsel(j+n,count) = eigenVectors(j+n,i);
                fPhi(j,count) = eigenVectors(j,i);
                double realvalabs = fabs(fPhi(j,count).real());
                if (realvalabs > maxvaleigenvec) {
                    maxvaleigenvec = realvalabs;
                }
            }
            eigvalsel[count] = eigenvalues[i];
            eigvalmat(0,count) = eigenvalues[i];
            for (int j=0; j<n; j++) {
                QVectors(j,count) /= maxvaleigenvec;
                eigvecsel(j,count) /= maxvaleigenvec;
                eigvecsel(j+n,count) /= maxvaleigenvec;
                fPhi(j,count) /= maxvaleigenvec;

            }
            count++;
        }
    }

    if (dim == 2)
    {
        int nstate = Connect(0).NState();
        if (nstate != 2 && nstate != 1) {
            DebugStop();
        }
        if(count != n-nstate) {
            DebugStop();
        }
        int ncon = fCondEl->NConnects();
        int eq=0;
        std::set<int64_t> cornercon;
        BuildCornerConnectList(cornercon);
        for (int ic=0; ic<ncon; ic++) {
            int64_t conindex = fCondEl->ConnectIndex(ic);
            if (cornercon.find(conindex) != cornercon.end())
            {
                fPhi(eq,count) = 1;
                eigvecsel(eq,count) = 1;
                if (nstate == 2)
                {
                    fPhi(eq+1,count+1) = 1;
                    eigvecsel(eq+1,count+1) = 1;
                }
            }
            eq += fCondEl->Connect(ic).NShape()*fCondEl->Connect(ic).NState();
        }
    }
    if(dim==3 && count != n)
    {
        DebugStop();
    }
    fEigenvalues = eigvalsel;
        
    TPZFMatrix<std::complex<double> > phicopy(fPhi);
    fPhiInverse.Redim(n, n);
    fPhiInverse.Identity();
    
    try
    {
        TPZVec<int> pivot;
        phicopy.Decompose_LU(pivot);
        phicopy.Substitution(&fPhiInverse, pivot);
    }
    catch(...)
    {
        exit(-1);
    }

    TPZFMatrix<std::complex<double> > ekloc;
    QVectors.Multiply(fPhiInverse, ekloc);
    if(0)
    {
        std::ofstream out("EigenProblem.nb");
        globmatkeep.Print("matrix = ",out,EMathematicaInput);
        eigvecsel.Print("eigvec =",out,EMathematicaInput);
        eigvalmat.Print("lambda =",out,EMathematicaInput);
        fPhi.Print("phi = ",out,EMathematicaInput);
        fPhiInverse.Print("phiinv = ",out,EMathematicaInput);
        QVectors.Print("qvec = ",out,EMathematicaInput);
    }

    ek.fMat.Resize(ekloc.Rows(),ekloc.Cols());
    
    TPZFMatrix<double> ekimag(ekloc.Rows(),ekloc.Cols());
    for (int i=0; i<ekloc.Rows(); i++) {
        for (int j=0; j<ekloc.Cols(); j++) {
            ek.fMat(i,j) = ekloc(i,j).real();
            ekimag(i,j) = ekloc(i,j).imag();
        }
    }
    
    for (auto cel : fElGroup)
    {
        TPZSBFemVolumeHdiv *sbfem = dynamic_cast<TPZSBFemVolumeHdiv *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->SetPhiEigVal(fPhi, fEigenvalues);
    }
    ofstream sout("cmeshmultiphysics.txt");
    Mesh()->Print(sout);
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

    auto ncon = fCondEl->NConnects();
    for (auto icon = 0; icon < ncon/2; icon++)
    {
        auto &c = fCondEl->Connect(icon);
        c.SetCondensed(true);
    }
    auto celref = fCondEl->ReferenceCompEl();
    auto elgr = dynamic_cast<TPZElementGroup * >(celref);
    if(!elgr)
    {
        DebugStop();
    }
    elgr->ReorderConnects();
    Mesh()->InitializeBlock();
    fCondEl->Resequence();
}

void TPZSBFemMultiphysicsElGroup::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef) const
{
	const int ncon = fCondEl->NConnects();
	int numeq = 0;
	ef.fBlock.SetNBlocks(ncon);
	ek.fBlock.SetNBlocks(ncon);
    
    for (int ic=0; ic<ncon; ic++)
    {
        TPZConnect &c = Connect(ic);
        int blsize = c.NShape()*c.NState();
        numeq += blsize;
        ef.fBlock.Set(ic, blsize);
        ek.fBlock.Set(ic, blsize);
    }
    const int numloadcases = 1;
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
	ef.fMat.Redim(numeq,numloadcases);
	ef.fConnect.Resize(ncon);

    ek.fMesh = Mesh();
    ek.fType = TPZElementMatrix::EK;
	ek.fMat.Redim(numeq,numeq);
	ek.fConnect.Resize(ncon);

	for(int i=0; i<ncon; i++)
    {
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	}
    std::map<int64_t,TPZOneShapeRestraint>::const_iterator it;
    for (it = fRestraints.begin(); it != fRestraints.end(); it++) 
    {
        ef.fOneRestraints.push_back(it->second);
        ek.fOneRestraints.push_back(it->second);
    }
}//void