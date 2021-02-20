//
//  TPZSBFemElementGroup.cpp
//  PZ
//
//  Created by Karolinne Coelho on 22/01/2021.
//
//
#include "TPZSBFemMultiphysicsElGroup.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzcondensedcompel.h"
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
    if (fElGroup.size() == 0 || condensedmatid.size() != 4)
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
    Mesh()->ComputeNodElCon();

    bool keepmatrix = false;
    auto cond = new TPZCondensedCompEl(fCondensedEls, keepmatrix);
    Mesh()->ExpandSolution();
}

void TPZSBFemMultiphysicsElGroup::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    InitializeElementMatrix(ek, ef);

    TPZElementMatrix E0, E1, E2;
    ComputeMatrices(E0, E1, E2);

}

void TPZSBFemMultiphysicsElGroup::ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2)
{
    TPZElementMatrix eklocal, eflocal;
    for (auto cel : fElGroup)
    {
#ifdef PZDEBUG
        if (!cel)
        {
            DebugStop();
        }
#endif
        auto sbfemvol = dynamic_cast<TPZSBFemVolumeHdiv * >(cel);
#ifdef PZDEBUG
        if (!sbfemvol)
        {
            DebugStop();
        }
#endif
        TPZMaterial * material = sbfemvol->Material();
#ifdef PZDEBUG
        if(!material){
            PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
            eklocal.Reset();
            eflocal.Reset();
            return;
        }
#endif
        
        InitializeElementMatrix(eklocal, eflocal);
        
        // TPZManVector<TPZMaterialData,6> datavec;

        // auto ElementVec = sbfem->ElementVec();
        // const int64_t nref = ElementVec.size();

        // datavec.resize(nref);
        // InitMaterialData(datavec);
        
        // TPZManVector<TPZTransform<> > trvec;
        // AffineTransform(trvec);
        
        // int dim = Dimension();
        // TPZAutoPointer<TPZIntPoints> intrule;
        
        // TPZManVector<REAL,4> intpointtemp(TGeometry::Dimension,0.);
        // REAL weight = 0.;
        
        // TPZManVector<int,4> ordervec;
        // //ordervec.resize(nref);
        // for (int64_t iref=0;  iref<nref; iref++)
        // {
        //     TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        //     int svec;
        //     if(msp)
        //     {
        //         ordervec.Resize(ordervec.size()+1);
        //         svec = ordervec.size();
        //     }
        //     else
        //     {
        //         continue;
        //     }
        //     datavec[iref].p = msp->MaxOrder();
        //     ordervec[svec-1] = datavec[iref].p;
        // }
        // int order = material->IntegrationRuleOrder(ordervec);
        
        // TPZGeoEl *ref = this->Reference();
        // intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, order);
        
        // TPZManVector<int,4> intorder(dim,order);
        // intrule->SetOrder(intorder);
        // int intrulepoints = intrule->NPoints();
        // if(intrulepoints > 1000) {
        //     DebugStop();
        // }
        
        // TPZFMatrix<REAL> jac, axe, jacInv;
        // REAL detJac;
        // for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
        // {
        //     intrule->Point(int_ind,intpointtemp,weight);
        //     ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
        //     weight *= fabs(detJac);
        //     for (int i = 0; i < fElementVec.size(); i++) {
        //         TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[i].Element());
        //         if (!msp) {
        //             continue;
        //         }
        //         datavec[i].intLocPtIndex = int_ind;
        //     }
            
        //     this->ComputeRequiredData(intpointtemp,trvec,datavec);
            
        //     material->Contribute(datavec,weight,ek.fMat,ef.fMat);
        // }//loop over integration points
        
        // CleanupMaterialData(datavec);

    }
    

}