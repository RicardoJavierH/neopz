//
//  TPZSBFemElementGroup.hpp
//  PZ
//
//  Created by Karolinne Coelho on 22/01/21.
//
//

#pragma once

#include <stdio.h>
#include "TPZSBFemVolumeHdiv.h"
#include "TPZSBFemElementGroup.h"
#include "pzcondensedcompel.h"

using namespace std;

class TPZSBFemMultiphysicsElGroup : public TPZSBFemElementGroup
{
    
private:
    
    TPZElementGroup * fCondensedEls;

    TPZCondensedCompEl * fCondEl;
    
    /// Matrix of eigenvectors which compose the stiffness matrix
    TPZFMatrix<std::complex<double>> fPhi;
    
    /// Inverse of the eigenvector matrix (transfers eigenvector coeficients to side shape coeficients)
    TPZFMatrix<std::complex<double>> fPhiInverse;
    
    /// Vector of eigenvalues of the SBFem analyis
    TPZManVector<std::complex<double>> fEigenvalues;
    
    /// Multiplying coefficients of each eigenvector
    TPZFMatrix<std::complex<double>> fCoef;
    
    TPZManVector<int64_t> fLocalindices;

public:
    
    TPZSBFemMultiphysicsElGroup() : TPZSBFemElementGroup()
    {
        
    }
    
    /// constructor
    TPZSBFemMultiphysicsElGroup(TPZCompMesh &mesh, int64_t &index) : TPZSBFemElementGroup(mesh,index)
    {
        
    }

    void AddElement(TPZCompEl *cel);

    virtual void Print(std::ostream &out) const;

    void GroupandCondense(set<int> & matidscondensed);


    void ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2, TPZElementMatrix &M0);
    
    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);


    void ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2);


    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override
    {
        for (auto sbfemvol : fElGroup)
        {
            sbfemvol->BuildCornerConnectList(connectindexes);
        }
    }

    void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef) const;

    void LoadSolution();
    
    TPZFMatrix<std::complex<double>> &PhiInverse()
    {
        return fPhiInverse;
    }

    void AdjustConnectivities();

    void SetLocalIndices(int64_t index);
    
    TPZManVector<double> EigenvaluesReal(TPZManVector<complex<double> > & eigenvalues)
    {
        int64_t nel = eigenvalues.NElements();
        TPZManVector<double> eig(nel);
        for(int64_t el=0; el<nel; el++)
        {
            eig[el] = eigenvalues[el].real();
        }
        return eig;
    }

    TPZFMatrix<double> EigenVectorsReal(TPZFMatrix<std::complex<double>> &eigvec)
    {
        auto nrows = eigvec.Rows();
        auto ncols = eigvec.Cols();
        TPZFMatrix<double> eig(nrows,ncols);

        for (auto ir = 0; ir < nrows; ir++)
        {
            for (auto ic = 0; ic < ncols; ic++)
            {
                eig(ir,ic) = eigvec(ir,ic).real();
            }
        }
        
        return eig;
    }
};