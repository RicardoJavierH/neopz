#ifndef PZRANDOMFIELD_H
#define PZRANDOMFIELD_H

#include "pzfunction.h"
#include <iostream>
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzerror.h"
#include "tpzverysparsematrix.h"
#include "pzsfulmat.h"


#include <math.h>
#include <complex>
#include <string>
#include <random>
#include <chrono>


template<class TVar>
public class TPZRandomField : public TPZFunction<TVar>
{

protected:
    
    void (*fFunc)(const TPZVec<REAL> &x, TPZVec<TVar> &f);
    void (*fFunc2)(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf);
    void (*fFunc3)(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf);
    void (*fFunc4)(const TPZVec<REAL> &f, int id);
    
    int fPorder;
    int fnSquareElements;
    int fmatsize;
    int fstochasticInclined;
    REAL fdirection;
    REAL finclination;
    REAL frw; // wellbore radius
    REAL frext; // external radius
    REAL fH; // cylinder total height
    REAL fh; // elements height (square elements for now)
    TPZGeoMesh* fgmesh; //geometric mesh
    TPZFMatrix<TVar> fK; // correlation matrix
    TPZFMatrix<TVar> fRand_U; //random distribution
    TPZFMatrix<TVar> fU; // random correlated* distribution
    TPZFMatrix<TVar> fM; //Decomposed matrix from Mathematica

    TPZFMatrix<TVar> fU_E; // random correlated* distribution of E
    TPZFMatrix<TVar> fU_nu; // random correlated* distribution of nu
   
    // normal ==1;
    // lognormal ==2;
    int fE_dist;
    int fnu_dist;
    
    // exponential ==1;
    // spherical ==2;
    int fE_funct;
    int fnu_funct;
    
    TPZFMatrix<TVar> fKE; // correlation matrix of E
    TPZFMatrix<TVar> fKnu; // correlation matrix of nu
    
    
public:
    
    /** @brief Class constructor */
    TPZRandomField(TPZGeoMesh* geometricMesh, int numSquareElems, int stochasticInclined, REAL direction,
                   REAL inclination, REAL rw, REAL rext, const TPZFMatrix<TVar> &M) : TPZFunction<TVar>(), fPorder(-1)
    {
        fFunc  = 0;
        fFunc2 = 0;
        fFunc3 = 0;
        fFunc4 = 0;
        
        //This should be given by the user
        fgmesh = geometricMesh;
        
        //Any way to get this information from the geometric mesh?
        fnSquareElements = numSquareElems;  // number of Square Elements
        fstochasticInclined = stochasticInclined;
        fdirection = direction;
        finclination = inclination;
        frw = rw;
        frext = rext;
        
        // this should not be here, delete after SVD implementation
        fM = M;
        
        //** Calculation of stochastic field should start here!! *****/////
        
        //This should be called by the object
        SetFieldDistribution(fM);
        
        //Set inclined parameters
        if(stochasticInclined==true) {
            SetInclinedField(stochasticInclined,direction, inclination);
        }
        
        //Chnage geometry in case it is an inclined well
        SetFieldGeometry();
        
        // this should be given by the user
        fE_dist = 1;
        fnu_dist = 1;
        fE_funct = 1;
        fnu_funct = 1;
        //SetYoungField(E_distribution, E_function)
        //SetPoissonField(nudistribution, nu_function)
            
        //Get fK
        if(fE_funct==fnu_funct){
            fKE = EvaluateCorrelation(fE_funct);
            fKnu = fKE;
            
            // Create function to decompose fK using SVD decomposition
            // fK = U sqt(S) V'
            // create glogal parameter for U and sqrt(S)
        }
        
        else {
        //Evaluate for E
        fKE = EvaluateCorrelation(fE_funct);
        //Evaluate for nu
        fKnu = EvaluateCorrelation(fnu_funct);
            
            // Create function to decompose fK using SVD decomposition
            // fKE = U sqt(S) V'
            // fKnu = U sqt(S) V'
            // create glogal parameter for U and sqrt(S)
        }
        
        
        //Write a function to get fU
        fU_E = GetCorrelatedVector(fE_dist);
        fU_nu = GetCorrelatedVector(fnu_dist);
        
        //GetStochasticField(fE, fnu); ***** Should use this one to get the random correlated and scaled field, in TPZMaterial use fE and fnu
        
        
        //Overwrtie fU just to test - This should be deleted
        // Multiplying decomposed Matrix M (U*Sqrt(S)) and random normal vector fRand_U
        fU = fM * fU_E; // Get correlated random distribution
        // In this function fM should be replaced by the left singular vetor U and the square root of the diagonal matrix S, then multiply by fRand_U
        
        
        //******* End
    }
    
    /** @brief Copy constructor */
    TPZRandomField &operator=(const TPZRandomField &cp)
    {
        fFunc = cp.fFunc;
        fFunc2 = cp.fFunc2;
        fFunc3 = cp.fFunc3;
        fFunc4 = cp.fFunc4;
        fPorder = cp.fPorder;
        
        fgmesh = cp.fgmesh;
        fnSquareElements = cp.fnSquareElements;  // number of Square Elements
        fstochasticInclined = cp.fstochasticInclined;
        fdirection = cp.fdirection;
        finclination = cp.finclination;
        
        // this should not be here, try to get from fgmesh
        frw = cp.frw;
        frext = cp.frext;
        fM = cp.fM;
        
        return *this;
    }
    
    /** @brief Class destructor */
    virtual ~TPZRandomField()
    {
        
    }
    
    /** @brief Set distribution type and function of E */
    void SetYoungField(int distribution, int function){
    }
    
    /** @brief Set distribution type and function of nu */
    void SetPoissonField(int distribution, int function){
    }
    
    /** @brief Set distribution type: normal or lognormal */
    void SetFieldDistribution(const TPZFMatrix<TVar> &M)
    {
    }
    
    /** @brief Set inclined parameters if wellbore is inclined */
    virtual void SetInclinedField(int stochasticInclined,REAL direction, REAL inclination)
    {
    }
    
    /** @brief Set geometric parameters  */
    virtual void SetFieldGeometry()
    {
    }
    
    /** @brief Set inclined geometry */
    virtual void InclinedFieldGeometry()
    {
    }
    
    /** @brief Calculates Correlation either vertical or inclined */
     TPZFMatrix<TVar>  EvaluateCorrelation(int function)
    {
    }
    
    
    
    /** @brief Calculates Correlation either vertical or inclined */
    TPZFMatrix<REAL> GetCorrelatedVector(int distribution)
    {
    }
    
    /** @brief Get vector of distribution type */
    TPZFMatrix<TVar> GetDistribution(int matrixSize, int distribution)
    {
    }
    
    /** @brief Calculates Correlation either vertical or inclined */
    virtual void GetStochasticField( TPZFMatrix<TVar> f_E, TPZFMatrix<TVar> f_nu)
    {
    }
    
    /** @brief Calculates correlation matrix for vertical well */
    TPZFMatrix<REAL> calcCorrelationMatrix(int function)
    {
    }
    
    
    /** @brief Calculates correlation matrix for inclined well */
    TPZFMatrix<REAL> calcCorrelationMatrixInclined(int function) {
    }
    
    
    /** @brief Print correlation matrix to be decomposed at Mathematica - This should be removed after SVD decomposition!!! */
    virtual void PrintCorrelation()
    {
    }
    
    
    /*////*/ /*////*/ /*////*/ /*////*/ /*////*/ /*////*/ /*////*/ /*////*/ /*////*/ /*////*//*////*/ /*////*/ /*////*/ /*////*/ /*////*/
    /*////*/
    /*////*/
    
    
    
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, TPZVec<TVar> &val))
    {
    }
    
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, TPZVec<TVar> &val, TPZFMatrix<TVar> &gradf))
    {
    }
    
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &val, TPZFMatrix<TVar> &gradf))
    {
    }
    
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &f, int id))
    {
    }
    
    TPZRandomField(const TPZRandomField &cp) : fFunc(cp.fFunc), fFunc2(cp.fFunc2), fFunc3(cp.fFunc3), fFunc4(cp.fFunc4), fPorder(cp.fPorder)
    {
        
    }
    
    
    /**
     * @brief Performs function computation
     * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
     * @param f function values
     * @param df function derivatives
     */
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df)
    {
    }
    
    /**
     * @brief Performs time dependent function computation
     * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
     * @param ftime  time to evaluate
     * @param f function values
     * @param gradf function derivatives
     */
    virtual void Execute(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf)
    {
    }
    
    /**
     * @brief Execute method receiving axes. It is used in shape functions
     * @note NOT IMPLEMENTED
     */
    virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<TVar> &f, TPZFMatrix<TVar> &df)
    {
    }
    
    /**
     * @brief Simpler version of Execute method which does not compute function derivatives
     * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
     * @param f function values
     */
    // This method gives random values of Young Modulus from a specific range
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f) {
        
    }
    
    // Call this method in the PZMatElasticity 2D (Gaussian Field)
    virtual void Execute(const TPZVec<TVar> &f, int id) {
    }
    
    // Call this method in the PZMatElasticity for Element Failure Area
    virtual void ExecuteArea(const TPZVec<TVar> &f, int id)
    {
    }
    
    // Calc Correlation Matrix
    TPZVec<TVar> calcStochasticField()
    {
    }
    
    
    /** @brief Returns number of functions. */
    virtual int NFunctions()
    {
    }
    
    void SetPolynomialOrder(int porder)
    {
    }
    
    /** @brief Polynomial order of this function. */
    /** In case of non-polynomial function it can be a reasonable approximation order. */
    virtual int PolynomialOrder()
    {
    }
    
    /** @brief Unique identifier for serialization purposes */
    virtual int ClassId() const
    {
    }
    
    /** @brief Saves the element data to a stream */
    virtual void Write(TPZStream &buf, int withclassid)
    {
    }
    
    /** @brief Reads the element data from a stream */
    virtual void Read(TPZStream &buf, void *context)
    {
    }
};

#endif
