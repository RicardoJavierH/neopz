#ifndef TPZCOHESIVEBCH
#define TPZCOHESIVEBCH
/**
 * @file
 * @brief Contains the TPZCohesiveBC class which implements a cohesive boundary condition
 */


#include <iostream>
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzmatwithmem.h"

const int TPZCohesiveBCID = 400;

/**
 * @ingroup material
 * @brief This class implements a cohesive bc, which has the stress dependent of the displacement.
 * Remembering that the std::pair memory holds the pair (DeltaT,SigmaT) at each time step.
 * @author Nathan Shauer in 01/04/2014
 */
class  TPZCohesiveBC : public TPZMatWithMem<TPZFMatrix<REAL> >
{
private:

	/// SigmaT (Maximum Traction) of the cohesive equation
	REAL fSigmaT;
	/// DeltaC (Critical Displacement) of the cohesive equation
	REAL fDeltaC;
	/// DeltaT (diplacement for maximum traction). The point (DeltaT,SigmaT)
	REAL fDeltaT;
	
public:
	
	/** @brief Creates a material object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
	TPZCohesiveBC(int id);
	
	/** @brief Default constructor */
	TPZCohesiveBC();
	
	/** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
	TPZCohesiveBC(const TPZCohesiveBC &mat);
	
	/** @brief Default destructor */
	virtual ~TPZCohesiveBC();
	
	/** @brief Sets the Cohesive Data */
	void SetCohesiveData(const REAL &SigmaT, const REAL &DeltaC, const REAL &DeltaT);
	
	/** @brief Calculates Sigma for determined solution */
	void CalculateSigma(TPZVec<TPZMaterialData> &datavec, REAL &sigma) const; 

	/** @brief Updates the cohesive curve acording to the calculated w of the time step */
	void UpdateCohesiveCurve(TPZVec<TPZMaterialData> &datavec);
	
	/** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
	{
		datavec[0].SetAllRequirements(false);
		datavec[0].fNeedsSol = true;
	}
	
	/** @brief Returns the name of the material */
	virtual std::string Name() { return "TPZCohesiveBC"; }
	
	/** @brief Returns the integrable dimension of the material */
	virtual int Dimension() {
		return 1;
	}
		
	/** @brief Prints out the data associated with the material */
	virtual void Print(std::ostream &out = std::cout);
	
public:

	
	/**
	 * @brief It computes a contribution to the residual vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ef [out] is the residual vector
	 * @since April 16, 2007
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
	 * @param datavec [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	
	/** @brief Unique identifier for serialization purposes */
	virtual int ClassId() const
	{
		return TPZCohesiveBCID;
	}
	
};

#endif