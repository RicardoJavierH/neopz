/* Generated by Together */
class TPZFMatrix;
class TPZBndCond;

#ifndef TPZLINEARCONVECTION_H
#define TPZLINEARCONVECTION_H
#include "pzmaterial.h"

/**
 * This class implements a linear scalar convection equation with modified 
 * SUPG difusion
 */
class TPZLinearConvection : public TPZMaterial {
public:  
    TPZLinearConvection(TPZLinearConvection & copy);

    TPZLinearConvection(int id,TPZVec<REAL> &conv) ;


    /**returns the integrable dimension of the material*/
    virtual int Dimension() ;

    /** returns the number of state variables associated with the material*/
    virtual int NStateVariables()  ;

    /** return the number of components which form the flux function*/
    virtual int NFluxes() {return 2;}

    /**
       * Contribute adds the contribution to the stiffness matrix
       **/
    virtual void Contribute(TPZMaterialData &data, REAL weight,
  			  TPZFMatrix &ek,TPZFMatrix &ef);


    virtual void ContributeBC(TPZMaterialData &data,REAL weight,
    			    TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

    /**
       * Contribute adds the contribution to the stiffness matrix
       **/
    virtual void Contribute(TPZMaterialData &data, REAL weight,
			  TPZFMatrix &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}


	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
					TPZFMatrix &ef,TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
    }

	/** print out the data associated with the material*/
    virtual void Print(std::ostream &out = std::cout);

    /**returns the variable index associated with the name*/
    virtual int VariableIndex(const std::string &name);

    /** returns the number of variables associated with the variable indexed by var.
   *       var is obtained by calling VariableIndex*/
    virtual int NSolutionVariables(int var);

protected:
	/**returns the solution associated with the var index based on the finite element approximation*/
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
public:
	/**returns the solution associated with the var index based on the finite element approximation*/
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<REAL> &Solout)
	{
        Solution(data.sol,data.dsol,data.axes,var,Solout);
    }

    /**compute the value of the flux function to be used by ZZ error estimator*/
    virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux) {}

    /**To create another material of the same type*/
    virtual TPZAutoPointer<TPZMaterial> NewMaterial();

    /**Read data of the material from a istream (file data)*/
    virtual void SetData(std::istream &data);

private:    
    REAL fConvect[2];
};
#endif //TPZLINEARCONVECTION_H
