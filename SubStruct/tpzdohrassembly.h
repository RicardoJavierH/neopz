/**
 * @file
 * @brief Contains the TPZDohrAssembly class which implements assembling using Dohrmann algorithm.
 */
/*
 *  tpzdohrassembly.h
 *  SubStruct
 *
 *  Created by Philippe Devloo on 04/03/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 *
 */

#ifndef TPZDOHRASSEMBLYH
#define TPZDOHRASSEMBLYH

#include "pzvec.h"

class TPZFMatrix;

/**
 * @ingroup substructure
 * @brief Assembling using Dohrmann algorithm. \ref substructure "Sub structure"
 */
class TPZDohrAssembly 
{
public:
	/**
	 * @brief For each substructure the equation numbering of the substructures
	 * 
	 * The order of the equations follows the ordering of the connects
	 */
	TPZVec< TPZVec< int > > fFineEqs;
	
	/** @brief For each substructure the equation numbering of the coarse equations */
	TPZVec< TPZVec< int > > fCoarseEqs;
	
	/** @brief Sum the values in the local matrix into the global matrix */
	void Assemble(int isub, const TPZFMatrix &local, TPZFMatrix &global);
	
	/** @brief Extract the values from the global matrix into the local matrix */
	void Extract(int isub, const TPZFMatrix &global, TPZFMatrix &local);
	
	/** @brief Sum the values in the local matrix into the global matrix */
	void AssembleCoarse(int isub, const TPZFMatrix &local, TPZFMatrix &global);
	
	/** @brief Extract the values from the global matrix into the local matrix */
	void ExtractCoarse(int isub, const TPZFMatrix &global, TPZFMatrix &local);
};

#endif
