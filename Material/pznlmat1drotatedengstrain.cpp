/**
 * @file
 * @brief Contains implementations of the TPZNLMat1dRotatedEngStrain methods.
 */

#include "pznlmat1drotatedengstrain.h"

TPZNLMat1dRotatedEngStrain::TPZNLMat1dRotatedEngStrain(int id) : TPZNLMat1d(id)
{}


TPZNLMat1dRotatedEngStrain::~TPZNLMat1dRotatedEngStrain()
{}



REAL TPZNLMat1dRotatedEngStrain::Eps( TPZVec<REAL> &sol,
									 TPZFMatrix<REAL> &axes,
									 TPZFMatrix<REAL> &dphi)
{
	double theta = atan(axes (0,1) / axes (0,0) );
	double l = 2./ dphi(0,0);
	double x21 = l * cos(theta);
	double z21 = l * sin(theta);
	double u1,u2,w1,w2;
	u1 = sol [0];
	u2 = sol [2];
	w1 = sol [1];
	w2 = sol [3];
	
	double alpha0 = l/2.;
	double eps = ((1./4.)/alpha0)*((-x21*u1 + x21*u2 -z21*w1 +z21*w2));
	eps += (0.5/alpha0)*(u1*(u1-u2) + u2*(u2-u1) + w1*(w1-w2) + w2*(w2-w1));
	return eps;
}
