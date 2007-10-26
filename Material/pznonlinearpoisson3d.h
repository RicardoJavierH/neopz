// -*- c++ -*-

//$Id: pznonlinearpoisson3d.h,v 1.5 2007-10-26 13:13:41 tiago Exp $

#ifndef MATNLPOISSON3DH
#define MATNLPOISSON3DH

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

class TPZNonLinearPoisson3d : public TPZMatPoisson3dReferred {

 protected:

  /**
   * Definitions of stabilization terms
   */
  enum EStabilizationType {ENoStabilization = 0, ESUPG = 1, EGradient = 2};

  /** Stabilization term definition
   */
  int fStabilizationType;

 public:

  TPZNonLinearPoisson3d(int nummat, int dim);
  
  TPZNonLinearPoisson3d(const TPZNonLinearPoisson3d &cp);
  
  virtual ~TPZNonLinearPoisson3d();
  
  bool IsReferred(){ return this->fIsReferred;}
  
  void SetReferred(bool Is){ this->fIsReferred = Is; }
  
  /** Define SUPG stabilization term.
   */
  void SetSUPGStab(REAL sd = 1.0);
  
  /** Define gradient stabilization term.
   */
  void SetGradientStab(REAL sd = 1.0);
  
  /** Define no stabilization term.
   */
  void SetNoStabilizationTerm();

  virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef);
               
  virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix &ek,
                              TPZFMatrix &ef,
                              TPZBndCond &bc);
                            
  virtual void ContributeInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ek,
                                     TPZFMatrix &ef);

  virtual void ContributeBCInterface(TPZMaterialData &data,
                                       REAL weight,
                                       TPZFMatrix &ek,
                                       TPZFMatrix &ef,
                                       TPZBndCond &bc);

  protected:
    bool fIsReferred;

};

#endif
