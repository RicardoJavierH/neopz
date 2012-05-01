/* Generated by Together */

#ifndef TPZADAPTMESH_H
#define TPZADAPTMESH_H

#include "pzcmesh.h"
#include "pzcclonemesh.h"
#include "pzvec.h"

class TPZInterpolatedElement;
class TPZTransfer;
class TPZTransform;
template<class T, class V>
class TPZAvlMap;
class TPZOneDRef;

/** @brief Interface to generate adapted meshes */
class TPZAdaptMesh {
 public:
	 void RemoveCloneBC(TPZCompMesh *mesh);

  /** @brief Simple constructor */
  TPZAdaptMesh();    
  
  /** @brief Simple destructor */
  ~TPZAdaptMesh();
  
  /** @brief Defines the computational reference mesh */
  void SetCompMesh(TPZCompMesh * mesh);

  /** @brief Defines the maximum p order of an element */
  void SetMaxP(int maxp);
  
  /**
   * @brief Public interface to get the optmally refined mesh 
   * @param error: returns the estimated error
   * @param truerror: returns the true error if analitical solution is provided
   * @param ervec: estimated element error for original mesh element vector
   * @param f: analitical solution
   * @param truervec: real element error at each orginal mesh element
   * @param effect: error estimator effectivity
   * @param use_trueerror: evaluates the error throgh the analitical solution provided by f
   */
  TPZCompMesh * GetAdaptedMesh(REAL &error,
			       REAL &truerror,
			       TPZVec<REAL> &ervec, 
			       void (*f)(const TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix<REAL> &deriv),
			       TPZVec<REAL> &truervec, 
			       TPZVec<REAL> &effect,
			       int use_trueerror = 0);


static void DeleteElements(TPZCompMesh *mesh);

  REAL UseTrueError(TPZInterpolatedElement *coarse, void (*f)(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<REAL> &deriv));

 protected:
  
  
  /** @brief Retrieves the geometric reference elements to create the patches */
  void GetReferenceElements();
  
  /**
   * @brief Builds the patch of all reference elements. \n
   * The patches are stored into patch vectors
   */   
  void BuildReferencePatch();
  
  /** @brief Fill the vector of clone meshes */
  void CreateClones();
  
  /** @brief Sorts the elements by the error vector vec, returning permutation vector */
  void Sort(TPZVec<REAL> &vec, TPZVec<int> &perm);
  
  /** @brief Sort */
  void HeapSort(TPZVec<REAL> &sol, TPZVec<int> &perm);
  
  /**
   * @brief Sorts the errvec returning the ordering indexes in perm param.
   * @param errvec vector of errors to sort
   * @param perm ordering indexes after sorting
   * @param errpercent is the percentual of the error that must be considered in returning minimum error
   */
  REAL SortMinError (TPZVec<REAL> errvec, TPZVec<int> perm, REAL errpercent);

  /**
   * @brief Creates an adpted computational mesh based on original mesh and in a hp refinement pattern also
   * @param mesh original mesh
   * @param gelstack h refinement pattern given by a list of an adapted geometric elements
   * @param porders p refinement pattern for each element of gelstack
   */
  TPZCompMesh* CreateCompMesh (TPZCompMesh *mesh,TPZVec<TPZGeoEl *> &gelstack,TPZVec<int> &porders);

  /**
   * @brief Verifies if one clone, specified by its index, must be analysed \n
   * This method only be called when the true solution is available and the 
   * option usetrueerror in void  GetAdaptedMesh is seted to 1.
   * @param clindex index of the clone to be verified
   * @param minerror minimum error to the clone be analysed
   * @param ervec vector containing the treu error 
   */
  int HasTrueError(int clindex, REAL &minerror, TPZVec<REAL> &ervec);


 private:   
  
  static TPZInterpolatedElement * LargeElement(TPZInterpolatedElement *cint);
  
  /** @brief Computational reference mesh */
  TPZCompMesh *fReferenceCompMesh;
  
  /** @brief Geometric reference elements vector */
  TPZStack < TPZGeoEl * > fGeoRef;
  
  /** @brief Patches vector */
  TPZStack < TPZGeoEl * > fPatch;
  
  /** @brief Maps the start position of each patch into patches vector */
  TPZStack < int > fPatchIndex;
  
  /** @brief Element error vector */
  TPZStack < REAL > fElementError;
  
  /** @brief Clone meshes vector */
  TPZStack<TPZCompCloneMesh *> fCloneMeshes ;
  
  /** @brief Refined clone meshes */
  TPZStack <TPZCompMesh *> fFineCloneMeshes ;
  
  /** @brief Delete temporary clone meshes from memory */
  void CleanUp();

  /** @brief Maximum p order of an element */
  int fMaxP;
  
};

#endif //TPZADAPTMESH_H
