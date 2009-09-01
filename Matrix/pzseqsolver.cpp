/* Generated by Together */

#include "pzseqsolver.h"
using namespace std;

TPZSequenceSolver::TPZSequenceSolver(TPZMatrix *refmat) : TPZMatrixSolver(refmat), fSolvers() {
}
TPZSequenceSolver::TPZSequenceSolver(const TPZSequenceSolver & copy)
    : TPZMatrixSolver(copy) {
    int nums = copy.fSolvers.NElements();
    int s;
    for(s=0; s<nums; s++) AppendSolver(*copy.fSolvers[s]);
}
class TPZSolver;
class TPZMatrixSolver;

TPZSolver * TPZSequenceSolver::Clone() const {
    return new TPZSequenceSolver(*this);
}
void TPZSequenceSolver::AppendSolver(TPZMatrixSolver & solve){
    fSolvers.Push((TPZMatrixSolver *) solve.Clone());
}
void TPZSequenceSolver::ResetSolver() {
    int nums = fSolvers.NElements();
    int s;
    for(s=0; s<nums; s++) delete fSolvers.Pop();
}
void TPZSequenceSolver::Solve(const TPZFMatrix &F, TPZFMatrix &result, TPZFMatrix *residual){
  if(!Matrix()) {
    cout << "TPZSequenceSolver::Solve called without a matrix pointer\n";
    exit(-1);
  }
  TPZAutoPointer<TPZMatrix> mat = Matrix();
  if(result.Rows() != mat->Rows() || result.Cols() != F.Cols()) {
    result.Redim(mat->Rows(),F.Cols());
  }

  fScratch = F;
  TPZFMatrix delu(result.Rows(),result.Cols(),0.);
  TPZFMatrix resloc(F.Rows(),F.Cols(),0.);
  result.Zero();
  if(fScratch.Rows() != result.Rows() || fScratch.Cols() != result.Cols()) {
	  fScratch.Redim(result.Rows(),result.Cols());
  }
    int nums = fSolvers.NElements();
    int s;
    for(s=0; s<nums; s++) {
        fSolvers[s]->Solve(fScratch,delu,&resloc);
        result += delu;
        mat->Residual(result,F,fScratch);
    }
    if(residual) *residual = fScratch;
}

  
  /**
  This method will reset the matrix associated with the solver
  This is useful when the matrix needs to be recomputed in a non linear problem
  */
void TPZSequenceSolver::ResetMatrix()
{
    int nums = fSolvers.NElements();
    int s;
    for(s=0; s<nums; s++) {
        fSolvers[s]->ResetMatrix();
    }
    TPZMatrixSolver::ResetMatrix();
}
  
  /**
  Updates the values of the preconditioner based on the values of the matrix
  */
void TPZSequenceSolver::UpdateFrom(TPZAutoPointer<TPZMatrix> matrix)
{
    int nums = fSolvers.NElements();
    int s;
    for(s=0; s<nums; s++) {
        fSolvers[s]->UpdateFrom(matrix);
    }
    TPZMatrixSolver::UpdateFrom(matrix);
}

void TPZSequenceSolver::Write(TPZStream &buf, int withclassid)
{
  TPZMatrixSolver::Write(buf, withclassid);
  int StackSz = fSolvers.NElements();
  buf.Write(&StackSz, 1);
  int i = 0;
  for(i = 0; i < StackSz; i++)
  {
    fSolvers[i]->Write(buf, 1);
  }

}
void TPZSequenceSolver::Read(TPZStream &buf, void *context)
{
  TPZMatrixSolver::Read(buf, context);
  int StackSz = 0;
  buf.Read(&StackSz, 1);
  fSolvers.Resize(StackSz);
  int i = 0;
  for(i = 0; i< StackSz; i++)
  {
    fSolvers[i] = dynamic_cast<TPZMatrixSolver *>(TPZSaveable::Restore(buf, context));
  }
}

template class TPZRestoreClass< TPZSequenceSolver, TPZSQUENCESOLVER_ID>;
