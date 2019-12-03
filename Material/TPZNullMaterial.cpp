/**
 * @file
 * @brief Contains implementations of the TPZMaterial methods.
 */

#include "TPZNullMaterial.h"
#include "TPZNullMaterialTranslator.h"
#include "pzmaterialdata.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzbndcond.h"
#include "pzreal.h"
#include "pzadmchunk.h"
#include "tpzintpoints.h"

#include "pzlog.h"

#ifdef LOG4CXX
#ifdef PZDEBUG
#define DEBUG2
#endif
static LoggerPtr logger(Logger::getLogger("pz.material"));
#endif



TPZNullMaterial::TPZNullMaterial() : TPZRegisterClassId(&TPZNullMaterial::ClassId),
TPZMaterial() {
    fDim    =-1;
    fNState = 1;
}

TPZNullMaterial::~TPZNullMaterial() {
    
}

void TPZNullMaterial::Print(std::ostream & out) {
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZMaterial::Print(out);
    
}

int TPZNullMaterial::VariableIndex(const std::string &name) {
    return TPZMaterial::VariableIndex(name );
}

int TPZNullMaterial::NSolutionVariables(int index) {
    return TPZMaterial::NSolutionVariables(index);
}

void TPZNullMaterial::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	this->Solution(data.sol[0], data.dsol[0], data.axes, var, Solout);
}

void TPZNullMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	DebugStop();
}

void TPZNullMaterial::Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout)
{
}

void TPZNullMaterial::Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl *left, TPZCompEl *right)
{
    
}

void TPZNullMaterial::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
						   TPZVec<STATE> &Solout){
    if(var == 0) Solout = Sol;

    else if (var==1){
        STATE val = 0.;
        for(int i=0; i<fDim; i++){
            val += DSol(i,i);
        }
        Solout[0] = val;
    }

    else
    {
        TPZMaterial::Solution(Sol, DSol,axes,var,Solout);
    }
}


void TPZNullMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
}

void TPZNullMaterial::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
}

void TPZNullMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
}




int TPZNullMaterial::ClassId() const{
    return Hash("TPZNullMaterial") ^ TPZMaterial::ClassId() << 1;
}

void TPZNullMaterial::Write(TPZStream &buf, int withclassid) const
{
	TPZMaterial::Write(buf,withclassid);
    if (fDim < 1 || fDim >3) {
        DebugStop();
    }
    buf.Write(&fDim);
    buf.Write(&fNState);
}

void TPZNullMaterial::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf,context);
    buf.Read(&fDim);
#ifdef PZDEBUG
    if (fDim < 1 || fDim >3) {
        DebugStop();
    }
#endif
    buf.Read(&fNState);
}

void TPZNullMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
}

void TPZNullMaterial::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){

}
void TPZNullMaterial::ErrorsHdiv(TPZMaterialData &data,TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values){
    
}

template class TPZRestoreClassWithTranslator<TPZNullMaterial,TPZNullMaterialTranslator>;
