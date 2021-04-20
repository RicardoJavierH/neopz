/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixTBBFlow methods.
 */

#include "pzstrmatrixflowtbb.h"
#include "TPZStructMatrix.h"
#include "TPZStructMatrixTBBFlowUtils.h"
#include "TPZGuiInterface.h"

#ifdef PZ_LOG
#include "pzlog.h"
static TPZLogger loggerel("pz.strmatrix.element");
#endif

#ifndef USING_TBB
#define NO_TBB \
    PZError<<"The class TPZStructMatrixTBBFlow depends on the TBB library.\n";\
    PZError<<"Please reconfigure the NeoPZ library using:\n";\
    PZError<<"USING_TBB=ON"<<std::endl;\
    DebugStop();
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif

template<class TVar>
TPZStructMatrixTBBFlow<TVar>::TPZStructMatrixTBBFlow() : TPZStrMatParInterface(), fFlowGraph(nullptr){
}

template<class TVar>
TPZStructMatrixTBBFlow<TVar>::TPZStructMatrixTBBFlow(const TPZStructMatrixTBBFlow &copy){
#ifndef USING_TBB
    NO_TBB
#else
        if (copy.fCompMesh) {
        fCompMesh = copy.fCompMesh;
    }
    fMaterialIds = copy.fMaterialIds;
    fNumThreads = copy.fNumThreads;    
    fFlowGraph = new TPZFlowGraph(*copy.fFlowGraph);
#endif
}

template<class TVar>
TPZStructMatrixTBBFlow<TVar>::~TPZStructMatrixTBBFlow(){}

#include "run_stats_table.h"

static RunStatsTable stat_ass_graph_tbb("-ass_graph_tbb", "Run statistics table for the graph creation, coloring and tbb::flow::graph TPZStructMatrixTBBFlow.");


static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

template<class TVar>
void TPZStructMatrixTBBFlow<TVar>::Assemble(TPZBaseMatrix & stiffness, TPZBaseMatrix & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    const auto &equationFilter =
        (dynamic_cast<TPZStructMatrix*>(this))->EquationFilter();
    ass_stiff.start();
    if (equationFilter.IsActive()) {
        int64_t neqcondense = equationFilter.NActiveEquations();
#ifdef PZDEBUG
        if (stiffness.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<STATE> rhsloc(neqcondense,rhs.Cols(),0.);
        this->MultiThread_Assemble(stiffness,rhsloc,guiInterface);
        equationFilter.Scatter(rhsloc, rhs);
    }
    else
    {
        this->MultiThread_Assemble(stiffness,rhs,guiInterface);
        
    }
    ass_stiff.stop();
}

template<class TVar>
void TPZStructMatrixTBBFlow<TVar>::Assemble(TPZBaseMatrix & rhs_base,TPZAutoPointer<TPZGuiInterface> guiInterface){
    const auto &equationFilter =
        (dynamic_cast<TPZStructMatrix*>(this))->EquationFilter();
    if(!dynamic_cast<TPZFMatrix<STATE>*>(&rhs_base)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" Incompatible types. Aborting...\n";
        DebugStop();
    }
    auto &rhs = dynamic_cast<TPZFMatrix<STATE> &>(rhs_base);
    ass_rhs.start();
    if(equationFilter.IsActive())
    {
        int64_t neqcondense = equationFilter.NActiveEquations();
        int64_t neqexpand = equationFilter.NEqExpand();
        //TODONORM
        if(rhs.Rows() != neqexpand || Norm(rhs) != 0.)
        {
            DebugStop();
        }
        TPZFMatrix<STATE> rhsloc(neqcondense,1,0.);
        this->MultiThread_Assemble(rhsloc,guiInterface);
        equationFilter.Scatter(rhsloc,rhs);
    }
    else
    {
        this->MultiThread_Assemble(rhs,guiInterface);
    }
    ass_rhs.stop();
}

template<class TVar>
TPZBaseMatrix * TPZStructMatrixTBBFlow<TVar>::CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
#ifndef USING_TBB
    NO_TBB
#else
    this->SetNumThreads(0);
    this->fFlowGraph = new TPZFlowGraph(this);
#endif
    auto *myself =
        dynamic_cast<TPZStructMatrix*>(this);
    const auto &equationFilter = myself->EquationFilter();
    TPZBaseMatrix *stiff = myself->Create();
    
    int64_t cols = MAX(1, rhs.Cols());
    rhs.Redim(equationFilter.NEqExpand(),cols);
    
    Assemble(*stiff,rhs,guiInterface);
    
#ifdef PZ_LOG2
    if(loggerel.isDebugEnabled())
    {
        std::stringstream sout;
        stiff->Print("Stiffness matrix",sout);
        rhs.Print("Right hand side", sout);
        LOGPZ_DEBUG(loggerel,sout.str())
    }
#endif
    return stiff;
    
}

template<class TVar>
void TPZStructMatrixTBBFlow<TVar>::MultiThread_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
#ifdef USING_TBB
    this->fFlowGraph->ExecuteGraph(&rhs, &mat);
#else
    std::cout << "To use the tbb flow graph assemble please compile the NeoPZ with USING_TBB." << std::endl;
    DebugStop();
#endif
}

template<class TVar>
void TPZStructMatrixTBBFlow<TVar>::MultiThread_Assemble(TPZBaseMatrix & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
#ifdef USING_TBB
    this->fFlowGraph->ExecuteGraph(&rhs);
#else
    std::cout << "To use the tbb flow graph assemble please compile the NeoPZ with USING_TBB." << std::endl;
    DebugStop();
#endif
}

template<class TVar>
int TPZStructMatrixTBBFlow<TVar>::ClassId() const{
    return Hash("TPZStructMatrixTBBFlow") ^ TPZStrMatParInterface::ClassId() << 1;
}

template<class TVar>
void TPZStructMatrixTBBFlow<TVar>::Read(TPZStream& buf, void* context) {
    TPZStrMatParInterface::Read(buf,context);
}

template<class TVar>
void TPZStructMatrixTBBFlow<TVar>::Write(TPZStream& buf, int withclassid) const {
    TPZStrMatParInterface::Write(buf,withclassid);
}


template class TPZStructMatrixTBBFlow<STATE>;
template class TPZRestoreClass<TPZStructMatrixTBBFlow<STATE>>;
