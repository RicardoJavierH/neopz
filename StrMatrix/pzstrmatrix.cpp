/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixOR methods.
 */
#include <thread>
#include <vector>
#include "pzstrmatrix.h"

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include "pzelementgroup.h"
#include "pzanalysis.h"
#include "pzsfulmat.h"

#include "pzgnode.h"
#include "TPZTimer.h"

#include "pzcheckmesh.h"
#include "pzcheckconsistency.h"
#include "TPZMaterial.h"
#include "run_stats_table.h"

using namespace std;

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.TPZStructMatrixOR");
static TPZLogger loggerel("pz.strmatrix.element");
static TPZLogger loggerel2("pz.strmatrix.elementinterface");
static TPZLogger loggerelmat("pz.strmatrix.elementmat");
static TPZLogger loggerCheck("pz.strmatrix.checkconsistency");
static TPZLogger loggerGlobStiff("pz.strmatrix.globalstiffness");
#endif

#ifdef PZ_LOG
static TPZLogger loggerCS("calcstiffTime");
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif


TPZStructMatrixOR::TPZStructMatrixOR(TPZCompMesh *mesh) : TPZStructMatrixBase(mesh) {
    
}

TPZStructMatrixOR::TPZStructMatrixOR(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrixBase(cmesh) {
    
}

TPZMatrix<STATE> *TPZStructMatrixOR::Create() {
    cout << "TPZStructMatrixOR::Create should never be called\n";
    return 0;
}

TPZStructMatrixOR *TPZStructMatrixOR::Clone() {
    cout << "TPZStructMatrixOR::Clone should never be called\n";
    DebugStop();
    return 0;
}

TPZMatrix<STATE> * TPZStructMatrixOR::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                          unsigned numthreads_assemble, unsigned numthreads_decompose) {
  SetNumThreads(numthreads_assemble);
  return CreateAssemble(rhs, guiInterface);
}

static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixOR::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    ass_stiff.start();
    if (fEquationFilter.IsActive()) {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
#ifdef PZDEBUG
        if (stiffness.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<STATE> rhsloc(neqcondense, rhs.Cols(), 0.);
        if (this->fNumThreads) {
            this->MultiThread_Assemble(stiffness, rhsloc, guiInterface);
        } else {
            this->Serial_Assemble(stiffness, rhsloc, guiInterface);
        }

        fEquationFilter.Scatter(rhsloc, rhs);
    } else {
        if (this->fNumThreads) {
            this->MultiThread_Assemble(stiffness, rhs, guiInterface);
        } else {
            this->Serial_Assemble(stiffness, rhs, guiInterface);
        }
    }
    ass_stiff.stop();
}

void TPZStructMatrixOR::Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    ass_rhs.start();
    if (fEquationFilter.IsActive()) {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
        int64_t neqexpand = fEquationFilter.NEqExpand();
        if (rhs.Rows() != neqexpand || Norm(rhs) != 0.) {
            DebugStop();
        }
        TPZFMatrix<STATE> rhsloc(neqcondense, 1, 0.);
        if (this->fNumThreads) {
            this->MultiThread_Assemble(rhsloc, guiInterface);
        } else {
            this->Serial_Assemble(rhsloc, guiInterface);
        }
        fEquationFilter.Scatter(rhsloc, rhs);
    } else {
        if (this->fNumThreads) {
            this->MultiThread_Assemble(rhs, guiInterface);
        } else {
            this->Serial_Assemble(rhs, guiInterface);
        }
    }
    ass_rhs.stop();
}
double calcstiffTime;
void TPZStructMatrixOR::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
#ifdef PZDEBUG
    TExceptionManager activateExceptions;
#endif
    if (!fMesh) {
        LOGPZ_ERROR(logger, "Serial_Assemble called without mesh")
        DebugStop();
    }
#ifdef PZ_LOG
    if (loggerelmat.isDebugEnabled()) {
        if (dynamic_cast<TPZSubCompMesh *> (fMesh)) {
            std::stringstream sout;
            sout << "AllEig = {};";
            LOGPZ_DEBUG(loggerelmat, sout.str())
        }
    }
#endif

#ifdef PZDEBUG
    if (rhs.Rows() != fEquationFilter.NActiveEquations()) {
        DebugStop();
    }
#endif

    int64_t iel;
    int64_t nelem = fMesh->NElements();
    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
#ifdef PZ_LOG
    bool globalresult = true;
    bool writereadresult = true;
#endif
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();

    int64_t count = 0;
    for (iel = 0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if (!el) continue;
        int matid = 0;
        TPZGeoEl *gel = el->Reference();
        if (gel) {
            matid = gel->MaterialId();
        }
        int matidsize = fMaterialIds.size();
        if(matidsize){
            if(!el->NeedsComputing(fMaterialIds)) continue;
        }

        count++;
        if (!(count % 1000)) {
            std::cout << '*';
            std::cout.flush();
        }
        if (!(count % 20000)) {
            std::cout << "\n";
        }
        
        calcstiff.start();
        ek.Reset();
        ef.Reset();
/*
#ifdef PZ_LOG
        TPZTimer timer;
        if(loggerCS.isDebugEnabled()){
            timer.start();}
#endif
*/
        el->CalcStiff(ek,ef);
/*
#ifdef PZ_LOG
        if(loggerCS.isDebugEnabled()){
        timer.stop();
        calcstiffTime += timer.seconds();
        }
#endif
*/
        if (guiInterface) if (guiInterface->AmIKilled()) {
                return;
            }

#ifdef PZ_LOG
        if (loggerelmat.isDebugEnabled()) {
            if (dynamic_cast<TPZSubCompMesh *> (fMesh)) {
                std::stringstream objname;
                objname << "Element" << iel;
                std::string name = objname.str();
                objname << " = ";
                std::stringstream sout;
                ek.fMat.Print(objname.str().c_str(), sout, EMathematicaInput);
                sout << "AppendTo[AllEig,Eigenvalues[" << name << "]];";

                LOGPZ_DEBUG(loggerelmat, sout.str())
                        /*		  if(iel == 133)
                         {
                         std::stringstream sout2;
                         el->Reference()->Print(sout2);
                         el->Print(sout2);
                         LOGPZ_DEBUG(logger,sout2.str())
                         }
                         */
            }
        }
#endif

#ifdef CHECKCONSISTENCY
        //extern TPZCheckConsistency stiffconsist("ElementStiff");
        stiffconsist.SetOverWrite(true);
        bool result;
        result = stiffconsist.CheckObject(ek.fMat);
        if (!result) {
            globalresult = false;
            std::stringstream sout;
            sout << "element " << iel << " computed differently";
            LOGPZ_ERROR(loggerCheck, sout.str())
        }

#endif

        calcstiff.stop();
        /*
        assemble.start();

        if (!ek.HasDependency()) {
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            //			TPZSFMatrix<STATE> test(stiffness);
            //			TPZFMatrix<STATE> test2(stiffness.Rows(),stiffness.Cols(),0.);
            //			stiffness.Print("before assembly",std::cout,EMathematicaInput);
            stiffness.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
#ifdef PZDEBUG
            REAL rhsnorm = Norm(ef.fMat);
            REAL eknorm = Norm(ek.fMat);
            if (std::isnan(rhsnorm) || std::isnan(eknorm)) {
                std::cout << "element " << iel << " has norm " << rhsnorm << std::endl;
                el->Print();
                ek.fMat.Print("ek",std::cout);
                ef.fMat.Print("ef",std::cout);
            }
#endif
            rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
            //			stiffness.Print("stiffness after assembly STK = ",std::cout,EMathematicaInput);
            //			rhs.Print("rhs after assembly Rhs = ",std::cout,EMathematicaInput);
            //			test2.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //			test -= stiffness;
            //			test.Print("matriz de rigidez diference",std::cout);
            //			test2.Print("matriz de rigidez interface",std::cout);
#ifdef PZ_LOG
            if (loggerel.isDebugEnabled()) {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                if (gel) {
                    TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                    gel->CenterPoint(gel->NSides() - 1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for computational element index " << el->Index() << " material id " << gel->MaterialId() << std::endl;
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                } else {
                    sout << "Stiffness for computational element without associated geometric element index " << el->Index() << "\n";
                }
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel, sout.str())
            }
#endif
        } else {
            // the element has dependent nodes
            ek.ApplyConstraints();
            ef.ApplyConstraints();
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            stiffness.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);

#ifdef PZ_LOG
            if (loggerel.isDebugEnabled()) {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                //                el->Print();
                //                int nc = el->NConnects();
                //                for (int ic=0; ic<nc; ic++) {
                //                    std::cout << "Index " << el->ConnectIndex(ic) << " ";
                //                    el->Connect(ic).Print(*fMesh);
                //                    fMesh->ConnectVec()[ic].Print(*fMesh);
                //                }
                if (gel) {
                    TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                    gel->CenterPoint(gel->NSides() - 1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                } else {
                    sout << "Stiffness for computational element index " << iel << std::endl;
                }
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel, sout.str())
            }
#endif
        }
        // tototototo
        //        GK.Multiply(Sol, GKSol);
        //        GKSol -= GF;
        //        GKSol.Transpose();
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "GKSol" << iel << " = ";
        //            GKSol.Print(str.str().c_str(),sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //        stiffness.Multiply(Sol, GKSol);
        //        GKSol -= rhs;
        //        GKSol.Transpose();
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "StiffSol" << iel << " = ";
        //            GKSol.Print(str.str().c_str(),sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "GK" << iel << " = ";
        //            GK.Print(str.str().c_str(),sout,EMathematicaInput);
        //            std::stringstream str2;
        //            str2 << "ST" << iel << " = ";
        //            stiffness.Print(str2.str().c_str(),sout,EMathematicaInput);
        //            sout << "GK-ST\n";
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //        
        //        stiffness.Zero();
        //        rhs.Zero();
        assemble.stop();
    */}//fim for iel
    if (count > 1000) std::cout << std::endl;

#ifdef PZ_LOG
    if (loggerCheck.isDebugEnabled()) {
        std::stringstream sout;
        sout << "The comparaison results are : consistency check " << globalresult << " write read check " << writereadresult;
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
    if (loggerGlobStiff.isDebugEnabled())
    {
        std::stringstream sout;
        stiffness.Print("GK = ",sout,EMathematicaInput);
        rhs.Print("GR = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerel,sout.str())
    }

#endif

}

void TPZStructMatrixOR::Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {

    int64_t iel;
    int64_t nelem = fMesh->NElements();

    TPZTimer calcresidual("Computing the residual vector");
    TPZTimer assemble("Assembling the residual vector");

    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();

    for (iel = 0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if (!el || el->IsInterface()) continue;

        int matid = 0;
        TPZGeoEl *gel = el->Reference();
        if (gel) {
            matid = gel->MaterialId();
        }
        int matidsize = fMaterialIds.size();
        if(matidsize){
            TPZMaterial * mat = el->Material();
            TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
            if (!mat)
            {
                if (!submesh) {
                    continue;
                }
                else if(submesh->NeedsComputing(fMaterialIds) == false) continue;
            }
            else
            {
                if (this->ShouldCompute(matid) == false) continue;
            }
        }
                
        TPZElementMatrix ef(fMesh, TPZElementMatrix::EF);

        calcresidual.start();

        el->CalcResidual(ef);

        calcresidual.stop();
#ifdef PZ_LOG
        if(loggerel.isDebugEnabled())
        {
            std::stringstream sout;
            TPZGeoEl *gel = el->Reference();
            if (gel) {
                TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                gel->CenterPoint(gel->NSides() - 1, center);
                gel->X(center, xcenter);
                sout << "Residual for computational element index " << el->Index() << " material id " << gel->MaterialId() << std::endl;
                sout << "Residual for geometric element " << gel->Index() << " center " << xcenter << std::endl;
            } else {
                sout << "Residual for computational element without associated geometric element index " << el->Index() << "\n";
            }
            LOGPZ_DEBUG(loggerel, sout.str())
        }
#endif
        assemble.start();

        if (!ef.HasDependency()) {
            ef.ComputeDestinationIndices();
            fEquationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
        } else {
            // the element has dependent nodes
            ef.ApplyConstraints();
            ef.ComputeDestinationIndices();
            fEquationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat, ef.fSourceIndex, ef.fDestinationIndex);
        }
#ifdef PZ_LOG
        if(loggerel.isDebugEnabled())
        {
            std::stringstream sout;
            ef.Print(sout);
            LOGPZ_DEBUG(loggerel, sout.str())
        }
#endif

        assemble.stop();

    }//fim for iel
#ifdef PZ_LOG
    {
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << calcresidual.processName() << " " << calcresidual << std::endl;
            sout << assemble.processName() << " " << assemble;
            LOGPZ_DEBUG(logger, sout.str().c_str());
        }
    }
#endif
    //std::cout << std::endl;
}

void TPZStructMatrixOR::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    ThreadData threaddata(this,mat,rhs,fMaterialIds,guiInterface);
    const int numthreads = this->fNumThreads;
    std::vector<std::thread> allthreads;
    int itr;
    if (guiInterface) {
        if (guiInterface->AmIKilled()) {
            return;
        }
    }
    for (itr = 0; itr < numthreads; itr++) {
      allthreads.push_back(std::thread(ThreadData::ThreadWork, &threaddata));
    }

    ThreadData::ThreadAssembly(&threaddata);

    for (itr = 0; itr < numthreads; itr++) {
      allthreads[itr].join();
    }

#ifdef PZ_LOG2
    if (loggerCheck.isDebugEnabled()) {
        std::stringstream sout;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        mat.Print("Matriz de Rigidez: ", sout, EMathematicaInput);
        rhs.Print("Right Handside", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
#endif
}

void TPZStructMatrixOR::MultiThread_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    ThreadData threaddata(this, rhs, fMaterialIds, guiInterface);
    const int numthreads = this->fNumThreads;
    std::vector<std::thread> allthreads;
    int itr;
    if (guiInterface) {
        if (guiInterface->AmIKilled()) {
            return;
        }
    }

    for (itr = 0; itr < numthreads; itr++) {
      allthreads.push_back(std::thread(ThreadData::ThreadWork, &threaddata));
    }

    ThreadData::ThreadAssembly(&threaddata);

    for (itr = 0; itr < numthreads; itr++) {
      allthreads[itr].join();
    }
}

TPZStructMatrixOR::ThreadData::ThreadData(TPZStructMatrixOR *strmat, TPZMatrix<STATE> &mat,
        TPZFMatrix<STATE> &rhs,
        std::set<int> &MaterialIds,
        TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat), fGuiInterface(guiInterface), fGlobMatrix(&mat), fGlobRhs(&rhs), fNextElement(0) {

}

int TPZStructMatrixOR::ClassId() const{
    return Hash("TPZStructMatrixOR") ^ TPZStructMatrixBase::ClassId() << 1;
}

void TPZStructMatrixOR::Read(TPZStream& buf, void* context) {
    TPZStructMatrixBase::Read(buf,context);
}

void TPZStructMatrixOR::Write(TPZStream& buf, int withclassid) const {
    TPZStructMatrixBase::Write(buf,withclassid);
}


TPZStructMatrixOR::ThreadData::ThreadData(TPZStructMatrixOR *strmat,
        TPZFMatrix<STATE> &rhs,
        std::set<int> &MaterialIds,
        TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat), fGuiInterface(guiInterface), fGlobMatrix(0), fGlobRhs(&rhs), fNextElement(0) {
}

TPZStructMatrixOR::ThreadData::~ThreadData() {
}

//#define DRY_RUN

void *TPZStructMatrixOR::ThreadData::ThreadWork(void *datavoid) {
#ifdef PZDEBUG
    TExceptionManager activateExceptions;
#endif
    ThreadData *data = (ThreadData *) datavoid;
    // compute the next element (this method is threadsafe)
    int64_t iel = data->NextElement();
    TPZCompMesh *cmesh = data->fStruct->fMesh;
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    int64_t nel = cmesh->NElements();
    while (iel < nel) {

        TPZAutoPointer<TPZElementMatrix> ek;
        TPZAutoPointer<TPZElementMatrix> ef = new TPZElementMatrix(cmesh, TPZElementMatrix::EF);
        if (data->fGlobMatrix) {
            ek = new TPZElementMatrix(cmesh, TPZElementMatrix::EK);
        } else {
            ek = ef;
        }

        TPZCompEl *el = cmesh->ElementVec()[iel];
        TPZElementMatrix *ekp = ek.operator->();
        TPZElementMatrix *efp = ef.operator->();
        TPZElementMatrix &ekr = *ekp;
        TPZElementMatrix &efr = *efp;

#ifndef DRY_RUN
        if (data->fGlobMatrix) {
            el->CalcStiff(ekr, efr);
        } else {
            el->CalcResidual(efr);
        }
#else
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (el);
            if (!intel) {
                DebugStop();
            }
            if (data->fGlobMatrix) {
                intel->InitializeElementMatrix(ekr, efr);
            } else {
                intel->InitializeElementMatrix(efr);
            }
        }
#endif

        if (guiInterface) if (guiInterface->AmIKilled()) {
                break;
            }

        if (!el->HasDependency()) {
            ek->ComputeDestinationIndices();

            if (data->fStruct->HasRange()) {
                data->fStruct->FilterEquations(ek->fSourceIndex, ek->fDestinationIndex);
            }
#ifdef PZ_LOG
            if (loggerel.isDebugEnabled()) {
                std::stringstream sout;
                sout << "Element index " << iel << std::endl;
                ek->fMat.Print("Element stiffness matrix", sout);
                ef->fMat.Print("Element right hand side", sout);
                LOGPZ_DEBUG(loggerel, sout.str())
            }
#endif
        } else {
            // the element has dependent nodes
            if (data->fGlobMatrix) {
                ek->ApplyConstraints();
            }
            ef->ApplyConstraints();
            ek->ComputeDestinationIndices();
            if (data->fStruct->HasRange()) {
                data->fStruct->FilterEquations(ek->fSourceIndex, ek->fDestinationIndex);
            }
#ifdef PZ_LOG
            if (loggerel2.isDebugEnabled() && el->Reference() && el->Reference()->MaterialId() == 1 && el->IsInterface()) {
                std::stringstream sout;
                el->Reference()->Print(sout);
                el->Print(sout);
                ek->Print(sout);
                //			ef->Print(sout);
                LOGPZ_DEBUG(loggerel2, sout.str())
            }
#endif
#ifdef PZ_LOG
            if (loggerel.isDebugEnabled()) {
                std::stringstream sout;
                sout << "Element index " << iel << std::endl;
                ek->fConstrMat.Print("Element stiffness matrix", sout);
                ef->fConstrMat.Print("Element right hand side", sout);
                LOGPZ_DEBUG(loggerel, sout.str())
            }
#endif
        }


        // put the elementmatrices on the stack to be assembled (threadsafe)
        data->ComputedElementMatrix(iel, ek, ef);
        // compute the next element (this method is threadsafe)
        iel = data->NextElement();
    }

   {
      std::scoped_lock lock(data->fMutexAccessElement);
      data->fAssembly.Post();
    }

    return 0;
}

// The function which will compute the assembly

void *TPZStructMatrixOR::ThreadData::ThreadAssembly(void *threaddata) {
    ThreadData *data = (ThreadData *) threaddata;
    TPZCompMesh *cmesh = data->fStruct->fMesh;
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    int64_t nel = cmesh->NElements();
    data->fMutexAccessElement.lock();
    int64_t nextel = data->fNextElement;
    int numprocessed = data->fProcessed.size();
    while (nextel < nel || numprocessed) {
        if (guiInterface) if (guiInterface->AmIKilled()) {
            break;//mutex will still be unlocked at the end of the function
            }
        std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > >::iterator itavail;
        std::set<int>::iterator itprocess;
        bool keeplooking = false;
        if (data->fSubmitted.size() && data->fProcessed.size()) {
            itavail = data->fSubmitted.begin();
            itprocess = data->fProcessed.begin();
            if (itavail->first == *itprocess) {
                // make sure we come back to look for one more element
                keeplooking = true;
                // Get a hold of the data
#ifdef PZ_LOG
                int iel = *itprocess;
#endif
                data->fProcessed.erase(itprocess);
                TPZAutoPointer<TPZElementMatrix> ek = itavail->second.first;
                TPZAutoPointer<TPZElementMatrix> ef = itavail->second.second;
                data->fSubmitted.erase(itavail);
#ifdef PZ_LOG
                if (logger.isDebugEnabled()) {
                    std::stringstream sout;
                    sout << "Assembling element " << iel;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
#ifdef CHECKCONSISTENCY
                //static TPZCheckConsistency stiffconsist("ElementStiff");
                stiffconsist.SetOverWrite(true);
                bool result;
                result = stiffconsist.CheckObject(ek->fMat);
                if (!result) {
                    globalresult = false;
                    std::stringstream sout;
                    sout << "element " << iel << " computed differently";
                    LOGPZ_ERROR(loggerCheck, sout.str())
                }
#endif

                // Release the mutex
                data->fMutexAccessElement.unlock();
#ifndef DRY_RUN
                auto globMatrix =
                    dynamic_cast<TPZMatrix<STATE> *>(data->fGlobMatrix);
                auto globRhs =
                    dynamic_cast<TPZFMatrix<STATE> *>(data->fGlobRhs);
                // Assemble the matrix
                if (!ek->HasDependency()) {
                    if (globMatrix) {
                        globMatrix->AddKel(ek->fMat, ek->fSourceIndex, ek->fDestinationIndex);
                    }
                    globRhs->AddFel(ef->fMat, ek->fSourceIndex, ek->fDestinationIndex);
                } else {
                    if (globMatrix) {
                        globMatrix->AddKel(ek->fConstrMat, ek->fSourceIndex, ek->fDestinationIndex);
                    }
                    globRhs->AddFel(ef->fConstrMat, ek->fSourceIndex, ek->fDestinationIndex);
                }
#endif
                // acquire the mutex
                data->fMutexAccessElement.lock();
            }
        }
        if (!keeplooking) {
          data->fMutexAccessElement.unlock();
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                LOGPZ_DEBUG(logger, "Going to sleep within assembly")
            }
#endif
            // wait for a signal
            data->fAssembly.Wait();
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                LOGPZ_DEBUG(logger, "Waking up for assembly")
            }
#endif
            data->fMutexAccessElement.lock();
        }
        nextel = data->fNextElement;
        numprocessed = data->fProcessed.size();

    }
    //	std::cout << std::endl;
#ifdef PZ_LOG
    if (loggerCheck.isDebugEnabled()) {
        std::stringstream sout;
        sout << "nextel = " << nextel << " numprocessed = " << numprocessed << " submitted " << data->fSubmitted.size() << std::endl;
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
#endif
#ifdef PZ_LOG
    if (loggerCheck.isDebugEnabled()) {
        std::stringstream sout;
        if (data->fGlobMatrix) {
            data->fGlobMatrix->Print("Matriz de Rigidez: ", sout, EMathematicaInput);
        }
        data->fGlobRhs->Print("Right hand side", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
#endif
    data->fMutexAccessElement.unlock();
    return 0;
}

int64_t TPZStructMatrixOR::ThreadData::NextElement() {
    fMutexAccessElement.lock();
    int64_t iel;
    int64_t nextel = fNextElement;
    TPZCompMesh *cmesh = fStruct->fMesh;
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    int64_t nel = elementvec.NElements();
    for (iel = fNextElement; iel < nel; iel++) {
        TPZCompEl *el = elementvec[iel];
        if (!el) continue;
        if (fStruct->fMaterialIds.size() == 0) break;
        if (el->NeedsComputing(fStruct->fMaterialIds)) break;
    }
    fNextElement = iel + 1;
    nextel = iel;
    if (iel < nel) fProcessed.insert(iel); //AQUIBORIN pelo que percebi, aqui tem que acontecer antes do unlock no caso sem Critical Section
    fMutexAccessElement.unlock();
#ifdef PZ_LOG
    {
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << " returning " << nextel << " fNextElement " << fNextElement;
            LOGPZ_DEBUG(logger, sout.str())
        }
    }
#endif
    return nextel;
}


// put the computed element matrices in the map

void TPZStructMatrixOR::ThreadData::ComputedElementMatrix(int64_t iel, TPZAutoPointer<TPZElementMatrix> &ek, TPZAutoPointer<TPZElementMatrix> &ef) {
    std::scoped_lock lock(fMutexAccessElement);
    std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > el(ek, ef);
    fSubmitted[iel] = el;
    fAssembly.Post();
}

template class TPZRestoreClass<TPZStructMatrixOR>;
