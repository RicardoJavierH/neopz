//
//  TRMBuildTransfers.cpp
//  PZ
//
//  Created by Omar on 10/27/15.
//
//


#include "TRMBuildTransfers.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


/** @brief Default constructor */
TRMBuildTransfers::TRMBuildTransfers(){
    
    fSimulationData = NULL;

}

/** @brief Default desconstructor */
TRMBuildTransfers::~TRMBuildTransfers(){

}

/** @brief Copy constructor $ */
TRMBuildTransfers::TRMBuildTransfers(const TRMBuildTransfers &copy)
{
    fSimulationData = copy.fSimulationData;
}

/** @brief Copy assignemnt operator $ */
TRMBuildTransfers & TRMBuildTransfers::operator=(const TRMBuildTransfers &other)
{
    if (this != & other) // prevent self-assignment
    {
        fSimulationData = other.fSimulationData;
    }
    return *this;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Transfers:: Iterative Coupling by Operator Splitting //////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Segregated Transfer methods (Gamma and Omega) :: Build methods Elliptic
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @brief bluid linear applications: u and grad_u to elliptic $ */
void TRMBuildTransfers::Build_elliptic_To_elliptic(TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!elliptic) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    elliptic->LoadReferences();    
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    
    fe_e_cindexes.resize(0);
    std::pair< long, long > chunk_geo_cel_indexes;
    
    // Step 1 :: Counting for valid elements (apply all the needed filters in this step)
    for (long i = 0; i < geometry->NElements(); i++) {
        
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        
        int mat_id = gel->MaterialId();
        if ((dim == 2 && (mat_id >= 8 && mat_id <= 11)) || (dim == 3 && (mat_id >= 8 && mat_id <= 13))) { // Filtering bc reservoir elements
            continue;
        }
        
        chunk_geo_cel_indexes.first = gel->Index();
        chunk_geo_cel_indexes.second = -1;
        fe_e_cindexes.Push(chunk_geo_cel_indexes);
    }
    
    long n_el = fe_e_cindexes.size();
    fu_dof_scatter.resize(n_el);
    fe_e_intp_indexes.resize(n_el);

    
    // Block size structue including (Omega and Gamma)
    std::pair<long, TPZVec<long> > chunk_intp_indexes;
    TPZVec< std::pair<long, long> > blocks_dimensions_phi(n_el);
    TPZVec< std::pair<long, long> > blocks_dimensions_grad_phi(n_el);
    
    // Step 2 :: filling linking vectors
    for (long iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        
#ifdef PZDEBUG
        if(!mf_e_cel)
        {
            DebugStop();
        }
#endif
        
        int gel_dim = gel->Dimension();
        
        // Geometry and cel link
        fe_e_cindexes[iel].second = e_cel->Index();
        
        // Getting local integration index
        TPZManVector<long> e_int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        mf_e_cel->GetMemoryIndices(e_int_point_indexes);
        
        chunk_intp_indexes.first = gel->Index();
        chunk_intp_indexes.second  = e_int_point_indexes;
        fe_e_intp_indexes.Push(chunk_intp_indexes);
        
        this->ElementDofIndexes(mf_e_cel, dof_indexes);
        fu_dof_scatter[iel] = dof_indexes;
        
        
        blocks_dimensions_phi[iel].first = e_int_point_indexes.size()*dim;
        blocks_dimensions_phi[iel].second = dof_indexes.size();
        
        blocks_dimensions_grad_phi[iel].first = e_int_point_indexes.size()*dim*gel_dim;
        blocks_dimensions_grad_phi[iel].second = dof_indexes.size();
    }
    
    // Step 3 :: Initialize the matrix
    fu_To_elliptic.Initialize(blocks_dimensions_phi);
    fgrad_u_To_elliptic.Initialize(blocks_dimensions_grad_phi);
    
    
    TPZManVector<long> e_int_point_indexes(0,0);
    TPZManVector<long> dof_indexes(0,0);
    
    // Step 4 :: Filling the matrix
    std::pair<long, long> block_dim;
    for (long iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        
        
#ifdef PZDEBUG
        if(!mf_e_cel)
        {
            DebugStop();
        }
#endif
        
        
        TPZInterpolationSpace * e_intel = dynamic_cast<TPZInterpolationSpace * >(mf_e_cel->Element(0));
        
        // Getting local integration index
        mf_e_cel->GetMemoryIndices(e_int_point_indexes);
        dof_indexes = fu_dof_scatter[iel];
        
        block_dim.first = e_int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_geomechanic = mf_e_cel->GetIntegrationRule();
        int np_cel = int_points_geomechanic.NPoints();
        
#ifdef PZDEBUG
        if (e_int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi, block_grad_phi;
        
        block_phi.Resize(block_dim.first*dim,block_dim.second);
        block_grad_phi.Resize(block_dim.first*dim*gel_dim,block_dim.second);
        
        // for derivatives in real space
        int nshape = e_intel->NShapeF();
        TPZFNMatrix<220> phi(nshape,1);
        TPZFNMatrix<660> dphi(gel_dim,nshape),dphix_axes(gel_dim,nshape);
        TPZFMatrix<double> dphidx;
        TPZFNMatrix<9,STATE> jacobian(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> jacinv(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> axes;
        REAL detjac;
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(gel_dim,0.0);
            STATE w;
            int_points_geomechanic.Point(ip, qsi, w);
            
            // Get the phi and dphix for H1 elasticity
            e_intel->Shape(qsi, phi, dphi);
            gel->Jacobian( qsi, jacobian, axes, detjac , jacinv);
            
            switch(gel_dim) {
                case 0:
                    break;
                case 1:
                    dphix_axes = dphi;
                    dphix_axes *= (1./detjac);
                    break;
                case 2:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
                    }
                    break;
                case 3:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
                        dphix_axes(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
                    }
                    break;
                default:
                    std::stringstream sout;
                    sout << "pzintel.c please implement the " << gel_dim << "d Jacobian and inverse\n";
                    LOGPZ_ERROR(logger,sout.str());
            }
            
            TPZAxesTools<STATE>::Axes2XYZ(dphix_axes, dphidx, axes);
            
#ifdef PZDEBUG
            if(block_dim.second != phi.Rows() * dim){
                DebugStop();
            }
#endif
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < dim; id++) {
                    block_phi(ip*dim+id,jp*dim+id) = phi(jp,0);
                }
            }
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < dim; id++) {
                    for (int jd = 0; jd < gel_dim; jd++) {
                        if(gel_dim == dim){
                            block_grad_phi(ip*dim*gel_dim + id*gel_dim + jd,jp*dim+id) = dphidx(jd,jp);
                        }
                        else{
                            block_grad_phi(ip*dim*gel_dim+id*gel_dim + jd,jp*dim+id) = 0.0;
                        }
                    }
                }
            }
            
        }
        
        fu_To_elliptic.SetBlock(iel, block_phi);
        fgrad_u_To_elliptic.SetBlock(iel, block_grad_phi);

    }
    
//    fu_To_elliptic.Print(" u_to_e ");
//    fgrad_u_To_elliptic.Print(" grad_u_to_e ");
    
    return;
    
}

void TRMBuildTransfers::space_To_elliptic(TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!elliptic || fu_To_elliptic.Rows() == 0 || fgrad_u_To_elliptic.Rows() == 0) {
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = elliptic->Dimension();
    int n_el = fe_e_cindexes.size();
    
    
//    long first_point_phi = 0;
//    long first_point_dphi = 0;
    std::pair<long, long> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    TPZFMatrix<double> block_phi;
    TPZFMatrix<double> block_grad_phi;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(fe_e_cindexes[iel].second);
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        
#ifdef PZDEBUG
        if(!mf_e_cel)
        {
            DebugStop();
        }
#endif
        
//        first_point_phi += b_size_phi.first;
//        first_point_dphi += b_size_dphi.first;
        b_size_phi = fu_To_elliptic.GetSizeofBlock(iel);
        b_size_dphi = fgrad_u_To_elliptic.GetSizeofBlock(iel);
        fu_To_elliptic.GetBlock(iel, block_phi);
        fgrad_u_To_elliptic.GetBlock(iel, block_grad_phi);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        TPZMaterial * material = elliptic->FindMaterial(matd_id);
        
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            int gel_dim = gel->Dimension();
            int nshapes = b_size_phi.second/dim;
            TPZFNMatrix<3,STATE> phi_u(nshapes,1,0.0);
            TPZFNMatrix<9,STATE> grad_phi_u(dim,nshapes);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int is = 0; is < nshapes; is++) {
                    phi_u(is,0) = block_phi(ip*dim,is*dim);
                    
                }
                
                for (int is = 0; is < nshapes; is++) {
                    for (int id = 0 ; id < dim; id++) {
                        grad_phi_u(id,is) = block_grad_phi(ip*dim*gel_dim + id,is*dim);
                    }
                    
                }
                
                associated_material->GetMemory()[ipos].Set_phi_u(phi_u);
                associated_material->GetMemory()[ipos].Set_grad_phi_u(grad_phi_u);
            }
            
            
        }
        else{
            
            TPZMatWithMem<TRMMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            int gel_dim = gel->Dimension();
            int nshapes = b_size_phi.second/dim;
            TPZFNMatrix<3,STATE> phi_u(nshapes,1,0.0);
            TPZFNMatrix<9,STATE> grad_phi_u(dim,nshapes);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int is = 0; is < nshapes; is++) {
                    phi_u(is,0) = block_phi(ip*dim,is*dim);
                    
                }
                
                for (int is = 0; is < nshapes; is++) {
                    for (int id = 0 ; id < dim; id++) {
                        grad_phi_u(id,is) = block_grad_phi(ip*dim*gel_dim + id,is*dim);
                    }
                    
                }
                
                associated_material->GetMemory()[ipos].Set_phi_u(phi_u);
                associated_material->GetMemory()[ipos].Set_grad_phi_u(grad_phi_u);
            }
            
        }

        
    }
    
}

void TRMBuildTransfers::elliptic_To_elliptic(TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!elliptic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_u(fu_To_elliptic.Cols(),1,0.0);
    int n = fe_e_cindexes.size();
    long pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fu_dof_scatter[i].size(); iequ++) {
            Scatter_u(pos,0) = elliptic->Solution()(fu_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> u_at_elliptic,grad_u_at_elliptic;
    fu_To_elliptic.Multiply(Scatter_u,u_at_elliptic);
    fgrad_u_To_elliptic.Multiply(Scatter_u, grad_u_at_elliptic);
    
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = elliptic->Dimension();
    int n_el = fe_e_cindexes.size();
    
    long iblock = 0;
    long first_point_phi = 0;
    long first_point_dphi = 0;
    std::pair<long, long> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(fe_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        
#ifdef PZDEBUG
        if(!mf_e_cel)
        {
            DebugStop();
        }
#endif
        
        first_point_phi += b_size_phi.first;
        first_point_dphi += b_size_dphi.first;
        b_size_phi = fu_To_elliptic.GetSizeofBlock(iblock);
        b_size_dphi = fgrad_u_To_elliptic.GetSizeofBlock(iblock);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMaterial * material = elliptic->FindMaterial(matd_id);
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            TPZFNMatrix<3,STATE> u(1,3,0.0);
            TPZFNMatrix<9,STATE> grad_u(3,3,0.0);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    u(0,id) = u_at_elliptic(first_point_phi + ip*dim + id,0);
                }
                
                for (int id = 0; id < dim ; id++) {
                    for (int jd = 0; jd < dim ; jd++) {
                        grad_u(jd,id)= grad_u_at_elliptic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                    }
                }
                
                if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                    associated_material->GetMemory()[ipos].Set_grad_u_0(grad_u);
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_u_n(u);
                    associated_material->GetMemory()[ipos].Set_grad_u_n(grad_u);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_u(u);
                    associated_material->GetMemory()[ipos].Set_grad_u(grad_u);
                }
                
            }
            
            
        }
        else{
            
            TPZMaterial * material = elliptic->FindMaterial(matd_id);
            TPZMatWithMem<TRMMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            TPZFNMatrix<3,STATE> u(1,3,0.0);
            TPZFNMatrix<9,STATE> grad_u(3,3,0.0);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    u(0,id) = u_at_elliptic(first_point_phi + ip*dim + id,0);
                }
                
                for (int id = 0; id < dim ; id++) {
                    for (int jd = 0; jd < dim ; jd++) {
                        grad_u(jd,id)= grad_u_at_elliptic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                    }
                }
                
                if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                    associated_material->GetMemory()[ipos].Set_grad_u_0(grad_u);
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_u_n(u);
                    associated_material->GetMemory()[ipos].Set_grad_u_n(grad_u);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_u(u);
                    associated_material->GetMemory()[ipos].Set_grad_u(grad_u);
                }
                
            }
            
        }
        
        iblock++;
        
    }
    
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Segregated Transfer methods (Gamma and Omega) :: Build methods Parabolic
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void TRMBuildTransfers::Build_parabolic_To_parabolic(TPZCompMesh * parabolic){
    
#ifdef PZDEBUG
    if (!parabolic) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    parabolic->LoadReferences();
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    
    fp_p_cindexes.resize(0);
    std::pair< long, long > chunk_geo_cel_indexes;
    
    // Step 1 :: Counting for valid elements (apply all the needed filters in this step)
    for (long i = 0; i < geometry->NElements(); i++) {
        
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        
        
        int mat_id = gel->MaterialId();
        if ( (dim == 2 && mat_id > 11) || (dim == 3 && mat_id > 13) ) { // Filtering bc reservoir elements
            continue;
        }
        
        chunk_geo_cel_indexes.first = gel->Index();
        chunk_geo_cel_indexes.second = -1;
        fp_p_cindexes.Push(chunk_geo_cel_indexes);

    }

    
    long n_el = fp_p_cindexes.size();
    fq_dof_scatter.resize(n_el);
    fp_dof_scatter.resize(n_el);
    fp_p_intp_indexes.resize(n_el);
    
    std::pair<long, TPZVec<long>  > chunk_intp_indexes;
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions_phi_q(n_el);
    TPZVec< std::pair<long, long> > blocks_dimensions_div_phi_q(n_el);
    TPZVec< std::pair<long, long> > blocks_dimensions_phi_p(n_el);
    
    int q_index = 0;
    int p_index = 1;
    
    int q_points = 0;
    int div_q_points = 0;
    int p_points = 0;
    
    for (long iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!p_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
#ifdef PZDEBUG
        if(!mf_p_cel)
        {
            DebugStop();
        }
#endif
        
        // Geometry and cel link
        fp_p_cindexes[iel].second = p_cel->Index();
        
        // Getting local integration index
        TPZManVector<long> p_int_point_indexes(0,0);
        
        TPZManVector<long> q_dof_indexes(0,0);
        TPZManVector<long> p_dof_indexes(0,0);
        
        int vec_dim = dim;
        mf_p_cel->GetMemoryIndices(p_int_point_indexes);
        q_points        = p_int_point_indexes.size();
        div_q_points    = p_int_point_indexes.size();
        p_points        = p_int_point_indexes.size();
        
        this->ElementDofIndexes(mf_p_cel, q_dof_indexes,q_index);

        if (gel->Dimension() == dim) {
            this->ElementDofIndexes(mf_p_cel, p_dof_indexes,p_index);

        }
        else{
            div_q_points = 0;
            p_points     = 0;
            vec_dim      = 1;
        }

        fq_dof_scatter[iel] = q_dof_indexes;
        fp_dof_scatter[iel] = p_dof_indexes;
        
        blocks_dimensions_phi_q[iel].first = q_points*vec_dim;
        blocks_dimensions_phi_q[iel].second = q_dof_indexes.size();
        
        blocks_dimensions_div_phi_q[iel].first = div_q_points;
        blocks_dimensions_div_phi_q[iel].second = q_dof_indexes.size();
        
        blocks_dimensions_phi_p[iel].first = p_points;
        blocks_dimensions_phi_p[iel].second = p_dof_indexes.size();
    
        
    }
    
    // Initialize the matrix
    fp_To_parabolic.Initialize(blocks_dimensions_phi_p);
    fq_To_parabolic.Initialize(blocks_dimensions_phi_q);
    fdiv_q_To_parabolic.Initialize(blocks_dimensions_div_phi_q);

    
    TPZManVector<long> p_int_point_indexes(0,0);
    
    std::pair<long, long> q_block_dim;
    std::pair<long, long> div_q_block_dim;
    std::pair<long, long> p_block_dim;
    
    // for velocity functions
    TPZMaterialData data;
    
    for (long iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fp_p_cindexes[iel].second);
        
#ifdef PZDEBUG
        if (!p_cel) {
            DebugStop();
        }
#endif
        
        
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
#ifdef PZDEBUG
        if(!mf_p_cel)
        {
            DebugStop();
        }
#endif
        
        
        // Getting local integration index
        mf_p_cel->GetMemoryIndices(p_int_point_indexes);

        q_block_dim     = fq_To_parabolic.GetSizeofBlock(iel);
        div_q_block_dim = fdiv_q_To_parabolic.GetSizeofBlock(iel);
        p_block_dim     = fp_To_parabolic.GetSizeofBlock(iel);
        
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points = mf_p_cel->GetIntegrationRule();
        int np_cel = int_points.NPoints();
        
#ifdef PZDEBUG
        if (p_int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi_q, block_div_phi_q, block_phi_p;
        
        block_phi_q.Resize(q_block_dim.first,q_block_dim.second);
        block_div_phi_q.Resize(div_q_block_dim.first, div_q_block_dim.second);
        block_phi_p.Resize(p_block_dim.first,p_block_dim.second);
        
        // Velocity functions
        TPZInterpolationSpace * q_intel = dynamic_cast<TPZInterpolationSpace * >(mf_p_cel->Element(0));
        if(q_intel)
        {
            
            // Computing over all integration points of the compuational element cel
            TPZFNMatrix<100,REAL> phi(q_intel->NShapeF(),1,0.0);
            int el_dim = gel->Dimension();
            TPZFNMatrix<300,REAL> dphidxi(el_dim,q_intel->NShapeF(),0.0);

            for (int ip = 0; ip < np_cel; ip++)
            {
                TPZManVector<REAL,3> qsi(el_dim,0.0);
                STATE w;
                int_points.Point(ip, qsi, w);
                // Get the vectorial phi
                q_intel->Shape(qsi, phi, dphidxi);
                q_intel->InitMaterialData(data);
                q_intel->ComputeRequiredData(data,qsi);
                
                TPZFMatrix<STATE> dphi       = data.dphi;
                REAL JacobianDet = data.detjac;
                
                TPZFMatrix<STATE> Qaxes = data.axes;
                TPZFMatrix<STATE> QaxesT;
                TPZFMatrix<STATE> Jacobian = data.jacobian;
                TPZFMatrix<STATE> JacobianInverse = data.jacinv;
                
                TPZFMatrix<STATE> GradOfX;
                TPZFMatrix<STATE> GradOfXInverse;
                TPZFMatrix<STATE> VectorOnMaster;
                TPZFMatrix<STATE> VectorOnXYZ(3,1,0.0);
                Qaxes.Transpose(&QaxesT);
                QaxesT.Multiply(Jacobian, GradOfX);
                JacobianInverse.Multiply(Qaxes, GradOfXInverse);

                if (el_dim == dim) {
                    for (int jp = 0; jp < q_block_dim.second; jp++) {
                        int vector_index = data.fVecShapeIndex[jp].first;
                        int shape_index = data.fVecShapeIndex[jp].second;
                        
                        for (int k = 0; k < dim; k++) {
                            VectorOnXYZ(k,0) = data.fNormalVec(k,vector_index);
                        }

                        GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
                        VectorOnMaster *= JacobianDet;
                        
                        for (int id = 0; id < dim; id++) {
                            block_phi_q(ip*dim+id,jp) = phi(shape_index,0)*VectorOnXYZ(id,0);
                            block_div_phi_q(ip,jp) += dphi(id,shape_index)*VectorOnMaster(id,0);
                        }
                    }
                }
                else{
                    for (int jp = 0; jp < q_block_dim.second; jp++) {
                        block_phi_q(ip,jp) = phi(jp,0);
                    }
                }
            }
            
        }
        
        
        // Pressure functions
        TPZInterpolationSpace * p_intel = dynamic_cast<TPZInterpolationSpace * >(mf_p_cel->Element(1));
        
        if(p_intel)
        {
            // for derivatives in real space
            int nshape = p_intel->NShapeF();
            TPZFNMatrix<220> phi(nshape,1);
            TPZFNMatrix<660> dphi(gel_dim,nshape);
            
            for (int ip = 0; ip < np_cel ; ip++)
            {
                TPZManVector<REAL,3> qsi(gel_dim,0.0);
                STATE w;
                int_points.Point(ip, qsi, w);
                p_intel->Shape(qsi, phi, dphi);
                
#ifdef PZDEBUG
                if(p_block_dim.second != phi.Rows()){
                    DebugStop();
                }
#endif
                for (int jp = 0; jp < phi.Rows(); jp++) {
                    block_phi_p(ip,jp) = phi(jp,0);
                }
                
            }
        }
        
        fq_To_parabolic.SetBlock(iel, block_phi_q);
        fdiv_q_To_parabolic.SetBlock(iel, block_phi_q);
        fp_To_parabolic.SetBlock(iel, block_phi_p);
        
    }

    
//    fq_To_parabolic.Print(" q_to_p ");
//    fdiv_q_To_parabolic.Print(" div_q_to_p ");
//    fp_To_parabolic.Print(" p_to_p ");
    
    return;

}

void TRMBuildTransfers::space_To_parabolic(TPZCompMesh * parabolic){
    
#ifdef PZDEBUG
    if (!parabolic || fq_To_parabolic.Rows() == 0 || fdiv_q_To_parabolic.Rows() == 0 || fp_To_parabolic.Rows() == 0) {
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = parabolic->Dimension();
    int n_el = fp_p_cindexes.size();
    
    std::pair<long, long> b_size_phi_q, b_size_div_phi_q, b_size_phi_p;

    b_size_phi_q.first = 0;
    b_size_phi_q.second = 0;
    
    b_size_div_phi_q.first = 0;
    b_size_div_phi_q.second = 0;
    
    b_size_phi_p.first = 0;
    b_size_phi_p.second = 0;
    
    TPZFMatrix<STATE> block_phi_q, block_div_phi_q, block_phi_p;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fp_p_cindexes[iel].second);
        
#ifdef PZDEBUG
        if (!p_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
#ifdef PZDEBUG
        if(!mf_p_cel)
        {
            DebugStop();
        }
#endif
        
        b_size_phi_q        = fq_To_parabolic.GetSizeofBlock(iel);
        b_size_div_phi_q    = fdiv_q_To_parabolic.GetSizeofBlock(iel);
        b_size_phi_p        = fp_To_parabolic.GetSizeofBlock(iel);

        fq_To_parabolic.GetBlock(iel, block_phi_q);
        fdiv_q_To_parabolic.GetBlock(iel, block_div_phi_q);
        fp_To_parabolic.GetBlock(iel, block_phi_p);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        TPZMaterial * material = parabolic->FindMaterial(matd_id);
        
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            int n_phi_q = b_size_phi_q.second;
            int n_phi_p = b_size_phi_p.second;
            
            TPZFNMatrix<3,STATE> phi_q(n_phi_q,dim,0.0);
            TPZFNMatrix<3,STATE> div_phi_q(n_phi_q,1,0.0);
            TPZFNMatrix<3,STATE> phi_p(n_phi_p,1,0.0);
            
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int is = 0; is < n_phi_q; is++) {
                    for (int d = 0; d < dim; d++) {
                        phi_q(is,d) = block_phi_q(ip*dim+d,is);
                    }
                    div_phi_q(is,0) = block_phi_q(ip,is);
                }
                
                for (int is = 0; is < n_phi_p; is++) {
                    phi_p(is,0) = block_phi_p(ip,is);
                }
                
                associated_material->GetMemory()[ipos].Set_phi_q(phi_q);
                associated_material->GetMemory()[ipos].Set_div_phi_q(div_phi_q);
                associated_material->GetMemory()[ipos].Set_phi_p(phi_p);
                
            }
            
        }
        else{
            
            TPZMatWithMem<TRMMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            int n_phi_q = b_size_phi_q.second;
            
            TPZFNMatrix<3,STATE> phi_q(n_phi_q,1,0.0);
            
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int is = 0; is < n_phi_q; is++) {
                    phi_q(is,0) = block_phi_q(ip,is);
                }
            
                associated_material->GetMemory()[ipos].Set_phi_q(phi_q);
                
            }
            
        }
        
        
    }
    
}

void TRMBuildTransfers::kappa_phi_To_parabolic(TPZCompMesh * parabolic){
    
    
#ifdef PZDEBUG
    if (!parabolic) {
        DebugStop();
    }
#endif
    
    int dim = parabolic->Dimension();
    
    
    TPZFMatrix<STATE> kappa, kappa_inv;
    TPZManVector<STATE, 10> vars;
    TPZManVector<STATE, 10> porosity;
    
    // Step one
    int n_elements = parabolic->NElements();
    TPZManVector<long, 30> indexes;
    for (int icel = 0; icel < n_elements; icel++) {
        TPZCompEl * cel = parabolic->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        
#ifdef PZDEBUG
        if (!mf_cel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension()!= dim) {
            continue;
        }
        
        const TPZIntPoints & int_points = mf_cel->GetIntegrationRule();
        int np = int_points.NPoints();
        GlobalPointIndexes(cel, indexes);
        
#ifdef PZDEBUG
        if (indexes.size() != np) {
            DebugStop();
        }
#endif
        
        int rockid = gel->MaterialId();
        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = parabolic->FindMaterial(rockid);
        TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
        
        TPZManVector<REAL,3> par_triplet(3,0.0);
        TPZManVector<REAL,3> x(3,0.0);
        REAL w;
        for (int ip = 0; ip<np; ip++) {
            int_points.Point(ip, par_triplet, w);
            gel->X(par_triplet, x);
            
            associated_material->GetMemory()[indexes[ip]].Set_x(x);
            //            fSimulationData->Map()->ComputePropertieSPE10Map(index, x, kappa, kappa_inv, porosity);
            fSimulationData->Map()->Kappa(x, kappa, kappa_inv, vars);
            fSimulationData->Map()->phi(x, porosity, vars);
            associated_material->GetMemory()[indexes[ip]].Set_K_0(kappa);
            associated_material->GetMemory()[indexes[ip]].Set_Kinv_0(kappa_inv);
            associated_material->GetMemory()[indexes[ip]].Set_phi_0(porosity[0]);
        }
    }
    
}


void TRMBuildTransfers::parabolic_To_parabolic(TPZCompMesh * parabolic){

#ifdef PZDEBUG
    if (!parabolic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_q(fq_To_parabolic.Cols(),1,0.0);
    TPZFMatrix<STATE> Scatter_p(fp_To_parabolic.Cols(),1,0.0);
    
    int n = fp_p_cindexes.size();
    long pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fq_dof_scatter[i].size(); iequ++) {
            Scatter_q(pos,0) = parabolic->Solution()(fq_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fp_dof_scatter[i].size(); iequ++) {
            Scatter_p(pos,0) = parabolic->Solution()(fp_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> q_at_parabolic,div_q_at_parabolic,p_at_parabolic;
    fq_To_parabolic.Multiply(Scatter_q,q_at_parabolic);
    fdiv_q_To_parabolic.Multiply(Scatter_q,div_q_at_parabolic);
    fp_To_parabolic.Multiply(Scatter_p,p_at_parabolic);

    
    
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = parabolic->Dimension();
    int n_el = fp_p_cindexes.size();
    
    long first_point_phi_q = 0;
    long first_point_div_phi_q = 0;
    long first_point_phi_p = 0;
    
    std::pair<long, long> b_size_phi_q, b_size_div_phi_q, b_size_phi_p;
    
    b_size_phi_q.first = 0;
    b_size_phi_q.second = 0;
    
    b_size_div_phi_q.first = 0;
    b_size_div_phi_q.second = 0;
    
    b_size_phi_p.first = 0;
    b_size_phi_p.second = 0;
    
    TPZFMatrix<STATE> block_phi_q, block_div_phi_q, block_phi_p;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fp_p_cindexes[iel].second);
        
#ifdef PZDEBUG
        if (!p_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
#ifdef PZDEBUG
        if(!mf_p_cel)
        {
            DebugStop();
        }
#endif
        
        first_point_phi_q     += b_size_phi_q.first;
        first_point_div_phi_q += b_size_div_phi_q.first;
        first_point_phi_p     += b_size_phi_p.first;
        
        b_size_phi_q        = fq_To_parabolic.GetSizeofBlock(iel);
        b_size_div_phi_q    = fdiv_q_To_parabolic.GetSizeofBlock(iel);
        b_size_phi_p        = fp_To_parabolic.GetSizeofBlock(iel);
        
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        TPZMaterial * material = parabolic->FindMaterial(matd_id);
        
        if(gel->Dimension() == dim){ // The volumetric ones!
        
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            TPZManVector<REAL,3> q(3,0.0);
            STATE div_q, p;
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    q[id] = q_at_parabolic(first_point_phi_q + ip*dim + id,0);
                }
                
                div_q   = div_q_at_parabolic(first_point_div_phi_q + ip,0);
                p       = p_at_parabolic(first_point_phi_p + ip,0);
                
                if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                    associated_material->GetMemory()[ipos].Set_p_0(p);
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_q_n(q);
                    associated_material->GetMemory()[ipos].Set_div_q_n(div_q);
                    associated_material->GetMemory()[ipos].Set_p_n(p);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_q(q);
                    associated_material->GetMemory()[ipos].Set_div_q(div_q);
                    associated_material->GetMemory()[ipos].Set_p(p);
                }
                
            }
            
            
        }
        else{
            

            TPZMatWithMem<TRMMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            TPZManVector<REAL,3> q(1,0.0);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                q[0]       = q_at_parabolic(first_point_phi_q + ip,0);
                
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_q_n(q);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_q(q);
                }
                
            }
            
        }
        
    }


}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Segregated Transfer methods (Gamma and Omega) :: Build methods Hyperbolic
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void TRMBuildTransfers::Build_hyperbolic_To_hyperbolic(TPZCompMesh * hyperbolic){
    
    DebugStop();
}

void TRMBuildTransfers::space_To_hyperbolic(TPZCompMesh * hyperbolic){
    
    DebugStop();
}

void TRMBuildTransfers::hyperbolic_To_hyperbolic(TPZCompMesh * hyperbolic){
    
    DebugStop();
}


////////////////////////// Transfers:: Iterative Coupling by Operator Splitting //////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////









































/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Matrices Initialization Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void TRMBuildTransfers::Initialize_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = cmesh_multiphysics->Reference()->Dimension(); // vectorial
    long element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fu_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions(nel);
    
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
#ifdef PZDEBUG
        if(!mf_cel)
        {
            DebugStop();
        }
#endif
        element_index = mf_cel->Index();
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        if(intel->Dimension() < n_var_dim){
            // there is boundary elements for normal flux where it is a scalar variable
//            mf_cel->GetMemoryIndices(int_point_indexes);
//            this->ElementDofIndexes(intel, dof_indexes);
//            fu_dof_scatter[element_index] = dof_indexes;
            blocks_dimensions[element_index].first = 0;
            blocks_dimensions[element_index].second = 0;
            fu_dof_scatter[element_index] = dof_indexes;
            continue;
        }
        
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel, dof_indexes);
        fu_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions[element_index].second = dof_indexes.size();
        fu_dof_scatter[element_index] = dof_indexes;
    }
    
    // Initialize the matrix
    fu_To_Mixed.Initialize(blocks_dimensions);
    
}


void TRMBuildTransfers::Initialize_p_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    long element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fp_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions(nel);
    
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
#ifdef PZDEBUG
        if(!mf_cel)
        {
            DebugStop();
        }
#endif
        element_index = mf_cel->Index();
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        if(!intel){
            // there is no boundary elements for pressure
            blocks_dimensions[element_index].first = 0*n_var_dim;
            blocks_dimensions[element_index].second = 0;
            fp_dof_scatter[element_index] = dof_indexes;
            continue;
        }
        
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel, dof_indexes);
        fp_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions[element_index].second = dof_indexes.size();
        fp_dof_scatter[element_index] = dof_indexes;
    }
    
    // Initialize the matrix
    fp_To_Mixed.Initialize(blocks_dimensions);
    
    
}


void TRMBuildTransfers::Initialize_s_To_Transport(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    long element_index = 0;
    int dimension = cmesh_multiphysics->Reference()->Dimension();
    // Compute destination index scatter by element (Omega and Gamma)
    if (!mesh_index) {
        fsa_dof_scatter.Resize(nel);
    }
    else{
        fsb_dof_scatter.Resize(nel);
    }

    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions(nel);
    
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
        if(gel->HasSubElement()){
            continue;
        }
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
       
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);        
        
#ifdef PZDEBUG
        if(!mf_cel)
        {
            if(!mf_int_cel){
                DebugStop();
            }
        }
#endif
        if (mf_int_cel) {
            continue;
        }
        
        element_index = mf_cel->Index();
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        if(intel->Dimension() < dimension){
            // there is no boundary elements for saturation
            blocks_dimensions[element_index].first = 0*n_var_dim;
            blocks_dimensions[element_index].second = 0;
            if (!mesh_index) {
                fsa_dof_scatter[element_index] = dof_indexes;
            }
            else{
                fsb_dof_scatter[element_index] = dof_indexes;
            }
            
            continue;
        }
        
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel, dof_indexes);

        if (!mesh_index) {
            fsa_dof_scatter[element_index] = dof_indexes;
        }
        else{
            fsb_dof_scatter[element_index] = dof_indexes;
        }
        
        blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions[element_index].second = dof_indexes.size();
    }
    
    // Initialize the matrix
    fs_To_Transport.Initialize(blocks_dimensions);
    
    
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Matrices Filling Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TRMBuildTransfers::Fill_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast operations and mesh structure, and  finally it initialize diagonal matrix blocks
    Initialize_u_To_Mixed(cmesh_multiphysics, mesh_index);
    
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = cmesh_multiphysics->Reference()->Dimension();; // vector
    long element_index = 0;
    
    TPZMaterialData data;
    
    std::pair<long, long> block_dim;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        dof_indexes = fu_dof_scatter[element_index];
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_cel->GetIntegrationRule();
        int np_cel = int_points_mixed.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel->NShapeF(),0.0);
        TPZFMatrix<double> block;
        
        if(intel->Dimension() < n_var_dim){ // lower dimensional elements
            block.Resize(block_dim.first,block_dim.second);        }
        else{
            block.Resize(block_dim.first*n_var_dim,block_dim.second);
        }
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(el_dim,0.0);
            STATE w;
            int_points_mixed.Point(ip, qsi, w);
            // Get the vectorial phi
            intel->Shape(qsi, phi, dphidxi);
            intel->InitMaterialData(data);
            intel->ComputeRequiredData(data,qsi);
            
            for (int id = 0; id < n_var_dim; id++) {
                for (int jp = 0; jp < block_dim.second; jp++) {
                    int vector_index = data.fVecShapeIndex[jp].first;
                    int shape_index = data.fVecShapeIndex[jp].second;
                    block(ip*n_var_dim+id,jp) = phi(shape_index,0)*data.fNormalVec(id,vector_index);
                }
            }
            
        }
        
        fu_To_Mixed.SetBlock(element_index, block);
        
    }
    
    return;
}


void TRMBuildTransfers::Fill_p_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_p_To_Mixed(cmesh_multiphysics, mesh_index);
    
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    long element_index = 0;
    
    std::pair<long, long> block_dim;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        if(!intel){
            continue;
        }
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        dof_indexes = fp_dof_scatter[element_index];
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_cel->GetIntegrationRule();
        int np_cel = int_points_mixed.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel->NShapeF(),0.0);
        TPZFMatrix<double> block(block_dim.first*n_var_dim,block_dim.second);
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(el_dim,0.0);
            STATE w;
            int_points_mixed.Point(ip, qsi, w);
            intel->Shape(qsi, phi, dphidxi);
            
            for (int id = 0; id < n_var_dim; id++) {
                for (int jp = 0; jp < block_dim.second; jp++) {
                    block(ip+id*n_var_dim,jp) = phi(jp,0);
                }
            }
            
        }
        
        fp_To_Mixed.SetBlock(element_index, block);
        
    }
    
    return;
}



void TRMBuildTransfers::Fill_s_To_Transport(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_s_To_Transport(cmesh_multiphysics, mesh_index);
    
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    long element_index = 0;
    int dimension = cmesh_multiphysics->Reference()->Dimension();
    std::pair<long, long> block_dim;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        
        TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);
        
        if (mf_int_cel) {
            continue;
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        if(intel->Dimension() < dimension){
            continue;
        }
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        if (!mesh_index) {
            dof_indexes = fsa_dof_scatter[element_index];
        }
        else{
            dof_indexes = fsb_dof_scatter[element_index];
        }
        
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_cel->GetIntegrationRule();
        int np_cel = int_points_mixed.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel->NShapeF(),0.0);
        TPZFMatrix<double> block(block_dim.first*n_var_dim,block_dim.second);
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(el_dim,0.0);
            STATE w;
            int_points_mixed.Point(ip, qsi, w);
            intel->Shape(qsi, phi, dphidxi);
            
            for (int id = 0; id < n_var_dim; id++) {
                for (int jp = 0; jp < block_dim.second; jp++) {
                    block(ip+id*n_var_dim,jp) = phi(jp,0);
                }
            }
            
        }
        
        fs_To_Transport.SetBlock(element_index, block);
        
    }
    return;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Transfer Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void TRMBuildTransfers::kappa_phi_To_Mixed_Memory(TPZCompMesh * cmesh_multiphysics){
    
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int dim = cmesh_multiphysics->Dimension();
    
    
    TPZFMatrix<STATE> kappa, kappa_inv;
    TPZManVector<STATE, 10> vars;
    TPZManVector<STATE, 10> porosity;

    // Step one
    int n_elements = cmesh_multiphysics->NElements();
    TPZManVector<long, 30> indexes;
    for (int icel = 0; icel < n_elements; icel++) {
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        
#ifdef PZDEBUG
        if (!mf_cel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension()!= dim) {
            continue;
        }
        
        const TPZIntPoints & int_points = mf_cel->GetIntegrationRule();
        int np = int_points.NPoints();
        GlobalPointIndexes(cel, indexes);
        
#ifdef PZDEBUG
        if (indexes.size() != np) {
            DebugStop();
        }
#endif
        
        int rockid = gel->MaterialId();
        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
        TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
        
        TPZManVector<REAL,3> par_triplet(3,0.0);
        TPZManVector<REAL,3> x(3,0.0);
        REAL w;
        for (int ip = 0; ip<np; ip++) {
            int_points.Point(ip, par_triplet, w);
            gel->X(par_triplet, x);
            
            associated_material->GetMemory()[indexes[ip]].Set_x(x);
//            fSimulationData->Map()->ComputePropertieSPE10Map(index, x, kappa, kappa_inv, porosity);
            fSimulationData->Map()->Kappa(x, kappa, kappa_inv, vars);
            fSimulationData->Map()->phi(x, porosity, vars);
            associated_material->GetMemory()[indexes[ip]].Set_K_0(kappa);
            associated_material->GetMemory()[indexes[ip]].Set_Kinv_0(kappa_inv);
            associated_material->GetMemory()[indexes[ip]].Set_phi_0(porosity[0]);
        }
    }
    
}


void TRMBuildTransfers::u_To_Mixed_Memory(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_multiphysics){
    

#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int nel = cmesh_multiphysics->NElements();
    int dim = cmesh_flux->Dimension();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterFlux(fu_To_Mixed.Cols(),1,0.0);
    long pos = 0;
    for (int el = 0; el < nel; el++) {
        for(int ip = 0; ip < fu_dof_scatter[el].size(); ip++) {
            ScatterFlux(pos,0) = cmesh_flux->Solution()(fu_dof_scatter[el][ip],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> Flux_at_intpoints;
    fu_To_Mixed.Multiply(ScatterFlux,Flux_at_intpoints);
    // Trasnfering values
    TPZManVector<STATE,3> u(dim,0.0);
    for(long i = 0; i <  np_cmesh; i++){
        for (int id = 0; id < dim ; id++) {
            u[id]= Flux_at_intpoints(i*dim+id,0);
        }

        if(fSimulationData->IsCurrentStateQ()){
            associated_material->GetMemory()[i].Set_u(u);
        }
        else{
            associated_material->GetMemory()[i].Set_u_n(u);
        }
        
    }
    
}

void TRMBuildTransfers::p_To_Mixed_Memory(TPZCompMesh * cmesh_pressure, TPZCompMesh * cmesh_multiphysics){

    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int nel = cmesh_multiphysics->NElements();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];

    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterPressure(fp_To_Mixed.Cols(),1,0.0);
    long pos = 0;
    for (int el = 0; el < nel; el++) {
        for(int ip = 0; ip < fp_dof_scatter[el].size(); ip++) {
            ScatterPressure(pos,0) = cmesh_pressure->Solution()(fp_dof_scatter[el][ip],0);
            pos++;
        }
    }
    
    // Step two
    TPZFNMatrix<30,STATE> Pressure_at_intpoints;
    fp_To_Mixed.Multiply(ScatterPressure,Pressure_at_intpoints);
    // Trasnfering values
    for(long i = 0; i <  np_cmesh; i++){
        if(fSimulationData->IsCurrentStateQ()){
            associated_material->GetMemory()[i].Set_p_n(Pressure_at_intpoints(i,0));
        }
        else{
            associated_material->GetMemory()[i].Set_p(Pressure_at_intpoints(i,0));
        }
    }
    
}

void TRMBuildTransfers::kappa_phi_To_Transport_Memory(TPZCompMesh * cmesh_multiphysics){
    
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int dim = cmesh_multiphysics->Dimension();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(material);
    

    TPZFMatrix<STATE> kappa, kappa_inv;
    TPZManVector<STATE, 10> vars;
    TPZManVector<STATE, 10> porosity;
//    REAL porosity;
//    long index = 0;
    // Step one
    int n_elements = cmesh_multiphysics->NElements();
    TPZManVector<long, 30> indexes;
    for (int icel = 0; icel < n_elements; icel++) {
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension()!= dim) {
            continue;
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        
#ifdef PZDEBUG
        if (!mf_cel) {
            DebugStop();
        }
#endif
        
        
        const TPZIntPoints & int_points = mf_cel->GetIntegrationRule();
        int np = int_points.NPoints();
        GlobalPointIndexes(cel, indexes);
        
#ifdef PZDEBUG
        if (indexes.size() != np) {
            DebugStop();
        }
#endif
        
        TPZManVector<REAL,3> par_triplet(3,0.0);
        TPZManVector<REAL,3> x(3,0.0);
        REAL w;
        for (int ip = 0; ip<np; ip++) {
            int_points.Point(ip, par_triplet, w);
            gel->X(par_triplet, x);
            
            associated_material->GetMemory()[indexes[ip]].Set_x(x);
            //            fSimulationData->Map()->ComputePropertieSPE10Map(index, x, kappa, kappa_inv, porosity);
            fSimulationData->Map()->Kappa(x, kappa, kappa_inv, vars);
            fSimulationData->Map()->phi(x, porosity, vars);
            associated_material->GetMemory()[indexes[ip]].Set_K_0(kappa);
            associated_material->GetMemory()[indexes[ip]].Set_Kinv_0(kappa_inv);
            associated_material->GetMemory()[indexes[ip]].Set_phi_0(porosity[0]);
        }
    }
    
}


void TRMBuildTransfers::s_To_Transport_Memory(TPZCompMesh * cmesh_saturation, TPZCompMesh * cmesh_multiphysics, int mesh_index){

    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int nel = cmesh_multiphysics->NElements();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterSaturation(fs_To_Transport.Cols(),1,0.0);
    long pos = 0;
    for (int el = 0; el < nel; el++) {
        
        if (!mesh_index) {
            for(int ip = 0; ip < fsa_dof_scatter[el].size(); ip++) {
                ScatterSaturation(pos,0) = cmesh_saturation->Solution()(fsa_dof_scatter[el][ip],0);
                pos++;
            }
        }
        else{
            for(int ip = 0; ip < fsb_dof_scatter[el].size(); ip++) {
                ScatterSaturation(pos,0) = cmesh_saturation->Solution()(fsb_dof_scatter[el][ip],0);
                pos++;
            }
        }
        
    }
    
    // Step two
    TPZFMatrix<STATE> Saturation_at_intpoints;
    fs_To_Transport.Multiply(ScatterSaturation,Saturation_at_intpoints);
    // Trasnfering values
    for(long i = 0; i <  np_cmesh; i++){
        
        if(!mesh_index){
            if(fSimulationData->IsCurrentStateQ()){
                associated_material->GetMemory()[i].Set_sa_n(Saturation_at_intpoints(i,0));
            }
            else{
                associated_material->GetMemory()[i].Set_sa(Saturation_at_intpoints(i,0));
            }
        }
        else{
            if(fSimulationData->IsCurrentStateQ()){
                associated_material->GetMemory()[i].Set_sb_n(Saturation_at_intpoints(i,0));
            }
            else{
                associated_material->GetMemory()[i].Set_sb(Saturation_at_intpoints(i,0));
            }
        }
        
    }
    
}

/** @brief Reciprocal (mixed <-> transpor) transfer average quantities to integration points of multiphysics meshes over volumetric elements */
void TRMBuildTransfers::Reciprocal_Memory_Transfer(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_trans){
    
    
#ifdef PZDEBUG
    if ( fmixed_transport_cindexes.size() == 0 ) {
        DebugStop();
    }

    if (!cmesh_mf_mixed || !cmesh_mf_trans) {
        DebugStop();
    }
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * mixed_material = cmesh_mf_mixed->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(mixed_material);
    
    TPZMaterial * trans_material = cmesh_mf_trans->FindMaterial(rockid);
    TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * trans_memory = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(trans_material);
    
    TPZManVector<long,30> p_point_indexes;
    TPZManVector<long,30> s_point_indexes;
    long nvolumes = fmixed_transport_cindexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
    TPZGeoEl  * gel = geometry->Element(fmixed_transport_cindexes[ivol].first);
    TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_cindexes[ivol].second.first);
    TPZCompEl * trans_cel = cmesh_mf_trans->Element(fmixed_transport_cindexes[ivol].second.second);
        
#ifdef PZDEBUG
        if (!mixed_cel || !trans_cel || !gel) {
            DebugStop();
        }
#endif

        REAL element_measure = DimensionalMeasure(mixed_cel->Reference());
    
        GlobalPointIndexes(mixed_cel, p_point_indexes);
        GlobalPointIndexes(trans_cel, s_point_indexes);
        
        TPZMultiphysicsElement * mf_mixed_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        TPZMultiphysicsElement * mf_trans_cel = dynamic_cast<TPZMultiphysicsElement * >(trans_cel);
        
#ifdef PZDEBUG
        if (!mf_mixed_cel || !mf_trans_cel) {
            DebugStop();
        }
#endif
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_mixed_cel->GetIntegrationRule();
        int np_mixed_cel = int_points_mixed.NPoints();
        
        const TPZIntPoints & int_points_trans = mf_trans_cel->GetIntegrationRule();
        int np_trans_cel = int_points_trans.NPoints();
        
#ifdef PZDEBUG
        if (np_mixed_cel != p_point_indexes.size() || np_trans_cel != s_point_indexes.size()) {
            DebugStop();
        }
#endif
        
        REAL w;
        TPZManVector<REAL,3> triplet(3,0.0);
        
        REAL detjac;
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        TPZFMatrix<REAL> jacinv;

        REAL p_avg      = 0.0;
        REAL p_avg_n    = 0.0;

        // Integrating pressure
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            int_points_mixed.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            p_avg_n += w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p_n()/element_measure;
            p_avg +=  w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p()/element_measure;
            
        }

        REAL sa      = 0.0;
        REAL sa_n    = 0.0;
        REAL sb      = 0.0;
        REAL sb_n    = 0.0;
        
        // Integrating Saturation
        for (int ip = 0; ip < np_trans_cel; ip++) {
            
            int_points_trans.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            sa_n += w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sa_n()/element_measure;
            sb_n += w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sb_n()/element_measure;
            sa +=  w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sa()/element_measure;
            sb +=  w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sb()/element_measure;

        }
        
//        std::cout << "p_avg_n = "<< p_avg_n << std::endl;
//        std::cout << "p_avg = "<< p_avg << std::endl;
//        
//        std::cout << "sa_n = "<< sa_n << std::endl;
//        std::cout << "sa = "<< sa << std::endl;
        
        // Inserting average pressure and saturation in mixed memory
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            if (fSimulationData->IsCurrentStateQ()) {
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg_n(p_avg_n);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sa_n(sa_n);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sb_n(sb_n);
            }
            else{
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg(p_avg);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sa(sa);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sb(sb);
            }

        }

        // Inserting average pressure in transport memory
        for (int ip = 0; ip < np_trans_cel; ip++) {
            
            if (fSimulationData->IsCurrentStateQ()) {
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_p_avg_n(p_avg_n);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa_n(sa_n);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa_n(sb_n);
            }
            else{
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_p_avg(p_avg);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa(sa);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa(sb);
            }

        }
        
    }
    
    return;
    
}

/** @brief Reciprocal (mixed <-> transpor) transfer average quantities to integration points of multiphysics meshes over volumetric elements */
void TRMBuildTransfers::Reciprocal_Memory_TransferII(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_trans){
    
    
#ifdef PZDEBUG
    if ( fmixed_transport_comp_indexes.size() == 0 ) {
        DebugStop();
    }
    
    if (!cmesh_mf_mixed || !cmesh_mf_trans) {
        DebugStop();
    }
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * mixed_material = cmesh_mf_mixed->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(mixed_material);
    
    TPZMaterial * trans_material = cmesh_mf_trans->FindMaterial(rockid);
    TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * trans_memory = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(trans_material);
    
    TPZManVector<long,30> p_point_indexes;
    TPZManVector<long,30> s_point_indexes;
    long nvolumes = fmixed_transport_comp_indexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
        TPZGeoEl  * gel = geometry->Element(fmixed_transport_comp_indexes[ivol].first);
        TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_comp_indexes[ivol].second.first);
        
        // for each transport subelement
        int n_subcels = fmixed_transport_comp_indexes[ivol].second.second.size();
        
        REAL avg_sa      = 0.0;
        REAL avg_sa_n    = 0.0;
        REAL avg_sb      = 0.0;
        REAL avg_sb_n    = 0.0;
        
        for (int icel = 0; icel < n_subcels; icel++) {
            TPZCompEl * trans_cel = cmesh_mf_trans->Element( fmixed_transport_comp_indexes[ivol].second.second[icel]);
            
#ifdef PZDEBUG
            if (!mixed_cel || !trans_cel || !gel) {
                DebugStop();
            }
#endif
            
            REAL element_measure_mixed = DimensionalMeasure(mixed_cel->Reference());
            REAL element_measure_transport = DimensionalMeasure(trans_cel->Reference());
            
            GlobalPointIndexes(mixed_cel, p_point_indexes);
            GlobalPointIndexes(trans_cel, s_point_indexes);
            
            TPZMultiphysicsElement * mf_mixed_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
            TPZMultiphysicsElement * mf_trans_cel = dynamic_cast<TPZMultiphysicsElement * >(trans_cel);
            
#ifdef PZDEBUG
            if (!mf_mixed_cel || !mf_trans_cel) {
                DebugStop();
            }
#endif
            
            // Computing the local integration points indexes
            const TPZIntPoints & int_points_mixed = mf_mixed_cel->GetIntegrationRule();
            int np_mixed_cel = int_points_mixed.NPoints();
            
            const TPZIntPoints & int_points_trans = mf_trans_cel->GetIntegrationRule();
            int np_trans_cel = int_points_trans.NPoints();
            
#ifdef PZDEBUG
            if (np_mixed_cel != p_point_indexes.size() || np_trans_cel != s_point_indexes.size()) {
                DebugStop();
            }
#endif
            
            REAL w;
            TPZManVector<REAL,3> triplet(3,0.0);
            
            REAL detjac;
            TPZFMatrix<REAL> jac;
            TPZFMatrix<REAL> axes;
            TPZFMatrix<REAL> jacinv;
            
            REAL p_avg      = 0.0;
            REAL p_avg_n    = 0.0;
            
            // Integrating pressure
            for (int ip = 0; ip < np_mixed_cel; ip++) {
                int_points_mixed.Point(ip, triplet, w);
                gel->Jacobian(triplet, jac, axes, detjac, jacinv);
                
                p_avg_n += w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p_n()/element_measure_mixed;
                p_avg +=  w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p()/element_measure_mixed;
                
            }
            
            TPZGeoEl  * gel_transport = mf_trans_cel->Reference();
            
            REAL sa      = 0.0;
            REAL sa_n    = 0.0;
            REAL sb      = 0.0;
            REAL sb_n    = 0.0;
            
            // Integrating Saturation
            for (int ip = 0; ip < np_trans_cel; ip++) {
                
                int_points_trans.Point(ip, triplet, w);
                gel_transport->Jacobian(triplet, jac, axes, detjac, jacinv);
                
                sa_n += w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sa_n()/element_measure_transport;
                sb_n += w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sb_n()/element_measure_transport;
                sa +=  w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sa()/element_measure_transport;
                sb +=  w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sb()/element_measure_transport;
                
            }
            
            if (fSimulationData->IsCurrentStateQ()) {
                avg_sa_n += sa_n*element_measure_transport/element_measure_mixed;
                avg_sb_n += sb_n*element_measure_transport/element_measure_mixed;
            }
            else{
                avg_sa += sa*element_measure_transport/element_measure_mixed;
                avg_sb += sb*element_measure_transport/element_measure_mixed;
            }
            
            // Inserting average pressure and saturation in mixed memory
            for (int ip = 0; ip < np_mixed_cel; ip++) {
                if (fSimulationData->IsCurrentStateQ()) {
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg_n(p_avg_n);
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sa_n(avg_sa_n);
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sb_n(avg_sb_n);
                }
                else{
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg(p_avg);
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sa(avg_sa);
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sb(avg_sb);
                }
                
            }
            
            // Inserting average pressure in transport memory
            for (int ip = 0; ip < np_trans_cel; ip++) {
                
                if (fSimulationData->IsCurrentStateQ()) {
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_p_avg_n(p_avg_n);
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa_n(sa_n);
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa_n(sb_n);
                }
                else{
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_p_avg(p_avg);
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa(sa);
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa(sb);
                }
                
            }
        }
        
    }
    
    return;
    
}

/** @brief Transfer average pressure to integration points of multiphysics mixed meshes over volumetric elements */
void TRMBuildTransfers::p_avg_Memory_Transfer(TPZCompMesh * cmesh_mf_mixed){
    
    
#ifdef PZDEBUG
    if ( fmixed_transport_cindexes.size() == 0 ) {
        DebugStop();
    }
    
    if (!cmesh_mf_mixed) {
        DebugStop();
    }
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    
    TPZManVector<long,30> p_point_indexes;
    long nvolumes = fmixed_transport_cindexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
        TPZGeoEl  * gel = geometry->Element(fmixed_transport_cindexes[ivol].first);
        TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_cindexes[ivol].second.first);
        
#ifdef PZDEBUG
        if (!mixed_cel || !gel) {
            DebugStop();
        }
#endif
        
        REAL element_measure = DimensionalMeasure(mixed_cel->Reference());
        
        GlobalPointIndexes(mixed_cel, p_point_indexes);
        TPZMultiphysicsElement * mf_mixed_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        
#ifdef PZDEBUG
        if (!mf_mixed_cel) {
            DebugStop();
        }
#endif
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_mixed_cel->GetIntegrationRule();
        int np_mixed_cel = int_points_mixed.NPoints();
        
        
#ifdef PZDEBUG
        if (np_mixed_cel != p_point_indexes.size()) {
            DebugStop();
        }
#endif
        
        int rockid = gel->MaterialId();
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * mixed_material = cmesh_mf_mixed->FindMaterial(rockid);
        TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(mixed_material);
        
        REAL w;
        TPZManVector<REAL,3> triplet(3,0.0);
        
        REAL detjac;
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        TPZFMatrix<REAL> jacinv;
        
        REAL p_avg      = 0.0;
        REAL p_avg_n    = 0.0;
        
        // Integrating pressure
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            int_points_mixed.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            p_avg_n +=  w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p_n()/element_measure;
            p_avg += w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p()/element_measure;
            
        }
        
        // Inserting average pressure in mixed memory
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            if (fSimulationData->IsCurrentStateQ()) {
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg_n(p_avg_n);
            }
            else{
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg(p_avg);

            }
            
        }
        
    }
        
}

/** @brief Transfer average pressure to integration points of multiphysics mixed meshes over volumetric elements */
void TRMBuildTransfers::p_avg_Memory_TransferII(TPZCompMesh * cmesh_mf_mixed){
    
    
#ifdef PZDEBUG
    if ( fmixed_transport_comp_indexes.size() == 0 ) {
        DebugStop();
    }
    
    if (!cmesh_mf_mixed) {
        DebugStop();
    }
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * mixed_material = cmesh_mf_mixed->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(mixed_material);
    
    TPZManVector<long,30> p_point_indexes;
    long nvolumes = fmixed_transport_comp_indexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
        TPZGeoEl  * gel = geometry->Element(fmixed_transport_comp_indexes[ivol].first);
        TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_comp_indexes[ivol].second.first);
        
#ifdef PZDEBUG
        if (!mixed_cel || !gel) {
            DebugStop();
        }
#endif
        
        REAL element_measure = DimensionalMeasure(mixed_cel->Reference());
        
        GlobalPointIndexes(mixed_cel, p_point_indexes);
        TPZMultiphysicsElement * mf_mixed_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        
#ifdef PZDEBUG
        if (!mf_mixed_cel) {
            DebugStop();
        }
#endif
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_mixed_cel->GetIntegrationRule();
        int np_mixed_cel = int_points_mixed.NPoints();
        
        
#ifdef PZDEBUG
        if (np_mixed_cel != p_point_indexes.size()) {
            DebugStop();
        }
#endif
        
        REAL w;
        TPZManVector<REAL,3> triplet(3,0.0);
        
        REAL detjac;
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        TPZFMatrix<REAL> jacinv;
        
        REAL p_avg      = 0.0;
        REAL p_avg_n    = 0.0;
        
        // Integrating pressure
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            int_points_mixed.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            p_avg_n +=  w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p_n()/element_measure;
            p_avg += w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p()/element_measure;
            
        }
        
        // Inserting average pressure in mixed memory
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            if (fSimulationData->IsCurrentStateQ()) {
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg_n(p_avg_n);
            }
            else{
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg(p_avg);
                
            }
            
        }
        
    }
    
}


/** @brief Initializate  diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Initialize_un_To_Transport(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ){
    
    
#ifdef PZDEBUG
    if (!flux_mesh || !transport_mesh) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    int mesh_index = 0;
    
    TPZGeoMesh * geometry = flux_mesh->Reference();
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<long,10> dof_indexes;
    
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;

    
    TPZManVector<long> indices;
    std::pair<long, long> duplet;
    TPZManVector<int,10> face_sides;
    long face_index;
    long n_interfaces;
    
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
        fun_dof_scatter_Gamma.Resize(n_interfaces);
        fun_To_Transport_Gamma.Resize(0, 0);
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
        fun_dof_scatter_gamma.Resize(n_interfaces);
        fun_To_Transport_gamma.Resize(0, 0);
    }

    
    // Block size structue (Gamma or gamma (Inner element interfaces))
    TPZVec< std::pair<long, long> > blocks_dimensions(n_interfaces);
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {


        if (IsBoundaryQ) {
            face_index  = finterface_g_indexes_Gamma[k_face];
            duplet      = fleft_right_g_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_g_indexes_gamma[k_face];
            duplet      = fleft_right_g_indexes_gamma[k_face];
        }
        
        face_gel = geometry->Element(face_index);
        
        if (!face_gel) {
            DebugStop();
        }
        
        long left_geo_index     = duplet.first;
        long right_geo_index    = duplet.second;
        
        left_gel    = geometry->Element(left_geo_index);
        right_gel   = geometry->Element(right_geo_index);
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        this->ComputeFaceIndex(left_gel,face_sides);
        
//        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(left_cel);
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
//        
//        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace *> (left_cel);
        
        
        int face_side     = -1;
        long connect_index = -1;

        if(!IdentifyFace(face_side,left_gel,face_gel)){
            std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
            DebugStop();
        }
        else{
            
            for(int ic = 0; ic < intel_vol->NConnects(); ic++){
                connect_index = intel_vol->SideConnectLocId(ic, face_side);
                if (connect_index != -1) {
                    break;
                }
            }
            
            if (connect_index == -1) {
                std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
                DebugStop();
            }
        }

        this->ElementDofFaceIndexes(connect_index, mf_cel, dof_indexes);
        
        int nshapes = left_cel->Connect(connect_index).NShape();
        
        blocks_dimensions[k_face].first = 1;
        blocks_dimensions[k_face].second = nshapes;
        if (IsBoundaryQ) {
            fun_dof_scatter_Gamma[k_face] = dof_indexes;
        }
        else{
            fun_dof_scatter_gamma[k_face] = dof_indexes;
        }

        
    }
    
    // Initialize the matrix
    
    if (IsBoundaryQ) {
        fun_To_Transport_Gamma.Initialize(blocks_dimensions);
    }
    else{
        fun_To_Transport_gamma.Initialize(blocks_dimensions);
    }
    
    
}

/** @brief Initializate  diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Initialize_un_To_TransportII(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ){
    
    
#ifdef PZDEBUG
    if (!flux_mesh || !transport_mesh) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    int mesh_index = 0;
    
    TPZGeoMesh * geometry = flux_mesh->Reference();
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<long,10> dof_indexes;
    int n_shapes;
    
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;
    
    
    TPZManVector<long> indices;
    std::pair<long, long> duplet;
    std::pair<long, std::pair< std::pair<long, long> , std::pair<long, long> > >   cint_ctransport_cmixed_duplet;
    TPZManVector<int,10> face_sides;
    long face_index;
    long n_interfaces;
    
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
        fun_dof_scatter_Gamma.Resize(n_interfaces);
        fun_To_Transport_Gamma.Resize(0, 0);
        fcinterface_ctransport_cmixed_indexes_Gamma.Resize(n_interfaces);
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
        fun_dof_scatter_gamma.Resize(n_interfaces);
        fun_To_Transport_gamma.Resize(0, 0);
        fcinterface_ctransport_cmixed_indexes_gamma.Resize(n_interfaces);
    }
    
    
    // Block size structue Gamma and gamma (inner element interfaces)
    TPZVec< std::pair<long, long> > blocks_dimensions(n_interfaces);
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {
        
        if (IsBoundaryQ) {
            face_index  = finterface_g_indexes_Gamma[k_face];
            duplet      = fleft_right_g_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_g_indexes_gamma[k_face];
            duplet      = fleft_right_g_indexes_gamma[k_face];
        }
        
        face_gel = geometry->Element(face_index);
        
#ifdef PZDEBUG
        if (!face_gel) {
            DebugStop();
        }
#endif
        
        long left_geo_index     = duplet.first;
        long right_geo_index    = duplet.second;
        
        left_gel    = geometry->Element(left_geo_index);
        right_gel   = geometry->Element(right_geo_index);
        
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        // Identify the father element in mixed mesh that touch or contain the left element in transport mesh
        
//        TPZStack<int> smallsides;
        int left_cel_index, father_left_cel_index_index;
        left_cel_index = left_cel->Index();
        TPZCompEl * left_mixed_cel = NULL;
        
        long n_data = fmixed_transport_comp_indexes.size();
        int n_cels;
        for (int i = 0; i < n_data; i++) {
            n_cels = fmixed_transport_comp_indexes[i].second.second.size();
            for (int j = 0; j < n_cels; j++) {
                if (fmixed_transport_comp_indexes[i].second.second[j] == left_cel_index) {
                    father_left_cel_index_index = fmixed_transport_comp_indexes[i].second.first;
                    left_mixed_cel = flux_mesh->Element(father_left_cel_index_index);
                    break;
                }
            }
            if(left_mixed_cel){
                break;
            }
        }
        
        int right_cel_index, father_right_cel_index_index;
        right_cel_index = right_cel->Index();
        TPZCompEl * right_mixed_cel = NULL;
        
        for (int i = 0; i < n_data; i++) {
            n_cels = fmixed_transport_comp_indexes[i].second.second.size();
            for (int j = 0; j < n_cels; j++) {
                if (fmixed_transport_comp_indexes[i].second.second[j] == right_cel_index) {
                    father_right_cel_index_index = fmixed_transport_comp_indexes[i].second.first;
                    right_mixed_cel = flux_mesh->Element(father_right_cel_index_index);
                    break;
                }
            }
            if(right_mixed_cel){
                break;
            }
        }
        
#ifdef PZDEBUG
        
        if (IsBoundaryQ) {
            
            if(!left_mixed_cel){
                DebugStop();
            }
            
        }
        else{
            
            if(!left_mixed_cel || !right_mixed_cel){
                DebugStop();
            }
            
        }

#endif
        
        // End
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // Left based element
        this->ComputeFaceIndex(left_mixed_cel->Reference(),face_sides);
        
        //        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(left_mixed_cel);
//        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
    
        this->ElementDofIndexes(mf_cel, dof_indexes);
        n_shapes = dof_indexes.size();
#ifdef PZDEBUG
        if (dof_indexes.size()==0) {
            DebugStop();
        }
#endif

        
        blocks_dimensions[k_face].first = 1;
        blocks_dimensions[k_face].second = n_shapes;
        
        cint_ctransport_cmixed_duplet.first = face_gel->Reference()->Index();
        cint_ctransport_cmixed_duplet.second.first.first = left_cel->Index();
        cint_ctransport_cmixed_duplet.second.first.second = right_cel->Index();

        
        if (IsBoundaryQ) {
            cint_ctransport_cmixed_duplet.second.second.first = left_mixed_cel->Index();
            cint_ctransport_cmixed_duplet.second.second.second = left_mixed_cel->Index(); // @omar:: left
            fun_dof_scatter_Gamma[k_face] = dof_indexes;
            fcinterface_ctransport_cmixed_indexes_Gamma[k_face] = cint_ctransport_cmixed_duplet;
        }
        else{
            cint_ctransport_cmixed_duplet.second.second.first = left_mixed_cel->Index();
            cint_ctransport_cmixed_duplet.second.second.second = right_mixed_cel->Index();
            fun_dof_scatter_gamma[k_face] = dof_indexes;            
            fcinterface_ctransport_cmixed_indexes_gamma[k_face] = cint_ctransport_cmixed_duplet;
        }
        
    }
    
    // Initialize the matrix
    // Initialize the matrix
    
    if (IsBoundaryQ) {
        fun_To_Transport_Gamma.Initialize(blocks_dimensions);
    }
    else{
        fun_To_Transport_gamma.Initialize(blocks_dimensions);
    }
    
}

/** @brief Initializate diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Fill_un_To_Transport(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ){
    
    
    if (fSimulationData->TransporResolution().first) {
        this->Fill_un_To_TransportII(flux_mesh,transport_mesh,IsBoundaryQ);
        return;
    }
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_un_To_Transport(flux_mesh,transport_mesh,IsBoundaryQ);
    
    int mesh_index = 0;
    TPZGeoMesh * geometry = flux_mesh->Reference();
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<long,10> dof_indexes;
    
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;
    
    TPZManVector<long> indices;
    std::pair<long, long> duplet;
    TPZManVector<int,10> face_sides;
    TPZFMatrix<REAL> normals;
    long face_index;
    long n_interfaces;
    
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
    }
    
    TPZFNMatrix<100,double> block;
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {
        
        if (IsBoundaryQ) {
            face_index  = finterface_g_indexes_Gamma[k_face];
            duplet      = fleft_right_g_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_g_indexes_gamma[k_face];
            duplet      = fleft_right_g_indexes_gamma[k_face];
        }

        face_gel = geometry->Element(face_index);

        
        if (!face_gel) {
            DebugStop();
        }
        
        long left_geo_index     = duplet.first;
        long right_geo_index    = duplet.second;
        
        left_gel    = geometry->Element(left_geo_index);
        right_gel   = geometry->Element(right_geo_index);
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        if(left_gel->HasSubElement() && right_gel->HasSubElement() && face_gel->HasSubElement()){
            DebugStop();
        }
        

        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(left_cel);
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        
        int face_side       = -1;
        int connect_index   = -1;
        
        if(!IdentifyFace(face_side,left_gel,face_gel)){
            std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
            DebugStop();
        }
        else{
            
            for(int ic = 0; ic < intel_vol->NConnects(); ic++){
                connect_index = intel_vol->SideConnectLocId(ic, face_side);
                if (connect_index != -1) {
                    break;
                }
            }
            
            if (connect_index == -1) {
                std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
                DebugStop();
            }
        }
        
        this->ElementDofFaceIndexes(connect_index, mf_cel, dof_indexes);
        TPZIntPoints *int_points   = left_gel->CreateSideIntegrationRule(face_side, left_cel->GetgOrder());
        
        int npoints = int_points->NPoints();
        int nshapes = left_cel->Connect(connect_index).NShape();
  
#ifdef PZDEBUG
        if (IsBoundaryQ) {
            if (1 != fun_To_Transport_Gamma.GetSizeofBlock(k_face).first || nshapes != fun_To_Transport_Gamma.GetSizeofBlock(k_face).second){
                DebugStop();
            }
        }
        else{
            if (1 != fun_To_Transport_gamma.GetSizeofBlock(k_face).first || nshapes != fun_To_Transport_gamma.GetSizeofBlock(k_face).second){
                DebugStop();
            }
            
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(nshapes,1,0.0);
        int el_dim = face_gel->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,nshapes,0.0);
        TPZFNMatrix<50,double> block(npoints,nshapes);
        TPZFNMatrix<50,double> block_integral(1,nshapes,0.0);
        
        REAL w;
        TPZManVector<STATE,2> par_duplet(el_dim,0.0);
        REAL ElementMeasure   = DimensionalMeasure(face_gel);

        REAL detjac;
        TPZFMatrix<REAL> jac,axes,jacinv;

        for (int ip = 0; ip < npoints; ip++) {
         
            // Get the vectorial phi
            int_points->Point(ip, par_duplet, w);
            face_gel->Jacobian(par_duplet, jac, axes, detjac, jacinv);
            intel_vol->SideShapeFunction(face_side, par_duplet, phi, dphidxi);
            
            for (int jp = 0; jp < nshapes; jp++) {
                block_integral(0,jp) +=  w * detjac * phi(jp,0)/ElementMeasure;
            }
            
        }
        
        if (IsBoundaryQ) {
            fun_To_Transport_Gamma.SetBlock(k_face, block_integral);
        }
        else{
            fun_To_Transport_gamma.SetBlock(k_face, block_integral);
        }
        

        
    }
    
    
}


/** @brief Initializate diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Fill_un_To_TransportII(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ){
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_un_To_TransportII(flux_mesh,transport_mesh,IsBoundaryQ);
    
    
    int mesh_index = 0;
    TPZGeoMesh * geometry = flux_mesh->Reference();
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<long,10> dof_indexes;
    
    TPZCompEl * face_cel;
    TPZCompEl * mixed_cel;
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;
    TPZGeoEl * mixed_gel;
    
    int int_order_interfaces = 1;
    
    TPZManVector<long> indices;
    std::pair<long, long> duplet;
    TPZManVector<int,10> face_sides;
    TPZFMatrix<REAL> normals;
    long face_index;
    long n_interfaces;

    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
    }
    
    TPZFNMatrix<100,double> block;
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {
        
        if (IsBoundaryQ) {
            face_index  = finterface_g_indexes_Gamma[k_face];
            duplet      = fleft_right_g_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_g_indexes_gamma[k_face];
            duplet      = fleft_right_g_indexes_gamma[k_face];
        }
        
        face_gel    = geometry->Element(face_index);
        
#ifdef PZDEBUG
        if (!face_gel) {
            DebugStop();
        }
#endif
        
        face_cel = face_gel->Reference();
        
#ifdef PZDEBUG
        if (!face_cel) {
            DebugStop();
        }
#endif
        
        long left_geo_index     = duplet.first;
        long right_geo_index    = duplet.second;
        
        left_gel    = geometry->Element(left_geo_index);
        right_gel   = geometry->Element(right_geo_index);
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        

        
        if (IsBoundaryQ) {
            mixed_cel = flux_mesh->Element(fcinterface_ctransport_cmixed_indexes_Gamma[k_face].second.second.first);
            mixed_gel = mixed_cel->Reference();
        }
        else{
            mixed_cel = flux_mesh->Element(fcinterface_ctransport_cmixed_indexes_gamma[k_face].second.second.first);
            mixed_gel = mixed_cel->Reference();
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        TPZMultiphysicsInterfaceElement * mf_face_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(face_cel);
        
        
        this->ElementDofIndexes(intel_vol, dof_indexes);
        TPZIntPoints *int_points   = face_gel->CreateSideIntegrationRule(face_gel->NSides()-1, int_order_interfaces);
        
        int npoints = int_points->NPoints();
        int nshapes = dof_indexes.size();
        
        // Computing over all integration points of the compuational element mf_cel
        TPZFNMatrix<100,REAL> phi_dot_n(nshapes,1,0.0);
        int face_gel_dim = face_gel->Dimension();
        int intel_vol_dim = intel_vol->Dimension();
        TPZFNMatrix<50,double> block(npoints,nshapes);
        TPZFNMatrix<50,double> block_integral(1,nshapes,0.0);
        
        REAL w;
        TPZManVector<STATE,2> par_duplet(face_gel_dim,0.0);
        TPZManVector<STATE,3> par_mixed_duplet(intel_vol_dim,0.0);
        
        REAL ElementMeasure   = DimensionalMeasure(face_gel);
        
        REAL detjac;
        TPZFMatrix<REAL> jac,axes,jacinv;
        TPZFNMatrix<3,REAL> n;
        TPZVec<int> vectorsides;
        TPZMaterialData data, face_data;
        face_data.fNeedsNormal = true;
        TPZFNMatrix<100,STATE> phi_qs;
        int nphiu,s_i,v_i;
        
        for (int ip = 0; ip < npoints; ip++) {
            
            // Get the vectorial phi
            int_points->Point(ip, par_duplet, w);
            face_gel->Jacobian(par_duplet, jac, axes, detjac, jacinv);

            mf_face_cel->ComputeRequiredData(face_data, par_duplet);
            ComputeTransformation(face_gel, left_gel, mixed_gel, par_duplet, par_mixed_duplet);
            
//            std::cout << "normal = " << face_data.normal <<  std::endl;
            
            intel_vol->InitMaterialData(data);
            intel_vol->ComputeRequiredData(data, par_mixed_duplet);
            
            phi_qs       = data.phi;
            nphiu       = data.fVecShapeIndex.NElements();

            for (int iu = 0; iu < nphiu; iu++)
            {
                
                v_i = data.fVecShapeIndex[iu].first;
                s_i = data.fVecShapeIndex[iu].second;
                
                for (int k = 0; k < intel_vol_dim; k++) {
                        phi_dot_n(iu,0) += 1.0 * phi_qs(s_i,0) * data.fNormalVec(k,v_i) * face_data.normal[k];
                }
                
            }
        
            for (int j = 0; j < nshapes; j++) {
                block_integral(0,j) +=  w * detjac * phi_dot_n(j,0)/ElementMeasure;
            }
            
        }
        
//        std::cout << "face index = " << face_gel->Index() <<  std::endl;
//        std::cout << "mat id = " << face_gel->MaterialId() <<  std::endl;
//        std::cout << "k_face = " << k_face <<  std::endl;
//        std::cout << "dof_indexes = " << dof_indexes <<  std::endl;
//        block_integral.Print(std::cout);
        
        if (IsBoundaryQ) {
            fun_To_Transport_Gamma.SetBlock(k_face, block_integral);
        }
        else{
            fun_To_Transport_gamma.SetBlock(k_face, block_integral);
        }
        
    }
    
    
}

/** @brief Transfer normal fluxes to integration points of transport meshes */
void TRMBuildTransfers::un_To_Transport_Mesh(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_transport, bool IsBoundaryQ){
  
#ifdef PZDEBUG
    if (!cmesh_flux || !cmesh_transport) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
#ifdef PZDEBUG
    if (this->SimulationData()->RawData()->fOmegaIds[0] != 4) {
        DebugStop();
    }
    
#endif
    
    TPZManVector<long,30>  point_index_trans;
    TPZManVector<long,30>  point_index_l;
    TPZManVector<long,30>  point_index_r;
    
    REAL p_avg_n_l = -1.0;
    REAL p_avg_n_r = -1.0;
    
    cmesh_transport->LoadReferences();
    TPZGeoMesh * geometry = cmesh_transport->Reference();
    long n_interfaces;
    int dimension = geometry->Dimension();
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
    }
    
    int rock_id = this->SimulationData()->RawData()->fOmegaIds[0];
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * rock_material = cmesh_flux->FindMaterial(rock_id);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin>  * material_mixe_mem = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(rock_material);

    
    if (IsBoundaryQ) {
        
        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int nbc = this->SimulationData()->RawData()->fGammaIds.size();
        
        for (int ibc = 0; ibc < nbc; ibc++) {
            
            int material_id = this->SimulationData()->RawData()->fGammaIds[ibc];
            //  Getting the total integration point of the destination cmesh
            TPZMaterial * material = cmesh_transport->FindMaterial(material_id);
            
            TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>  * material_bc_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond> *>(material);
            
            if (!material_bc_mem) {
                DebugStop();
            }
            
            // Step one
            TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_Gamma.Cols(),1,0.0);
            long pos = 0;
            for (int iface = 0; iface < n_interfaces; iface++) {
                for(int iflux = 0; iflux < fun_dof_scatter_Gamma[iface].size(); iflux++) {
                    ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_Gamma[iface][iflux],0);
                    pos++;
                }
            }
            
            // Step two
            TPZFMatrix<STATE> un_at_intpoints;
            fun_To_Transport_Gamma.Multiply(ScatterFluxes,un_at_intpoints);
            
            // Step three
            // Trasnfering integrated normal fluxes values
            int counter = 0;
            int i = 0;
            for (int iface = 0; iface < n_interfaces; iface++) {
                
                TPZGeoEl *gel    = geometry->Element(finterface_g_indexes_Gamma[iface]);
                TPZGeoEl *gel_l  = geometry->Element(fleft_right_g_indexes_Gamma[iface].first);
                TPZGeoEl *gel_r  = geometry->Element(fleft_right_g_indexes_Gamma[iface].second);

                
#ifdef PZDEBUG
                if (!gel || !gel_l || !gel_r) {
                    DebugStop();
                }
                
                if (gel_l->Dimension() != dimension) {
                    DebugStop();
                }
#endif
                TPZCompEl * mixed_cel_l = gel_l->Reference();
                TPZCompEl * mixed_cel_r = gel_r->Reference();

                
                if(gel->MaterialId() != material_bc_mem->Id()){
                    counter++;
                    continue;
                }
                
                
#ifdef PZDEBUG
                if (!mixed_cel_l || !mixed_cel_r) {
                    DebugStop();
                }
                
#endif
                
                GlobalPointIndexes(mixed_cel_l, point_index_l);
                
                p_avg_n_l = material_mixe_mem->GetMemory()[point_index_l[0]].p_avg_n();
                material_bc_mem->GetMemory()[i].Set_p_avg_n_l(p_avg_n_l);
                material_bc_mem->GetMemory()[i].Set_un(un_at_intpoints(counter,0));
                i++;
                counter++;
                
            }
        }

        
    }
    else{
        
        
        int material_id = this->SimulationData()->InterfacesMatId();        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = cmesh_transport->FindMaterial(material_id);
        
        TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin>  * material_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin> *>(material);
        
        if (!material_mem) {
            DebugStop();
        }

        int np_cmesh = material_mem->GetMemory().NElements();

        // Step one
        TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_gamma.Cols(),1,0.0);
        long pos = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {
            for(int iflux = 0; iflux < fun_dof_scatter_gamma[iface].size(); iflux++) {
                ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_gamma[iface][iflux],0);
                pos++;
            }
        }

        // Step two
        TPZFMatrix<STATE> un_at_intpoints;
        fun_To_Transport_gamma.Multiply(ScatterFluxes,un_at_intpoints);
        
        // Step three
        // Trasnfering integrated normal fluxes values
        for(long i = 0; i < np_cmesh; i++){
            material_mem->GetMemory()[i].Set_un(un_at_intpoints(i,0));
        }

        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int i = 0;
        int left_mixed_g_index, right_mixed_g_index;        
        for (int iface = 0; iface < n_interfaces; iface++) {

            TPZGeoEl *gel_l  = geometry->Element(fleft_right_g_indexes_gamma[iface].first);
            TPZGeoEl *gel_r  = geometry->Element(fleft_right_g_indexes_gamma[iface].second);
            
#ifdef PZDEBUG
            if (!gel_l || !gel_r) {
                DebugStop();
            }
            
            if (gel_l->Dimension() != dimension || gel_r->Dimension() != dimension) {
                DebugStop();
            }
#endif
            
            TPZCompEl * mixed_cel_l = gel_l->Reference();
            TPZCompEl * mixed_cel_r = gel_r->Reference();
            
            
#ifdef PZDEBUG
            if (!mixed_cel_l || !mixed_cel_r) {
                DebugStop();
            }
            
#endif

            GlobalPointIndexes(mixed_cel_l, point_index_l);
            GlobalPointIndexes(mixed_cel_r, point_index_r);
            
            p_avg_n_l = material_mixe_mem->GetMemory()[point_index_l[0]].p_avg_n();
            p_avg_n_r = material_mixe_mem->GetMemory()[point_index_r[0]].p_avg_n();
            
            material_mem->GetMemory()[i].Set_p_avg_n_l(p_avg_n_l);
            material_mem->GetMemory()[i].Set_p_avg_n_r(p_avg_n_r);
            i++;
            
        }
       
    }
    
    
    return;
    
}


/** @brief Transfer normal fluxes to integration points of transport meshes */
void TRMBuildTransfers::un_To_Transport_MeshII(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_transport, bool IsBoundaryQ){
    
#ifdef PZDEBUG
    if (!cmesh_flux || !cmesh_transport) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
#ifdef PZDEBUG
    if (this->SimulationData()->RawData()->fOmegaIds[0] != 4) {
        DebugStop();
    }
    
#endif
    
    TPZManVector<long,30>  point_index_trans;
    TPZManVector<long,30>  point_index_l;
    TPZManVector<long,30>  point_index_r;
    
    REAL p_avg_n_l = -1.0;
    REAL p_avg_n_r = -1.0;
    
    cmesh_transport->LoadReferences();
    TPZGeoMesh * geometry = cmesh_transport->Reference();
    long n_interfaces;
    int dimension = geometry->Dimension();
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
    }
    
    int rock_id = this->SimulationData()->RawData()->fOmegaIds[0];
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * rock_material = cmesh_flux->FindMaterial(rock_id);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin>  * material_mixe_mem = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(rock_material);
    
    
    if (IsBoundaryQ) {
        
        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int nbc = this->SimulationData()->RawData()->fGammaIds.size();
        
        for (int ibc = 0; ibc < nbc; ibc++) {
            
            int material_id = this->SimulationData()->RawData()->fGammaIds[ibc];
            //  Getting the total integration point of the destination cmesh
            TPZMaterial * material = cmesh_transport->FindMaterial(material_id);
            
            TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>  * material_bc_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond> *>(material);
            
            if (!material_bc_mem) {
                DebugStop();
            }
            
            // Step one
            TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_Gamma.Cols(),1,0.0);
            long pos = 0;
            for (int iface = 0; iface < n_interfaces; iface++) {
                for(int iflux = 0; iflux < fun_dof_scatter_Gamma[iface].size(); iflux++) {
                    ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_Gamma[iface][iflux],0);
                    pos++;
                }
            }
            
            // Step two
            TPZFMatrix<STATE> un_at_intpoints;
            fun_To_Transport_Gamma.Multiply(ScatterFluxes,un_at_intpoints);
            
            // Step three
            // Trasnfering integrated normal fluxes values
            int counter = 0;
            int i = 0;
            int face_g_index, left_mixed_g_index, right_mixed_g_index;
            for (int iface = 0; iface < n_interfaces; iface++) {
                
                //                TPZGeoEl *gel    = geometry->Element(finterface_g_indexes_Gamma[iface]);
                //                TPZGeoEl *gel_l  = geometry->Element(fleft_right_g_indexes_Gamma[iface].first);
                //                TPZGeoEl *gel_r  = geometry->Element(fleft_right_g_indexes_Gamma[iface].second);
                
                face_g_index = fcinterface_ctransport_cmixed_indexes_Gamma[iface].first;
                left_mixed_g_index = fcinterface_ctransport_cmixed_indexes_Gamma[iface].second.second.first;
                right_mixed_g_index = fcinterface_ctransport_cmixed_indexes_Gamma[iface].second.second.second;
                
                TPZGeoEl *gel    = cmesh_transport->Element(face_g_index)->Reference();
                TPZGeoEl *gel_l  = cmesh_flux->Element(left_mixed_g_index)->Reference();
                TPZGeoEl *gel_r  = cmesh_flux->Element(right_mixed_g_index)->Reference();
                
                
#ifdef PZDEBUG
                if (!gel || !gel_l || !gel_r) {
                    DebugStop();
                }
                
                if (gel_l->Dimension() != dimension) {
                    DebugStop();
                }
#endif
                TPZCompEl * mixed_cel_l = gel_l->Reference();
                TPZCompEl * mixed_cel_r = gel_r->Reference();
                
                
                if(gel->MaterialId() != material_bc_mem->Id()){
                    counter++;
                    continue;
                }
                
                
#ifdef PZDEBUG
                if (!mixed_cel_l || !mixed_cel_r) {
                    DebugStop();
                }
                
#endif
                
                GlobalPointIndexes(mixed_cel_l, point_index_l);
                
                p_avg_n_l = material_mixe_mem->GetMemory()[point_index_l[0]].p_avg_n();
                material_bc_mem->GetMemory()[i].Set_p_avg_n_l(p_avg_n_l);
                material_bc_mem->GetMemory()[i].Set_un(un_at_intpoints(counter,0));
                i++;
                counter++;
                
            }
        }
        
        
    }
    else{
        
        
        int material_id = this->SimulationData()->InterfacesMatId();
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = cmesh_transport->FindMaterial(material_id);
        
        TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin>  * material_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin> *>(material);
        
        if (!material_mem) {
            DebugStop();
        }
        
        int np_cmesh = material_mem->GetMemory().NElements();
        
        // Step one
        TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_gamma.Cols(),1,0.0);
        long pos = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {
            for(int iflux = 0; iflux < fun_dof_scatter_gamma[iface].size(); iflux++) {
                ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_gamma[iface][iflux],0);
                pos++;
            }
        }
        
        // Step two
        TPZFMatrix<STATE> un_at_intpoints;
        fun_To_Transport_gamma.Multiply(ScatterFluxes,un_at_intpoints);
        
        // Step three
        // Trasnfering integrated normal fluxes values
        for(long i = 0; i < np_cmesh; i++){
            material_mem->GetMemory()[i].Set_un(un_at_intpoints(i,0));
        }
        
        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int i = 0;
        int left_mixed_g_index, right_mixed_g_index;
        for (int iface = 0; iface < n_interfaces; iface++) {
            
            //            TPZGeoEl *gel_l  = geometry->Element(fleft_right_g_indexes_gamma[iface].first);
            //            TPZGeoEl *gel_r  = geometry->Element(fleft_right_g_indexes_gamma[iface].second);
            
            left_mixed_g_index = fcinterface_ctransport_cmixed_indexes_gamma[iface].second.second.first;
            right_mixed_g_index = fcinterface_ctransport_cmixed_indexes_gamma[iface].second.second.second;
            
            TPZGeoEl *gel_l  = cmesh_flux->Element(left_mixed_g_index)->Reference();
            TPZGeoEl *gel_r  = cmesh_flux->Element(right_mixed_g_index)->Reference();
            
#ifdef PZDEBUG
            if (!gel_l || !gel_r) {
                DebugStop();
            }
            
            if (gel_l->Dimension() != dimension || gel_r->Dimension() != dimension) {
                DebugStop();
            }
#endif
            
            TPZCompEl * mixed_cel_l = gel_l->Reference();
            TPZCompEl * mixed_cel_r = gel_r->Reference();
            
            
#ifdef PZDEBUG
            if (!mixed_cel_l || !mixed_cel_r) {
                DebugStop();
            }
            
#endif
            
            GlobalPointIndexes(mixed_cel_l, point_index_l);
            GlobalPointIndexes(mixed_cel_r, point_index_r);
            
            p_avg_n_l = material_mixe_mem->GetMemory()[point_index_l[0]].p_avg_n();
            p_avg_n_r = material_mixe_mem->GetMemory()[point_index_r[0]].p_avg_n();
            
            material_mem->GetMemory()[i].Set_p_avg_n_l(p_avg_n_l);
            material_mem->GetMemory()[i].Set_p_avg_n_r(p_avg_n_r);
            i++;
            
        }
        
    }
    
    
    return;
    
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Utility Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @brief Get Global integration point indexes associaded  */
void TRMBuildTransfers::GlobalPointIndexes(TPZCompEl * cel, TPZManVector<long,30> &int_point_indexes){
    
    TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
    
#ifdef PZDEBUG
    if(!mf_cel)
    {
        DebugStop();
    }
#endif
    
    mf_cel->GetMemoryIndices(int_point_indexes);
    
}

/** @brief Get Global integration point indexes associaded  */
void TRMBuildTransfers::GlobalPointIndexesInterface(TPZCompEl * int_cel, TPZManVector<long,30> &int_point_indexes){
    
    TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(int_cel);
    
#ifdef PZDEBUG
    if(!mf_int_cel)
    {
        DebugStop();
    }
#endif
    
    mf_int_cel->GetMemoryIndices(int_point_indexes);
    
}

bool TRMBuildTransfers::IdentifyFace(int &side, TPZGeoEl * vol, TPZGeoEl * face){
    
    int volu_nsides = vol->NSides();
    int face_nsides = face->NSides();
    side = -1;
    TPZGeoElSide face_itself =  face->Neighbour(face_nsides-1);
    bool IsMembershipQ = false;
    
    for (int iside = 0; iside < volu_nsides; iside++) {
        IsMembershipQ = bool(vol->NeighbourExists(iside, face_itself));
        if (IsMembershipQ) {
            side = iside;
            break;
        }
    }
    
    TPZGeoElSide vol_itself =  vol->Neighbour(volu_nsides-1);
    
    if(!IsMembershipQ){
        TPZGeoElSide neigh = face->Neighbour(face_nsides-1);
        
        if(!neigh.Element()){
            DebugStop();
        }
        
        TPZGeoElSide neigh_father = neigh.Father2();
        
        if(!neigh_father.Element()){
            DebugStop();
        }
        bool IsNeighQ = false;
        for (int iside = vol_itself.Element()->NNodes(); iside < volu_nsides; iside++) {
 
            IsNeighQ = bool(vol_itself.Element()->NeighbourExists(iside, neigh_father));
            
            if (IsNeighQ) {
                side = iside;
                IsMembershipQ = IsNeighQ;
                break;
            }
        }
    }
    
    return IsMembershipQ;
}

/** @brief Compute parametric transformation form origin to tarfet (xinverse based) */
void TRMBuildTransfers::ComputeTransformation(TPZGeoEl * face_gel_origin, TPZGeoEl * gel_origin , TPZGeoEl * gel_target, TPZVec<REAL> & origin, TPZVec<REAL> & target){
    
#ifdef PZDEBUG
    if (!face_gel_origin || !gel_origin || !gel_target) {
        DebugStop();
    }
    
    if (gel_origin->Dimension() != gel_target->Dimension()) {
        DebugStop();
    }
    
#endif
    
    bool IsmemberQ;
    int face_side;
    TPZManVector<REAL,3> origin_vol(gel_origin->Dimension());
    target.Resize(gel_origin->Dimension(),0.0);
    
    // Transformation step one
    TPZTransform<REAL> tr_face_to_vol(gel_origin->Dimension());
    IsmemberQ =  IdentifyFace(face_side, gel_origin, face_gel_origin);
    
    
#ifdef PZDEBUG
    if (!IsmemberQ) {
        DebugStop();
    }
#endif
    
    tr_face_to_vol = gel_origin->SideToSideTransform(face_side, gel_origin->NSides()-1);
    tr_face_to_vol.Apply(origin, origin_vol);
    gel_origin->TransformSonToFather(gel_target, origin_vol, target);
    
//    std::cout << "origin = "        << origin << std::endl;
//    std::cout << "origin_vol = "    << origin_vol << std::endl;
//    std::cout << "target = "        << target << std::endl;
    
}

/** @brief Compute indices associated to faces on 3D topologies */
void TRMBuildTransfers::ComputeFaceIndex(TPZGeoEl * gel , TPZVec<int> &sides){
    
    
    switch (gel->Type()) {
        case ECube:
        {
            int nfaces = 6;
            sides.Resize(nfaces);
            sides[0] = 20;
            sides[1] = 21;
            sides[2] = 22;
            sides[3] = 23;
            sides[4] = 24;
            sides[5] = 25;
            
        }
            break;
        case ETetraedro:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 10;
            sides[1] = 11;
            sides[2] = 12;
            sides[3] = 13;
            
        }
            break;
        case EQuadrilateral:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 4;
            sides[1] = 5;
            sides[2] = 6;
            sides[3] = 7;
            
        }
            break;
        case ETriangle:
        {
            int nfaces = 3;
            sides.Resize(nfaces);
            sides[0] = 3;
            sides[1] = 4;
            sides[2] = 5;
            
        }
            break;
        default:
        {
            std::cout << "Element not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Compute sides associated to faces on 3D topologies */
void TRMBuildTransfers::ComputeFaceNormals(TPZGeoEl * gel , TPZVec<int> &sides, TPZFMatrix<STATE> &normals){
    
    //  @omar:: Just for linear mapping
    
    TPZFMatrix<REAL> mat_normals;
    TPZVec<int> v_sides;
    gel->ComputeNormals(mat_normals, v_sides);
    
    switch (gel->Type()) {
        case ECube:
        {
            int nfaces = 6;
            sides.Resize(nfaces);
            sides[0] = 20;
            sides[1] = 21;
            sides[2] = 22;
            sides[3] = 23;
            sides[4] = 24;
            sides[5] = 25;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
        break;
        case ETetraedro:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 10;
            sides[1] = 11;
            sides[2] = 12;
            sides[3] = 13;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
        break;
        case EQuadrilateral:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 4;
            sides[1] = 5;
            sides[2] = 6;
            sides[3] = 7;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
        case ETriangle:
        {
            int nfaces = 3;
            sides.Resize(nfaces);
            sides[0] = 3;
            sides[1] = 4;
            sides[2] = 5;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
            
        default:
        {
            std::cout << "Element not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Compute left and right geometric element indexes */
void TRMBuildTransfers::ComputeLeftRight(TPZCompMesh * transport_mesh){
    
    if (fSimulationData->TransporResolution().first) {
        this->ComputeLeftRightII(transport_mesh);
        return;
    }
    
    fleft_right_g_indexes_Gamma.Resize(0);
    fleft_right_g_indexes_gamma.Resize(0);
    
#ifdef PZDEBUG
    if (!transport_mesh) {
        std::cout << "There is no computational transport mesh, transport_mesh = Null." << std::endl;
        DebugStop();
    }
#endif
    
    long nel = transport_mesh->NElements();
    long face_index;
    std::pair <long,long> duplet;
    transport_mesh->LoadReferences();
    int dimension = transport_mesh->Reference()->Dimension();
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = transport_mesh->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsInterfaceElement * interface = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);
        
        if (!interface) {
            continue;
        }
        
        TPZCompEl * left_cel = interface->LeftElement();
        TPZCompEl * right_cel = interface->RightElement();
        
#ifdef PZDEBUG
        
        if(!left_cel || !right_cel){
            DebugStop();
        }
#endif
        
        if(interface->Reference()->HasSubElement()) {
            continue;
        }
        
        face_index  = interface->Reference()->Index();
        duplet      = std::make_pair(left_cel->Reference()->Index(), right_cel->Reference()->Index());
        
        if(left_cel->Dimension() != dimension ||  right_cel->Dimension() != dimension){
            
            fleft_right_g_indexes_Gamma.Push(duplet);
            finterface_g_indexes_Gamma.Push(face_index);
            continue;
        }
        
        fleft_right_g_indexes_gamma.Push(duplet);
        finterface_g_indexes_gamma.Push(face_index);
        
    }
    
//    std::cout << " on Gamma " << std::endl;
//    for (int k = 0; k < fleft_right_g_indexes_Gamma.size(); k++) {
//        std::cout << " volume k : " << k << std::endl;
//        std::cout << " volume left : " << fleft_right_g_indexes_Gamma[k].first << std::endl;
//        std::cout << " volume ritgh : " << fleft_right_g_indexes_Gamma[k].second <<std::endl;
//    }
//    
//    std::cout << " on gamma " << std::endl;
//    for (int k = 0; k < fleft_right_g_indexes_gamma.size(); k++) {
//        std::cout << " volume k : " << k << std::endl;
//        std::cout << " volume left : " << fleft_right_g_indexes_gamma[k].first << std::endl;
//        std::cout << " volume ritgh : " << fleft_right_g_indexes_gamma[k].second <<std::endl;
//    }
    
#ifdef PZDEBUG
    if (finterface_g_indexes_Gamma.size() == 0) {
        DebugStop();
    }
    if (finterface_g_indexes_gamma.size() == 0) {
        std::cout << "Warning:: No inner interfaces were found" << std::endl;
    }
#endif
    
}


/** @brief Compute left and right geometric element indexes */
void TRMBuildTransfers::ComputeLeftRightII(TPZCompMesh * transport_mesh){
    
    
    fleft_right_g_indexes_Gamma.Resize(0);
    fleft_right_g_indexes_gamma.Resize(0);
    
#ifdef PZDEBUG
    if (!transport_mesh) {
        std::cout << "There is no computational transport mesh, transport_mesh = Null." << std::endl;
        DebugStop();
    }
#endif
    
    long nel = transport_mesh->NElements();
    long face_index;
    std::pair <long,long> duplet;
    transport_mesh->LoadReferences();
    int dimension = transport_mesh->Reference()->Dimension();
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = transport_mesh->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsInterfaceElement * interface = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);
        
        if (!interface) {
            continue;
        }
        
        TPZCompEl * left_cel = interface->LeftElement();
        TPZCompEl * right_cel = interface->RightElement();
        
#ifdef PZDEBUG
        
        if(!left_cel || !right_cel){
            DebugStop();
        }
#endif
        
        if(interface->Reference()->HasSubElement()) {
            continue;
        }
        
        face_index  = interface->Reference()->Index();
        duplet      = std::make_pair(left_cel->Reference()->Index(), right_cel->Reference()->Index());
        
        if(left_cel->Dimension() != dimension ||  right_cel->Dimension() != dimension){
            
            fleft_right_g_indexes_Gamma.Push(duplet);
            finterface_g_indexes_Gamma.Push(face_index);
            continue;
        }
        
        fleft_right_g_indexes_gamma.Push(duplet);
        finterface_g_indexes_gamma.Push(face_index);
        
    }
    
    
//    std::cout << " on Gamma " << std::endl;
//    for (int k = 0; k < fleft_right_g_indexes_Gamma.size(); k++) {
//        std::cout << " volume k : " << k << std::endl;
//        std::cout << " volume left : " << fleft_right_g_indexes_Gamma[k].first << std::endl;
//        std::cout << " volume ritgh : " << fleft_right_g_indexes_Gamma[k].second <<std::endl;
//    }
//    
//    std::cout << " on gamma " << std::endl;
//    for (int k = 0; k < fleft_right_g_indexes_gamma.size(); k++) {
//        std::cout << " volume k : " << k << std::endl;
//        std::cout << " volume left : " << fleft_right_g_indexes_gamma[k].first << std::endl;
//        std::cout << " volume ritgh : " << fleft_right_g_indexes_gamma[k].second <<std::endl;
//    }
    
#ifdef PZDEBUG
    if (finterface_g_indexes_Gamma.size() == 0) {
        DebugStop();
    }
    if (finterface_g_indexes_gamma.size() == 0) {
        std::cout << "Warning:: No inner interfaces were found" << std::endl;
    }
#endif
    
}


/** @brief Dimensionla Measure of the elemnt */
REAL TRMBuildTransfers::DimensionalMeasure(TPZGeoEl * gel){
    
#ifdef PZDEBUG
    if (!gel) {
        DebugStop();
    }
#endif
    REAL measure = 0.0;
    int order = 10;
    int element_itself  = gel->NSides() - 1;
    TPZIntPoints * int_points = gel->CreateSideIntegrationRule(element_itself, order);
    REAL detjac, w;
    TPZVec<REAL> par(gel->Dimension(),0.0);
    TPZFMatrix<REAL> jac;
    TPZFMatrix<REAL> axes;
    TPZFMatrix<REAL> jacinv;
    for (int i = 0; i < int_points->NPoints(); i++) {
        int_points->Point(i, par, w);
        gel->Jacobian(par, jac, axes, detjac, jacinv);
        measure += w * detjac;
    }
    
    return measure;
    
}


void TRMBuildTransfers::ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<long> index(0,0);
    int nconnect = intel->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = intel->Connect(icon);
        long seqnumber = con.SequenceNumber();
        long position = intel->Mesh()->Block().Position(seqnumber);
        int nshape = con.NShape();
        for (int ish=0; ish < nshape; ish++) {
            index.Push(position+ ish);
        }
    }
    
    dof_indexes = index;
    return;
}

void TRMBuildTransfers::ElementDofIndexes(TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!m_el) {
        DebugStop();
    }
#endif
    
    TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(m_el->Element(0));
    
#ifdef PZDEBUG
    if (!intel_vol) {
        DebugStop();
    }
#endif
    
//    TPZStack<long> index(0,0);
//    int nconnect = intel_vol->NConnects();
//    for (int icon = 0; icon < nconnect; icon++) {
//        TPZConnect  & con = m_el->Connect(icon);
//        long seqnumber = con.SequenceNumber();
//        long position = m_el->Mesh()->Block().Position(seqnumber);
//        int nshape = con.NShape();
//        for (int ish=0; ish < nshape; ish++) {
//            index.Push(position+ ish);
//        }
//    }
//    
//    dof_indexes = index;
//    return;
    
    
    TPZStack<long> index(0,0);
    int nconnect = intel_vol->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = m_el->Connect(icon);
        long seqnumber = con.SequenceNumber();
        long position = m_el->Mesh()->Block().Position(seqnumber);
        int b_size = m_el->Mesh()->Block().Size(seqnumber);
        for (int ib=0; ib < b_size; ib++) {
            index.Push(position+ ib);
        }
    }
    
    dof_indexes = index;
    return;

}

void TRMBuildTransfers::ElementDofFaceIndexes(int connect_index,TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<long> index(0,0);
    TPZConnect  & con = intel->Connect(connect_index);
    long seqnumber = con.SequenceNumber();
    long position = intel->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}

void TRMBuildTransfers::ElementDofIndexes(TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes, int el_index){
    
#ifdef PZDEBUG
    if (!m_el) {
        DebugStop();
    }
#endif
    
    TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(m_el->Element(el_index));
    
#ifdef PZDEBUG
    if (!intel_vol) {
        DebugStop();
    }
#endif
    
    int start = 0;
    int end = intel_vol->NConnects();
    if (el_index == 1) {
        
        TPZInterpolationSpace * intel_vol_q = dynamic_cast<TPZInterpolationSpace * >(m_el->Element(0));
        
#ifdef PZDEBUG
        if (!intel_vol_q) {
            DebugStop();
        }
#endif
        start = intel_vol_q->NConnects();
        end = m_el->NConnects();
    }
    
    TPZStack<long> index(0,0);
    int nconnect = end;
    for (int icon = start; icon < nconnect; icon++) {
        TPZConnect  & con = m_el->Connect(icon);
        long seqnumber = con.SequenceNumber();
        long position = m_el->Mesh()->Block().Position(seqnumber);
        int b_size = m_el->Mesh()->Block().Size(seqnumber);
        for (int ib=0; ib < b_size; ib++) {
            index.Push(position+ ib);
        }
    }
    
    dof_indexes = index;
    return;    
}


void TRMBuildTransfers::ElementDofFaceIndexes(int connect_index, TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!m_el && connect_index > 4) {
        DebugStop();
    }
#endif
    

    TPZStack<long> index(0,0);
    TPZConnect  & con = m_el->Connect(connect_index);
    long seqnumber = con.SequenceNumber();
    long position = m_el->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}


/** @brief Compute compuational mesh pair (mixed, transport) indexed by geometric volumetic element index */
void TRMBuildTransfers::FillComputationalElPairs(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_transport){
    
    
    if (fSimulationData->TransporResolution().first) {
        this->FillComputationalElPairsII(cmesh_mf_mixed,cmesh_mf_transport);
        return;
    }

    fmixed_transport_cindexes.Resize(0);
    
#ifdef PZDEBUG
    if (!cmesh_mf_mixed) {
        DebugStop();
    }
    
    if (!fSimulationData->IsOnePhaseQ() && !cmesh_mf_transport) {
        DebugStop();
    }
    
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    int dimension = geometry->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<long, std::pair <long,long> > gel_indexes;
    
    for (long i = 0; i < geometry->NElements(); i++) {
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->Dimension() != dimension || gel->HasSubElement() || gel->MaterialId() > 7) {
            continue;
        }
        gel_indexes.first = gel->Index();
        gel_indexes.second.first = -1;
        gel_indexes.second.second = -1;
        fmixed_transport_cindexes.Push(gel_indexes);

    }
    
    // counting volumetric elements
    long nvol_elements = fmixed_transport_cindexes.size();
    fmixed_transport_cindexes.Resize(nvol_elements);
    
    // inserting mixed indexes
    cmesh_mf_mixed->LoadReferences();
    for(long ivol = 0; ivol < nvol_elements; ivol++){

        TPZCompEl * mixed_cel = geometry->Element(fmixed_transport_cindexes[ivol].first)->Reference();
        
#ifdef PZDEBUG
        if (!mixed_cel) {
            //continue;
            DebugStop();
        }
#endif
        fmixed_transport_cindexes[ivol].second.first = mixed_cel->Index();
        
    }
    
    if(fSimulationData->IsOnePhaseQ()){
        
//        for (int k = 0; k < fmixed_transport_cindexes.size(); k++) {
//            std::cout << " volume k : " << k <<std::endl;
//            std::cout << " volume gel : " << fmixed_transport_cindexes[k].first <<std::endl;
//            std::cout << " volume cmixed : " << fmixed_transport_cindexes[k].second.first <<std::endl;
//            std::cout << " volume ctransport : " << fmixed_transport_cindexes[k].second.second <<std::endl;
//        }
        
        return;
    }
    
    // inserting transport indexes
    cmesh_mf_transport->LoadReferences();
    for(long ivol = 0; ivol < nvol_elements; ivol++){
        
        TPZCompEl * trans_cel = geometry->Element(fmixed_transport_cindexes[ivol].first)->Reference();
        
#ifdef PZDEBUG
        if (!trans_cel) {
            DebugStop();
        }
#endif
        fmixed_transport_cindexes[ivol].second.second = trans_cel->Index();
        
    }
    
//    for (int k = 0; k < fmixed_transport_cindexes.size(); k++) {
//        std::cout << " volume k : " << k <<std::endl;
//        std::cout << " volume gel : " << fmixed_transport_cindexes[k].first <<std::endl;
//        std::cout << " volume cmixed : " << fmixed_transport_cindexes[k].second.first <<std::endl;
//        std::cout << " volume ctransport : " << fmixed_transport_cindexes[k].second.second <<std::endl;
//    }
    
}


/** @brief Compute compuational mesh pair (mixed, transport) indexed by geometric volumetic element index */
void TRMBuildTransfers::FillComputationalElPairsII(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_transport){
    
    fmixed_transport_comp_indexes.Resize(0);
    
#ifdef PZDEBUG
    if (!cmesh_mf_mixed) {
        DebugStop();
    }
    
    if (!fSimulationData->IsOnePhaseQ() && !cmesh_mf_transport) {
        DebugStop();
    }
    
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    int dimension = geometry->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<long, std::pair <long, std::vector<long> > > gel_cel_indexes;
    
    for (long icel = 0; icel < cmesh_mf_mixed->NElements(); icel++) {
        
        TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(icel);
        
#ifdef PZDEBUG
        if (!mixed_cel) {
            DebugStop();
        }
#endif
        TPZGeoEl * gel = mixed_cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->Dimension() != dimension) {
            continue;
        }
        gel_cel_indexes.first = gel->Index();
        gel_cel_indexes.second.first = mixed_cel->Index();
        gel_cel_indexes.second.second.resize(0);
        fmixed_transport_comp_indexes.Push(gel_cel_indexes);
        
    }
    
    // counting volumetric elements
    int nvol_elements = fmixed_transport_comp_indexes.size();
    
    if(fSimulationData->IsOnePhaseQ()){
        
//        for (int k = 0; k < fmixed_transport_comp_indexes.size(); k++) {
//            std::cout << " volume k : " << k <<std::endl;
//            std::cout << " volume gel : " << fmixed_transport_comp_indexes[k].first <<std::endl;
//            std::cout << " volume cmixed : " << fmixed_transport_comp_indexes[k].second.first <<std::endl;
//        }
        
        return;
    }
    
    // inserting transport indexes
    cmesh_mf_transport->LoadReferences();
    TPZVec<TPZGeoEl *> n_refined_sons;
    int cel_index;
    for(long ivol = 0; ivol < nvol_elements; ivol++){
        
        TPZGeoEl * father_gel = geometry->Element(fmixed_transport_comp_indexes[ivol].first);
        
#ifdef PZDEBUG
        if (!father_gel) {
            DebugStop();
        }
#endif
        n_refined_sons.resize(0);
        father_gel->GetHigherSubElements(n_refined_sons);
        
        if (n_refined_sons.size() == 0) {
            n_refined_sons.resize(1);
            n_refined_sons[0] = father_gel;
        }
        
        for (int igel = 0; igel < n_refined_sons.size(); igel++) {
            cel_index = geometry->Element(n_refined_sons[igel]->Index())->Reference()->Index();
            fmixed_transport_comp_indexes[ivol].second.second.push_back(cel_index);
        }
        
    }
    
//    for (int k = 0; k < fmixed_transport_comp_indexes.size(); k++) {
//        std::cout << " volume k : " << k <<std::endl;
//        std::cout << " volume gel : " << fmixed_transport_comp_indexes[k].first <<std::endl;
//        std::cout << " volume cmixed : " << fmixed_transport_comp_indexes[k].second.first <<std::endl;
//        for (int igel = 0; igel < n_refined_sons.size(); igel++) {
//            int index = fmixed_transport_comp_indexes[k].second.second[igel]; ;
//            std::cout << " volume ctransport : " << index <<std::endl;
//        }
//    }
    
}