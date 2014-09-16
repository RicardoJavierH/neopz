/**
 * @file
 * @brief Contains the implementation of the TPZShapeLinear methods.
 */

#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzerror.h"
#include "pzreal.h"
using namespace std;

namespace pzshape {
	
	
	void TPZShapeLinear::Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
		// Quadratic or higher shape functions
		if(num <= 0) return;
		phi.Put(0,0,1.0);
		dphi.Put(0,0, 0.0);
		if(num == 1) return;
		phi.Put(1,0, x);
		dphi.Put(0,1, 1.0);
		int ord;
		for(ord = 2;ord<num;ord++) {
			phi.Put(ord,0, 2.0*x*phi(ord-1,0) - phi(ord-2,0));
			dphi.Put(0,ord, 2.0*x*dphi(0,ord-1) + 2.0*phi(ord-1,0) - dphi(0,ord-2));
		}
	}
	
	void TPZShapeLinear::Legendre(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
		
		// Quadratic or higher shape functions
		if (num <= 0) return;
		phi.Put(0, 0, 1.0);
		dphi.Put(0, 0, 0.0);
		
		if (num == 1) return;
		
		phi.Put(1, 0, x);
		dphi.Put(0, 1, 1.0);
		
		int ord;
		//Aqui fica diferente do Chebyshev
		REAL ord_real, value;
		for (ord = 2; ord < num; ord++)
		{
			//casting int ord to REAL ord_real
			ord_real = (REAL)ord;
			//computing the ord_th function
			value    = ( (2.0 * (ord_real - 1.0) + 1.0) * x * phi(ord - 1, 0) - (ord_real - 1.0) * phi(ord - 2 , 0) ) / (ord_real);
			phi.Put(ord, 0, value);
			
			//computing the ord_th function's derivative
			value    = (2.0 * (ord_real - 1.0) + 1.0) * phi(ord - 1, 0) + dphi(0, ord - 2);
			dphi.Put(0, ord, value);
		}
		
#ifdef DEBUG
		int printing = 0;
		if (printing){
			cout << "Legendre" << endl;
			for(ord = 0; ord < num; ord++)
			{
				cout << "x = " << x << endl;
				cout << "phi(" << ord << ", 0) = " << phi(ord, 0) << endl;
				cout << "dphi(0, " << ord << " = " << dphi(0, ord) << endl;
				cout << endl;
			}
		}
#endif
		
	} //end of method
	
	void TPZShapeLinear::Legendre(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, int nderiv){
		
		// Quadratic or higher shape functions
		if (num <= 0) return;
		phi.Put(0, 0, 1.0);
		dphi.Put(0, 0, 0.0);
		
		int ideriv;
		for (ideriv = 1; ideriv < nderiv; ideriv++)
			dphi.Put(ideriv, 0, 0.0);
		
		
		if (num == 1) return;
		
		phi.Put(1, 0, x);
		dphi.Put(0, 1, 1.0);
		
		for (ideriv = 1; ideriv < nderiv; ideriv++)
			dphi.Put(ideriv, 1, 0.0);
		
		int ord;
		//Aqui fica diferente do Chebyshev
		REAL ord_real, value;
		for (ord = 2; ord < num; ord++)
		{
			//casting int ord to REAL ord_real
			ord_real = (REAL)ord;
			//computing the ord_th function
			value    = ( (2.0 * (ord_real - 1.0) + 1.0) * x * phi(ord - 1, 0) - (ord_real - 1.0) * phi(ord - 2 , 0) ) / (ord_real);
			phi.Put(ord, 0, value);
			
			//computing the ord_th function's derivative
			value    = (2.0 * (ord_real - 1.0) + 1.0) * phi(ord - 1, 0) + dphi(0, ord - 2);
			dphi.Put(0, ord, value);
			
			for (ideriv = 1; ideriv < nderiv; ideriv++){
				value = (2.0 * (ord_real - 1.0) + 1.0) * dphi(ideriv - 1, ord - 1) + dphi(ideriv, ord - 2);
				dphi.Put(ideriv, ord, value);	    	 
			}
			
		}
		
#ifdef DEBUG
		int printing = 0;
		if (printing){
			cout << "Legendre" << endl;
			for(ord = 0; ord < num; ord++)
			{
				cout << "x = " << x << endl;
				cout << "phi(" << ord << ", 0) = " << phi(ord, 0) << endl;
				cout << "dphi(0, " << ord << " = " << dphi(0, ord) << endl;
				cout << endl;
			}
		}
#endif
		
	} //end of method
	
	void TPZShapeLinear::Jacobi(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	}
	
	// Setting Chebyshev polynomials as orthogonal sequence generating shape functions
	void (*TPZShapeLinear::fOrthogonal)(REAL, int, TPZFMatrix<REAL> &, TPZFMatrix<REAL> &) = TPZShapeLinear::Chebyshev;
	
	/**
	 * Computes the generating shape functions for a quadrilateral element
	 * @param pt (input) point where the shape function is computed
	 * @param phi (input/output) value of the  shape functions
	 * @param dphi (input/output) value of the derivatives of the shape functions holding the derivatives in a column
	 */
	void TPZShapeLinear::ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
	{
		
		phi(2,0) = phi(0,0)*phi(1,0);
		dphi(0,2) = dphi(0,0)*phi(1,0)+phi(0,0)*dphi(0,1);
		
		phi(2,0) *= 4.;
		dphi(0,2) *= 4.;
		
	}
	
	void TPZShapeLinear::Shape(TPZVec<REAL> &x,TPZVec<long> &id, TPZVec<int> &order,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		//	num = number of functions to compute
#ifndef NODEBUG
		if ( order[0] < 0 ) {
			PZError << "Compelbas::shape --> Invalid dimension for arguments: order = " << order[0]
			<< " phi.Rows = " << (int) phi.Rows() << " dphi.Cols = " << (int) dphi.Cols() << "\n";
			return;
		}
		if(phi.Rows() < order[0]+1) {
			PZError << "TPZShapeLinear::shape --> Invalid dimension for argument phi " << endl;
			phi.Resize(order[0]+1, phi.Cols());
		}
		if(dphi.Cols() < order[0]+1) {
			PZError << "TPZShapeLinear::shape --> Invalid dimension for argument dphi " << endl;
			dphi.Resize(dphi.Rows(),order[0]+1);
		}
#endif
		
		if ( order[0] == 0) 
		{
			phi(0,0) = 1.;
			dphi(0,0) = 0.;
		} else 
		{		// Linear shape functions
			phi(0,0) = (1-x[0])/2.;
			phi(1,0) = (1+x[0])/2.;
			dphi(0,0) = -0.5;
			dphi(0,1)= 0.5;
		}
		
		int is,d;
		TPZFNMatrix<100> phiblend(NSides,1),dphiblend(Dimension,NSides);
		for(is=0; is<NCornerNodes; is++)
		{
			phiblend(is,0) = phi(is,0);
			for(d=0; d<Dimension; d++)
			{
				dphiblend(d,is) = dphi(d,is);
			}
		}
		ShapeGenerating(x,phiblend,dphiblend);
		// Quadratic or higher shape functions
		int num2 = order[0]-1;
		int transformationindex = GetTransformId1d(id);
		TPZFNMatrix<10> phiint(num2,1),dphiint(1,num2);
		if(num2 > 0)
		{
			ShapeInternal(x,order[0],phiint,dphiint,transformationindex);
		}
		int ord;
		for (ord = 2; ord < order[0]+1; ord++) {    // even functions
			dphi(0,ord) = dphiint(0,ord-2)*phiblend(2,0)+dphiblend(0,2)*phiint(ord-2,0);
			phi(ord,0) = phiint(ord-2,0)*phiblend(2,0);
		}
	}
	
	void TPZShapeLinear::SideShape(int side, TPZVec<REAL> &pt, TPZVec<long> &id, TPZVec<int> &order,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		switch(side) {
			case 0:
			case 1:
				TPZShapePoint::Shape(pt,id,order,phi,dphi);
				break;
			case 2:
				Shape(pt,id,order,phi,dphi);
				break;
		}
	}
	
    void TPZShapeLinear::ShapeOrder(TPZVec<long> &id, TPZVec<int> &order, TPZGenMatrix<int> &shapeorders, TPZVec<long> &sides)
    {
        DebugStop();
    }
    
    
    void TPZShapeLinear::SideShapeOrder(int side,  TPZVec<long> &id, int order, TPZGenMatrix<int> &shapeorders)
    {
        DebugStop();
    }
    
	void TPZShapeLinear::ShapeInternal(TPZVec<REAL>  &x, int ord,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,int transformation_index){
		// Quadratic or higher shape functions
		int num = ord-1;
		if(num <= 0) return;
		REAL y;
		TransformPoint1d(transformation_index,x[0],y);
		fOrthogonal(y,num,phi,dphi);
		TransformDerivative1d(transformation_index,num,dphi);
	}
	
	void TPZShapeLinear::TransformPoint1d(int transid,REAL in,REAL &out) {
		if (!transid) out =  in;
		else          out = -in;
	}
	
	void TPZShapeLinear::TransformDerivative1d(int transid,int num,TPZFMatrix<REAL> &in) {
		
		if(transid == 0) return;
		int i;
		for(i=0;i<num;i++) in(0,i) = -in(0,i);
	}
	
	int TPZShapeLinear::GetTransformId1d(TPZVec<long> &id) {
		if (id[1] < id[0]) return 1;
		else               return 0;
	}
	
	int TPZShapeLinear::NConnectShapeF(int side, int order) {
		if(side<2) return 1;//0 a 4
		if(side<3) return (order-1);//6 a 14
		PZError << "TPZShapeLinear::NConnectShapeF, bad parameter side " << side << endl;
		return 0;
	}
	
	int TPZShapeLinear::NShapeF(TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}
	
#ifdef _AUTODIFF
	void TPZShapeLinear::ShapeInternal(FADREAL & x,int num,TPZVec<FADREAL> & phi,int transformation_index){
		// Quadratic or higher shape functions
		if(num <= 0) return;
		FADREAL y;
		TransformPoint1d(transformation_index,x,y);
		FADfOrthogonal(y,num,phi);
		//  TransformDerivative1d(transformation_index,num,phi);
	}
	
	void TPZShapeLinear::TransformPoint1d(int transid,FADREAL & in,FADREAL &out) {
		if (!transid) out =  in;
		else          out = -in;
	}
	
	void TPZShapeLinear::Chebyshev(FADREAL & x,int num,TPZVec<FADREAL> &phi){
		// Quadratic or higher shape functions
		if(num <= 0) return;
		//phi.Put(0,0,1.0);
		//dphi.Put(0,0, 0.0);
		phi[0] = 1.0; // <!> Remark: the derivatives other than the 0th are set to null
		if(num == 1) return;
		//phi.Put(1,0, x);
		//dphi.Put(0,1, 1.0);
		phi[1] = x;
		//phi[1].fastAccessDx(0)=1.0; // <!> Remark: the derivatives other than the 0th aren't set to null
		//just ensuring the derivatives are properly initialized and that FAD objects of more than
		// one derivative are used
		int ord;
		for(ord = 2;ord<num;ord++) {
			//phi.Put(ord,0, 2.0*x*phi(ord-1,0) - phi(ord-2,0));
			//dphi.Put(0,ord, 2.0*x*dphi(0,ord-1) + 2.0*phi(ord-1,0) - dphi(0,ord-2));
			phi[ord] = x * phi[ord-1] * 2.0 - phi[ord-2];
		}
	}
	
	void (*TPZShapeLinear::FADfOrthogonal)(FADREAL&,int ,TPZVec<FADREAL> &) =  TPZShapeLinear::Chebyshev/*(FADREAL&, int, TPZVec<FADREAL>&)*/;//Chebyshev;

#endif
	
};
