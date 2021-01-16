//tensor�̓����

#ifndef _TENSOR_UTILITY_H_
#define _TENSOR_UTILITY_H_

#include "../Tensor3x3.h"
#include "../Tensor3x3x3.h"
#include "../Tensor3x3x3x3.h"
#include "../solver/LapackFunctions.h"

#define EPS 1e-10

namespace tensor_utility{

inline double delta(int i,int j){
	if(i==j)return 1.0;
	else return 0.0;
}

inline Tensor3x3 delta(){
	return Tensor3x3(
		1.0,0.0,0.0,
		0.0,1.0,0.0,
		0.0,0.0,1.0
	);
}

inline double epsilon(int i,int j,int k){
	if(i==j || j==k || k==i)return 0;
	if(i==0)
		if(j==1)return 1.;
		else return -1.;
	if(i==1)
		if(j==0)return -1.;
		else return 1.;
	if(j==1)return -1.;
	else return 1.;
}

inline Tensor3x3x3 epsilon(){
	return Tensor3x3x3(
		0.0,0.0,0.0,
		0.0,0.0,1.0,
		0.0,-1.0,0.0,//
		0.0,0.0,-1.0,
		0.0,0.0,0.0,
		1.0,0.0,0.0,//
		0.0,1.0,0.0,
		-1.0,0.0,0.0,
		0.0,0.0,0.0//
 	);
}

//3rd-rank tensor convert from normal tensor to tilde tensor
inline Tensor3x3x3
Gijk2Gkji(const Tensor3x3x3& G){
	Tensor3x3x3 H;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				H(i,j,k)=G.get(k,j,i);
	return H;
}

//convert from hydro-sub-tensor ABC from abc (vice versas)
inline void
convertABC2abc(const Tensor3x3& A,const Tensor3x3& B,const Tensor3x3& C,
			   Tensor3x3& a,Tensor3x3& b,Tensor3x3& c){
	double w[6];
	double x[]={
		A.x.x,A.x.y,A.x.z,B.x.x,B.y.x,B.z.x,
		A.y.x,A.y.y,A.y.z,B.x.y,B.y.y,B.z.y,
		A.z.x,A.z.y,A.z.z,B.x.z,B.y.z,B.z.z,
		B.x.x,B.x.y,B.x.z,C.x.x,C.x.y,C.x.z,
		B.y.x,B.y.y,B.y.z,C.y.x,C.y.y,C.y.z,
		B.z.x,B.z.y,B.z.z,C.z.x,C.z.y,C.z.z
	};
	//lap::get_eigenvalues_symmetric(6,x,w);
	lap::get_eigenvalues_and_vectors_symmetric(6,x,w);
 	if(fabs(w[0])<EPS || fabs(w[1])<EPS || fabs(w[2])<EPS || fabs(w[3])<EPS || fabs(w[4])<EPS || fabs(w[5])<EPS){
		//if eigenvalues are zero, the corresponding eigen vectors are neglected.
		double y[]={0.0,0.0,0.0, 0.0,0.0,0.0,
		            0.0,0.0,0.0, 0.0,0.0,0.0,
		            0.0,0.0,0.0, 0.0,0.0,0.0,
		            0.0,0.0,0.0, 0.0,0.0,0.0,
		            0.0,0.0,0.0, 0.0,0.0,0.0,
		            0.0,0.0,0.0, 0.0,0.0,0.0};
		for(int aa=0;aa<6;aa++){
			if(fabs(w[aa])>EPS){
				for(int i=0;i<6;i++){
					for(int j=0;j<6;j++){
						y[6*i+j]+=x[6*aa+i]*x[6*aa+j]/w[aa];
					}
				}
			}
		}
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				a(i,j)=y[6*i+j];
				b(i,j)=y[6*i+j+3];
				c(i,j)=y[6*(i+3)+j+3];
			}
		}
	}else{
		Tensor3x3 Bt=B.transpose();
		Tensor3x3 C_=C.invert();
		a=(A-Bt*C_*B).invert();
		b=(-C_*B*a).transpose();
		c=(C-B*A.invert()*Bt).invert();
	}
}

inline void 
convertResistance2Mobility(const Tensor3x3& A,const Tensor3x3& B,const Tensor3x3& C,
                           const Tensor3x3x3& G,const Tensor3x3x3& H,const Tensor3x3x3x3& K,
                           Tensor3x3& a,Tensor3x3& b,Tensor3x3& c,
                           Tensor3x3x3& g,Tensor3x3x3& h,Tensor3x3x3x3& k){
	convertABC2abc(A,B,C,a,b,c);
	g=G*a+H*b;
	h=G*b.transpose()+H*c;
	k=-K+g*Gijk2Gkji(G)+h*Gijk2Gkji(H);
}

//translation theorem for resistance tensor
//below functions express the traslation of hydrodynamic center(not object)
//Tensor3x3 translation_A(); invariant for translation
inline Tensor3x3 
translation_B(const Tensor3x3& A0,const Tensor3x3& B0,const Vector3d& r01)
{
	return B0-(r01^A0);
}

inline Tensor3x3 
translation_C(const Tensor3x3& A0,const Tensor3x3& B0,const Tensor3x3& C0,
              const Vector3d& r01)
{
	return C0-(r01^(A0^r01))+(B0^r01)-(r01^(B0.transpose()));
}

inline Tensor3x3x3
translation_G(const Tensor3x3& A0,const Tensor3x3x3& G0,const Vector3d& r01)
{
	return G0-0.5*(dyad(r01,A0)+diamond(A0,r01));
}

inline Tensor3x3x3 
translation_H(const Tensor3x3& A0,const Tensor3x3& B0,const Tensor3x3x3& G0,
              const Tensor3x3x3& H0,const Vector3d& r01)
{
	Tensor3x3 t=(A0^r01)+B0.transpose();
	return H0+(G0^r01)-0.5*(diamond(r01,t)+dyad(r01,t));
}

inline Tensor3x3x3x3 
translation_K(const Tensor3x3& A0,const Tensor3x3x3& G0,const Tensor3x3x3x3& K0,
              const Vector3d& r01)
{
	Tensor3x3x3 G0_=Gijk2Gkji(G0);
	return K0-0.5*(dyad(G0,r01)+diamond(G0,r01)+dyad(r01,G0_)+diamond(r01,G0_))
	     +0.25*(dyad(r01,dyad(A0,r01))+dyad(r01,diamond(r01,A0))
	           +diamond(r01,dyad(A0,r01))+diamond(diamond(r01,A0),r01));
}

//translation theorem for mobility tensor
inline Tensor3x3
translation_a(const Tensor3x3& a0,const Tensor3x3& b0,const Tensor3x3& c0,
              const Vector3d& r01)
{
	return a0+(b0.transpose()^r01)-(r01^b0)-(r01^(c0^r01));
}

inline Tensor3x3 
translation_b(const Tensor3x3& b0,const Tensor3x3& c0,const Vector3d& r01)
{
	return b0+(c0^r01);
}

//void translation_c(); invariant for translation
inline Tensor3x3x3 
translation_g(const Tensor3x3x3& g0,const Tensor3x3x3& h0,const Vector3d& r01)
{
	return g0+(h0^r01)-0.5*(dyad(r01,delta())+diamond(r01,delta()));
}

//void translation_h(); invariant for translation
//void translation_k(); invariant for translation

//b and c are the sub mobility tensors
//return values is relative position from particle's center
inline Vector3d 
calcCenterOfMobility(const Tensor3x3& b,const Tensor3x3& c){
	Tensor3x3 t=b-b.transpose();
	if(fabs(c.det())<EPS)return Vector3d(0.0,0.0,0.0);
	return 0.5*(c.trace()*delta()-c).invert()*dot2(epsilon(),t);
}

//A and B are the sub resistance tensors
//return values is relative position from particle's center
inline Vector3d
calcCenterOfResistance(const Tensor3x3& A,const Tensor3x3& B){
	Tensor3x3 t=B-B.transpose();
	if(fabs(A.det())<EPS)return Vector3d(0.0,0.0,0.0);
	return 0.5*(A.trace()*delta()-A).invert()*dot2(epsilon(),t);
}


}//end of namespace tensor_utility

#endif

