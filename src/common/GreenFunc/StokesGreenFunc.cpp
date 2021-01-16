#include "StokesGreenFunc.h"

#define GPI (0.125/3.14159285358979)

namespace stokes_green_func{

Tensor3x3 delta(1.0,0.0,0.0,
                0.0,1.0,0.0,
                0.0,0.0,1.0);

/////////////////////////////////////////////////////////////
Tensor3x3 getG_F(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();
	return (delta+dyad(dr,dr)*l2)*sqrt(l2);
}

Tensor3x3 getG_D(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();
	double l=sqrt(l2);double l3=l2*l;double l5=l3*l2;
	Tensor3x3 a(l3-3.*dr.x*dr.x*l5,  -3.*dr.x*dr.y*l5,    3.*dr.x*dr.z*l5,
	              -3.*dr.y*dr.x*l5,l3-3.*dr.y*dr.y*l5,    3.*dr.y*dr.z*l5,
	              -3.*dr.z*dr.x*l5,  -3.*dr.z*dr.y*l5,-l3+3.*dr.z*dr.z*l5);
	return a;
}

Tensor3x3 getG_Q(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();double l=sqrt(l2);double l3=l2*l;
	Tensor3x3 a(0.0,     0.0,    -dr.x*l3,
	            0.0,     0.0,    -dr.y*l3,
	            -dr.x*l3,-dr.y*l3,0.0);
	return dr.z*getG_D(r,rp)+a;
}
//
Vector3d getP_F(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();double l=sqrt(l2);
	return 2.*dr*l2*l;
}

Vector3d getP_Q(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();double l=sqrt(l2);double l5=l2*l2*l;
	return Vector3d(-6.*dr.z*dr.x*l5,
	                -6.*dr.z*dr.y*l5,
	                -2.*(dr.x*dr.x+dr.y*dr.y-2.*dr.z*dr.z)*l5);	
}
//
Tensor3x3x3 getT_F(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();double l=sqrt(l2);double l5=l2*l2*l;
	return -6.*dyad(dr,dyad(dr,dr))*l5;
}

Tensor3x3x3 getT_D(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();double l=sqrt(l2);double l5=l2*l2*l;double l7=l5*l2;
	return Tensor3x3x3(
		-18.*dr.x*l5+30.*dr.x*dr.x*dr.x*l7, -6.*dr.y*l5+30.*dr.x*dr.x*dr.y*l7,-6.*dr.z*l5+30.*dr.x*dr.x*dr.z*l7,
		 -6.*dr.y*l5+30.*dr.x*dr.x*dr.y*l7, -6.*dr.x*l5+30.*dr.x*dr.y*dr.y*l7,            30.*dr.x*dr.y*dr.z*l7,
		  6.*dr.z*l5-30.*dr.x*dr.x*dr.z*l7,            -30.*dr.x*dr.y*dr.z*l7, 6.*dr.x*l5-30.*dr.x*dr.z*dr.z*l7,//
		
		-6.*dr.y*l5+30.*dr.y*dr.x*dr.x*l7,  -6.*dr.x*l5+30.*dr.y*dr.x*dr.y*l7,             30.*dr.y*dr.x*dr.z*l7,
		-6.*dr.x*l5+30.*dr.y*dr.y*dr.x*l7, -18.*dr.y*l5+30.*dr.y*dr.y*dr.y*l7, -6.*dr.z*l5+30.*dr.y*dr.y*dr.z*l7,
		           -30.*dr.y*dr.z*dr.x*l7,   6.*dr.z*l5-30.*dr.y*dr.z*dr.y*l7,  6.*dr.y*l5-30.*dr.y*dr.z*dr.z*l7,//
		
		-6.*dr.z*l5+30.*dr.z*dr.x*dr.x*l7,            30.*dr.z*dr.x*dr.y*l7, -6.*dr.x*l5+30.*dr.z*dr.x*dr.z*l7,
		            30.*dr.z*dr.y*dr.x*l7,-6.*dr.z*l5+30.*dr.z*dr.y*dr.y*l7, -6.*dr.y*l5+30.*dr.z*dr.y*dr.z*l7,
		 6.*dr.x*l5-30.*dr.z*dr.z*dr.x*l7, 6.*dr.y*l5-30.*dr.z*dr.z*dr.y*l7, 18.*dr.z*l5-30.*dr.z*dr.z*dr.z*l7//
	);
}

Tensor3x3x3 getT_Q(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();double l=sqrt(l2);double l5=6.*l2*l2*l;
	Tensor3x3x3 a(
		dr.x*dr.z*l5,             0.0,          0.0,
		dr.y*dr.z*l5,             0.0,          0.0,
		(dr.x*dr.x-dr.z*dr.z)*l5, dr.x*dr.y*l5, dr.x*dr.z*l5,//
		
		0.0,          dr.x*dr.z*l5,            0.0,
		0.0,          dr.y*dr.z*l5,            0.0,
		dr.x*dr.y*l5, (dr.y*dr.y-dr.z*dr.z)*l5, dr.y*dr.z*l5,//
		
		0.0,          0.0,          dr.x*dr.z*l5,
		0.0,          0.0,          dr.y*dr.z*l5,
		dr.x*dr.z*l5, dr.y*dr.z*l5, 0.0
	);
	return dr.z*getT_D(r,rp)+a;
}

//
//////////////////////////////

Tensor3x3 get_free_G(const Vector3d& r,const Vector3d& rp){
	return GPI*getG_F(r,rp);
}

Tensor3x3 get_wall_G(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi(rp.x,rp.y,-rp.z);
	Tensor3x3 G=getG_F(r,rp);
	G-=getG_F(r,rpi);
	G+=2.*rp.z*rp.z*getG_D(r,rpi);
	G-=2.*rp.z*getG_Q(r,rpi);
	return GPI*G;
}

Vector3d get_free_P(const Vector3d& r,const Vector3d& rp){
	return GPI*getP_F(r,rp);
}

Vector3d get_wall_P(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi(rp.x,rp.y,-rp.z);
	Vector3d P=getP_F(r,rp);
	P-=getP_F(r,rpi);
	P-=2.*rp.z*getP_Q(r,rpi);
	return GPI*P;
}

Tensor3x3x3 get_free_T(const Vector3d& r,const Vector3d& rp){
	return GPI*getT_F(r,rp);
}

Tensor3x3x3 get_wall_T(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi(rp.x,rp.y,-rp.z);
	Tensor3x3x3 T=getT_F(r,rp);
	T-=getT_F(r,rpi);
	T+=2.*rp.z*rp.z*getT_D(r,rpi);
	T-=2.*rp.z*getT_Q(r,rpi);
	return GPI*T;
}

///////////////////////////////////////////////////////

Tensor3x3x3 get_grad_free_G(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();double l=sqrt(l2);
	Tensor3x3x3 t=-diamond(delta,dr)+dyad(delta,dr)+dyad(dr,delta)-3.*dyad(dr,dyad(dr,dr))*l2;
	return  GPI*t*l2*l;
}


}

#undef GPI
