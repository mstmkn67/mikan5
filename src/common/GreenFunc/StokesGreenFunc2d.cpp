#include "StokesGreenFunc.h"

#define GPI (0.25/3.14159285358979)
#define PPI (0.5/3.14159285358979)
#define TPI (-1.0/3.14159285358979)

namespace stokes_green_func2d{

Tensor3x3 delta(1.0,0.0,0.0,
                0.0,1.0,0.0,
                0.0,0.0,1.0);

Tensor3x3 get_free_G(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=dr.length2();
	return GPI*(-log(sqrt(l2))*delta+dyad(dr,dr)/l2);
}

Vector3d get_free_P(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=dr.length2();
	return PPI*dr/l2;
}

Tensor3x3x3 get_free_T(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=dr.length2();
	return TPI*dyad(dr,dyad(dr,dr))/l2/l2;
}

//
Tensor3x3x3 get_grad_free_G(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();
	Tensor3x3x3 t=-diamond(delta,dr)+dyad(delta,dr)+dyad(dr,delta)-2.*dyad(dr,dyad(dr,dr))*l2;
	return GPI*l2*t;
}

}

#undef GPI
#undef PPI
#undef TPI
