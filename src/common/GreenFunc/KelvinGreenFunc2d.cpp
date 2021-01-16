#include "StokesGreenFunc.h"

#define GPI (0.125/3.14159285358979)
#define TPI (0.25/3.14159285358979)

namespace kelvin_green_func2d{

Tensor3x3 delta(1.0,0.0,0.0,
                0.0,1.0,0.0,
                0.0,0.0,1.0);

Tensor3x3 get_free_G(const Vector3d& r,const Vector3d& rp,double p){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();
	return GPI/(1.-p)*((3.-4.*p)*log(sqrt(l2))*delta+dyad(dr,dr)*l2);
}

Tensor3x3x3 get_free_T(const Vector3d& r,const Vector3d& rp,double p){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();double l=sqrt(l2);
	return -TPI/(1.-p)*l2*((1.-2.*p)*(diamond(delta,dr)+dyad(dr,delta)-dyad(delta,dr))*l
	                       +2.*l2*dyad(dr,dyad(dr,dr)));
}

}

#undef GPI
#undef TPI
