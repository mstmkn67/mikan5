#include "KelvinGreenFunc.h"

#define GPI (0.0625/3.14159285358979)
#define TPI (0.1250/3.14159265358979)

namespace kelvin_green_func{

Tensor3x3 delta(1.0,0.0,0.0,
                0.0,1.0,0.0,
                0.0,0.0,1.0);


Tensor3x3 get_free_G(const Vector3d& r,const Vector3d& rp,double p){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();
	return GPI/(1.-p)*((3.-4.*p)*delta+dyad(dr,dr)*l2)*sqrt(l2);
}


Tensor3x3x3 get_free_T(const Vector3d& r,const Vector3d& rp,double p){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();double l=sqrt(l2);
	return -TPI/(1.-p)*l*l2*(
	         (1.-2.*p)*(dyad(delta,dr)+dyad(dr,delta)-diamond(dr,delta))
	        +3.*l2*dyad(dr,dyad(dr,dr)));

}

}

#undef TPI
#undef GPI
