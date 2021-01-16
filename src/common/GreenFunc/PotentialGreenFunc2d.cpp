#include "PotentialGreenFunc2d.h"

#define GPI (0.5/3.14159285358979)

namespace potential_green_func2d{

double get_free_G(const Vector2d& r,const Vector2d& rp){
	Vector2d dr=r-rp;
	double l=dr.length();
	return -GPI*log(l);
}

Vector2d get_free_T(const Vector2d& r,const Vector2d& rp){
	Vector2d dr=r-rp;
	double l2=1./dr.length2();
	return -GPI*l2*dr;
}


//double get_zero_potential_on_wall_G(const Vector2d& r,const Vector2d& rp){
//	Vector2d rpi=Vector3d(rp.x,rp.y,-rp.z);
//	return get_free_G(r,rp)-get_free_G(r,rpi);
//}
//
//Vector2d get_zero_potential_on_wall_T(const Vector2d& r,const Vector2d& rp){
//	Vector2d rpi=Vector2d(rp.x,rp.y,-rp.z);
//	return get_free_T(r,rp)-get_free_T(r,rpi);
//}
//
//double get_no_flux_on_wall_G(const Vector2d& r,const Vector2d& rp){
//	Vector3d rpi=Vector2d(rp.x,rp.y,-rp.z);
//	return get_free_G(r,rp)+get_free_G(r,rpi);
//}
//
//Vector2d get_no_flux_on_wall_T(const Vector2d& r,const Vector2d& rp){
//	Vector2d rpi=Vector2d(rp.x,rp.y,-rp.z);
//	return get_free_T(r,rp)+get_free_T(r,rpi);
//}

}

#undef GPI
