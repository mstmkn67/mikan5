#include "PotentialGreenFunc.h"

#define GPI (0.25/3.14159285358979)

namespace potential_green_func{

double get_free_G(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l=1./dr.length();
	return GPI*l;
}

Vector3d get_free_T(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();
	return -GPI*l2*sqrt(l2)*dr;
}


double get_zero_potential_on_wall_G(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi=Vector3d(rp.x,rp.y,-rp.z);
	return get_free_G(r,rp)-get_free_G(r,rpi);
}

Vector3d get_zero_potential_on_wall_T(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi=Vector3d(rp.x,rp.y,-rp.z);
	return get_free_T(r,rp)-get_free_T(r,rpi);
}

double get_no_flux_on_wall_G(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi=Vector3d(rp.x,rp.y,-rp.z);
	return get_free_G(r,rp)+get_free_G(r,rpi);
}

Vector3d get_no_flux_on_wall_T(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi=Vector3d(rp.x,rp.y,-rp.z);
	return get_free_T(r,rp)+get_free_T(r,rpi);
}

}

