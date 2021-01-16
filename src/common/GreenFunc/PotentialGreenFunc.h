#ifndef _POTENTIAL_GREEN_FUNC_H_
#define _POTENTIAL_GREEN_FUNC_H_

#include "../Vector3d.h"

namespace potential_green_func{

double get_free_G(const Vector3d& r,const Vector3d& rp);
Vector3d get_free_T(const Vector3d& r,const Vector3d& rp);

//The wall is located on z=0.
double get_zero_potential_on_wall_G(const Vector3d& r,const Vector3d& rp);
Vector3d get_zero_potential_on_wall_T(const Vector3d& r,const Vector3d& rp);

double get_no_flux_on_wall_G(const Vector3d& r,const Vector3d& rp);
Vector3d get_no_flux_on_wall_T(const Vector3d& r,const Vector3d& rp);

}


#endif // _POTENTIAL_GREEN_FUNC_H_
