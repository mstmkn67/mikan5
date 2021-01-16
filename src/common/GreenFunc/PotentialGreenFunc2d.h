#ifndef _POTENTIAL_GREEN_FUNC_2D_H_
#define _POTENTIAL_GREEN_FUNC_2D_H_

#include "../Vector2d.h"

namespace potential_green_func2d{

double get_free_G(const Vector2d& r,const Vector2d& rp);
Vector2d get_free_T(const Vector2d& r,const Vector2d& rp);

//The wall is located on z=0.
//double get_zero_potential_on_wall_G(const Vector2d& r,const Vector2d& rp);
//Vector2d get_zero_potential_on_wall_T(const Vector2d& r,const Vector2d& rp);

//double get_no_flux_on_wall_G(const Vector2d& r,const Vector2d& rp);
//Vector2d get_no_flux_on_wall_T(const Vector2d& r,const Vector2d& rp);

}


#endif // _POTENTIAL_GREEN_FUNC__2D_H_
