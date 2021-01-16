#ifndef _KELVIN_GREEN_FUNC2D_H_
#define _KELVIN_GREEN_FUNC2D_H_

#include "../Tensor3x3x3.h"

namespace kelvin_green_func2d{

Tensor3x3 get_free_G(const Vector3d& r,const Vector3d& rp,double poisson_ratio);
Tensor3x3x3 get_free_T(const Vector3d& r,const Vector3d& rp,double poisson_ratio);

}

#endif // _KELVIN_GREEN_FUNC2D_H_
