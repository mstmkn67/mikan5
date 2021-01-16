#ifndef _STOKES_GREEN_FUNC2D_H_
#define _STOKES_GREEN_FUNC2D_H_

#include "../Tensor3x3x3.h"

namespace stokes_green_func2d{

Tensor3x3 get_free_G(const Vector3d& r,const Vector3d& rp);
Vector3d get_free_P(const Vector3d& r,const Vector3d& rp);
Tensor3x3x3 get_free_T(const Vector3d& r,const Vector3d& rp);
//
Tensor3x3x3 get_grad_free_G(const Vector3d& r,const Vector3d& rp);

}

#endif // _STOKES_GREEN_FUNC2D_H_
