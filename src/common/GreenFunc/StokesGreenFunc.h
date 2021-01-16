#ifndef _STOKES_GREEN_FUNC_H_
#define _STOKES_GREEN_FUNC_H_

#include "../Tensor3x3x3.h"

namespace stokes_green_func{

Tensor3x3 get_free_G(const Vector3d& r,const Vector3d& rp);
Tensor3x3 get_wall_G(const Vector3d& r,const Vector3d& rp);
Vector3d get_free_P(const Vector3d& r,const Vector3d& rp);
Vector3d get_wall_P(const Vector3d& r,const Vector3d& rp);
Tensor3x3x3 get_free_T(const Vector3d& r,const Vector3d& rp);
Tensor3x3x3 get_wall_T(const Vector3d& r,const Vector3d& rp);


Tensor3x3x3 get_grad_free_G(const Vector3d& r,const Vector3d& rp);

//////////////////////////////

Tensor3x3 getG_F(const Vector3d& r,const Vector3d& rp);
Tensor3x3 getG_D(const Vector3d& r,const Vector3d& rp);
Tensor3x3 getG_Q(const Vector3d& r,const Vector3d& rp);
Vector3d getP_F(const Vector3d& r,const Vector3d& rp);
Vector3d getP_Q(const Vector3d& r,const Vector3d& rp);
Tensor3x3x3 getT_F(const Vector3d& r,const Vector3d& rp);
Tensor3x3x3 getT_D(const Vector3d& r,const Vector3d& rp);
Tensor3x3x3 getT_Q(const Vector3d& r,const Vector3d& rp);
}

#endif // _STOKES_GREEN_FUNC_H_
