#ifndef _BEADS_INTERACTION_H_
#define _BEADS_INTERACTION_H_

#include "../Tensor3x3.h"

namespace beads_interaction{


//RotnePragerYamakawaTensor
//Jens Rotne and Stephen Prager
//
//Journal of Chemical Physics 
//Volume 50, 4831(1969)
//Hiromi Yamakawa
//Journal of Chemical Physics 
//Volume 53, 436(1970)
//
// length is reduced by radius of bead
// return value is reduced by 1/(vicosity x length unit)
Tensor3x3 getRPYTensor(const Vector3d& ri,const Vector3d& rj);

//
//
//
Tensor3x3 getRPYTensor(double ai,const Vector3d& ri,double aj,const Vector3d& rj);

}

#endif // _BEADS_INTERACTION_H_
