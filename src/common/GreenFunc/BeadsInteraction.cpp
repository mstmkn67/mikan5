#include "BeadsInteraction.h"

#define ZETA1 0.053051647697298445256294587790838
#define ZETA2 0.039788735772973833942220940843129

namespace beads_interaction{

Tensor3x3 delta(1.0,0.0,0.0,
                0.0,1.0,0.0,
                0.0,0.0,1.0);

Tensor3x3 getRPYTensor(const Vector3d& ri,const Vector3d& rj){
	Vector3d p=rj-ri;
	double l2=p.length2();
	double l=sqrt(l2);
	if(l>2.0){
		p/=l;
		return ZETA2/l*((1.+2./3./l2)*delta+(1.-2./l2)*dyad(p,p));
	}else{
		return ZETA1*((1.-9.*l/32.)*delta+3./32.*dyad(p,p)/l);
	}
}


Tensor3x3 getRPYTensor(double ai,const Vector3d& ri,double aj,const Vector3d& rj){
	Vector3d p=rj-ri;
	double l2=p.length2();
	double l=sqrt(l2);
	return ZETA2/l*(delta+dyad(p,p)/l2
	             +(ai*ai+aj*aj)/l2*(delta/3.0-dyad(p,p)/l2));
}


}

#undef ZETA1
#undef ZETA2
