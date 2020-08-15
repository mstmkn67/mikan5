#include "ConstantForce.h"

ConstantForce::ConstantForce(
	const vector<Particle*>& p,const Vector3d& f,const Vector3d& t,const Vector3d& a,bool pf)
	:Force(p),force(f),torque(t),actPoint(a),particleFrameFlag(pf){}

ConstantForce::~ConstantForce(){}

void ConstantForce::update()
{
	vector<Particle*>::iterator p=particles.begin();
	for(;p!=particles.end();p++){
		EulerAngle& e=(*p)->eulerAngle;
		Vector3d rp=e.rotate_(actPoint);
		if(particleFrameFlag){
			Vector3d ft=e.rotate_(force);
			(*p)->force+=ft;
			(*p)->torque+=e.rotate_(torque)+(rp^ft);
		}else{
			(*p)->force+=force;
			(*p)->torque+=(torque+(rp^force));
		}
	}
}
