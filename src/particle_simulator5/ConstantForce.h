//ConstantForce

#ifndef _CONSTANT_FORCE_H_
#define _CONSTANT_FORCE_H_

#include "Force.h"

class ConstantForce:public Force{
public:
	ConstantForce(const vector<Particle*>& particles,const Vector3d& force,const Vector3d& torque,
		const Vector3d& actPoint,bool frameFlag);//true --> particleFrame, false-->simulationFrame
	virtual ~ConstantForce();
	virtual void update();
private:
	bool beadFlag,particleFrameFlag;
	Vector3d actPoint,force,torque;
};

#endif // _CONSTANT_FORCE_H_
