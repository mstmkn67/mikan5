#ifndef _DYNAMIC_FLOW_H_
#define _DYNAMIC_FLOW_H_

#include "VelocityField.h"

class DynamicShearFlow:public VelocityField{
public:
    DynamicShearFlow(double gamma,double omega);
    virtual ~DynamicShearFlow();
	virtual void initial();
	virtual void update();
	virtual Vector3d getVelocity(const Vector3d& coord);
	virtual Vector3d getAngularVelocity(const Vector3d& coord);
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);
private:
    double gw;
    double w;
};

#endif // _DYNAMIC_FLOW_H_
