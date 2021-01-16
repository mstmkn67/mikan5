//EllipticRotationField.h

#ifndef _ELLIPTIC_ROTATION_FIELD_H_
#define _ELLIPTOC_ROTATION_FIELD_H_

#include "ElectricMagneticField.h"
using namespace std;

class EllipticRotationField:public ElectricMagneticField{
public:
	EllipticRotationField(double angular_velocity_x,const Vector3d& Bx,
                                double angular_velocity_y,const Vector3d& By,double dt);
	virtual ~EllipticRotationField();
	virtual void initial();
	virtual void update();
	virtual Vector3d getField(const Vector3d& position);
private:
	double ax,ay;
	double dt;
	Vector3d Bx;
	Vector3d By;
};

#endif // _ELLIPTIC_ROTATION_FIELD_H_
