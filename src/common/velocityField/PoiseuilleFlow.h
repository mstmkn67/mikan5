#ifndef _POISEUILLE_FLOW_H_
#define _POISEUILLE_FLOW_H_

#include "VelocityField.h"

class PlanePoiseuilleFlow:public VelocityField{
public:
	PlanePoiseuilleFlow(double h,double v):hm2(1.0/(h*h)),vmax(v){};
	virtual ~PlanePoiseuilleFlow(){};
	virtual void initial(){};
	virtual void update(){};
	virtual Vector3d getVelocity(const Vector3d& coord){
		double y2=coord.y*coord.y;
		return Vector3d(vmax*(1.-y2*hm2),0.0,0.0);
	};
	virtual Vector3d getAngularVelocity(const Vector3d& coord){
		return Vector3d(0.0,0.0,vmax*coord.y*hm2);
	};
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord){
		double a=-vmax*coord.y*hm2;
		return Tensor3x3(0.0,a,0.0, a,0.0,0.0, 0.0,0.0,0.0);
	};
private:
	double vmax;
	double hm2;
};

class CylinderPoiseuilleFlow:public VelocityField{
public:
	CylinderPoiseuilleFlow(double radius,double v):Rm2(1./(radius*radius)),vmax(v){};
	virtual ~CylinderPoiseuilleFlow(){};
	virtual void initial(){};
	virtual void update(){};
	virtual Vector3d getVelocity(const Vector3d& coord){
		double r2=coord.y*coord.y+coord.z*coord.z;
		return Vector3d(vmax*(1.-r2*Rm2),0.0,0.0);
	};
	virtual Vector3d getAngularVelocity(const Vector3d& coord){
		return Vector3d(0.0,-vmax*coord.z*Rm2,vmax*coord.y*Rm2);
	};
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord){
		double a=-vmax*Rm2;
		return Tensor3x3(0.0,a*coord.y,a*coord.z,
		                 a*coord.y,0.0,0.0,
		                 a*coord.z,0.0,0.0);
	};
private:
	double vmax;
	double Rm2;
};

#endif // _POISEUILLE_FLOW_H_
