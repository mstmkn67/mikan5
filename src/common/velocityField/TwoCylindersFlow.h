//TwoCylindersFlow.h
//二つの円筒内における流れ場

#ifndef _TWO_CYLINDERS_FLOW_H_
#define _TWO_CYLINDERS_FLOW_H_

#include "../Vector3d.h"
#include "../Tensor3x3.h"
using namespace std;

class TwoCylindersFlow:public VelocityField{
public:
	TwoCylindersFlow(double r_out,double r_in,double omega){
		A=omega/(1.0-r_in*r_in/r_out/r_out);
		B=r_in*r_in;
	};
	virtual ~TwoCylindersFlow(){};
	virtual void initial(){};
	virtual void update(){};
	virtual Vector3d getVelocity(const Vector3d& coord);//場の速度
	virtual Vector3d getAngularVelocity(const Vector3d& coord);//場の角速度
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);//場の変形速度テンソル
protected:
private:
	double A,B;
};

inline
Vector3d TwoCylindersFlow::getVelocity(const Vector3d& coord){
	double r2=coord.x*coord.x+coord.y*coord.y;
	double c=A*(1.-B/r2);
	return Vector3d(c*coord.y,-c*coord.x,0.0);
}
inline
Vector3d TwoCylindersFlow::getAngularVelocity(const Vector3d& coord){
	return Vector3d(0.0,0.0,-A);
}
inline
Tensor3x3 TwoCylindersFlow::getRateOfStrainTensor(const Vector3d& coord){
	double r2=coord.x*coord.x+coord.y*coord.y;
	double dvxdx= 2.*A*B*coord.x*coord.y/r2/r2;
	double dvxdy=    A*B*(coord.y*coord.y-coord.x*coord.x)/r2/r2;
	double dvydy=-2.*A*B*coord.x*coord.y/r2/r2;
	return Tensor3x3(dvxdx,dvxdy,0.0,
	                 dvxdy,dvydy,0.0,
	                 0.0,  0.0,  0.0);
}
#endif // _TWO_CYLINDERS_FLOW_H_
