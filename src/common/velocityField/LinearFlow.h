//LinearFlow.h
//üŒ`—¬‚ê
//ŽÀsŽž‚É—¬‘Ì‚ð‰ð‚©‚È‚¢

#ifndef _LINEAR_FLOW_H_
#define _LINEAR_FLOW_H_

#include "VelocityField.h"

/////////////////NonFlow/////////////////////////////////
class NonFlow:public VelocityField{
public:
	NonFlow(){};
	virtual ~NonFlow(){};
	virtual void initial(){};
	virtual void update(){};
	virtual Vector3d getVelocity(const Vector3d& coord);
	virtual Vector3d getAngularVelocity(const Vector3d& coord);
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);
};
inline
Vector3d NonFlow::getVelocity(const Vector3d& coord){
	return Vector3d(0.0,0.0,0.0);
}
inline
Vector3d NonFlow::getAngularVelocity(const Vector3d& coord){
	return Vector3d(0.0,0.0,0.0);
}
inline
Tensor3x3 NonFlow::getRateOfStrainTensor(const Vector3d& coord){
	return Tensor3x3(0.0,0.0,0.0,
	                 0.0,0.0,0.0,
	                 0.0,0.0,0.0);
}
/////////////////SimpleShearFlow/////////////////////////////////
class SimpleShearFlow:public VelocityField{
public:
	SimpleShearFlow(double gd):gammaDot(gd){};
	virtual ~SimpleShearFlow(){};
	virtual void initial(){};
	virtual void update(){};
	virtual Vector3d getVelocity(const Vector3d& coord);
	virtual Vector3d getAngularVelocity(const Vector3d& coord);
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);
private:
	double gammaDot;
};
inline
Vector3d SimpleShearFlow::getVelocity(const Vector3d& coord){
	return Vector3d(gammaDot*coord.y,0.0,0.0);
}
inline
Vector3d SimpleShearFlow::getAngularVelocity(const Vector3d& coord){
	return Vector3d(0.0,0.0,-0.5*gammaDot);
}
inline
Tensor3x3 SimpleShearFlow::getRateOfStrainTensor(const Vector3d& coord){
	return Tensor3x3(0.0,         0.5*gammaDot,0.0,
	                 0.5*gammaDot,0.0,         0.0,
	                 0.0,         0.0,         0.0);
}
/////////////////ElongationalFlow////////////////////////////////////////////
//ˆêŽ²L’·(k=0,epsilonDot>0)
//‚QŽ²L’·(k=0,epsilonDot<0)
//•½–ÊL’·(k=1)
class ElongationalFlow:public VelocityField{
public:
	ElongationalFlow(double _k,double e):k(_k),epsilonDot(e){};
	virtual ~ElongationalFlow(){};
	virtual void initial(){};
	virtual void update(){};
	virtual Vector3d getVelocity(const Vector3d& coord);
	virtual Vector3d getAngularVelocity(const Vector3d& coord);
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);
private:
	double epsilonDot;//strain rate
	double k;//parameter 0<=k<=1
};
inline
Vector3d ElongationalFlow::getVelocity(const Vector3d& coord){
	return Vector3d(epsilonDot*coord.x,-0.5*epsilonDot*(1+k)*coord.y,-0.5*epsilonDot*(1-k)*coord.z);
}
inline
Vector3d ElongationalFlow::getAngularVelocity(const Vector3d& coord){
	return Vector3d(0.0,0.0,0.0);
}
inline
Tensor3x3 ElongationalFlow::getRateOfStrainTensor(const Vector3d& coord){
	return Tensor3x3(epsilonDot,0.0,                  0.0,
	                 0.0,       -0.5*epsilonDot*(1+k),0.0,
	                 0.0,       0.0,                  -0.5*epsilonDot*(1-k));
}
/////////////////LinearFlow(generalize)/////////////////////////////
//v_i=a_i+b_{ij}r_j
class LinearFlow:public VelocityField
{
public:
	LinearFlow(const Vector3d& _a,const Tensor3x3& _b):a(_a),b(_b){};
	virtual ~LinearFlow(){};
	virtual void initial(){};
	virtual void update(){};
	virtual Vector3d getVelocity(const Vector3d& coord);
	virtual Vector3d getAngularVelocity(const Vector3d& coord);
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);
private:
	const Vector3d a;const Tensor3x3 b;
};
inline 
Vector3d LinearFlow::getVelocity(const Vector3d& coord)
{
	return a+b*coord;
}
inline 
Vector3d LinearFlow::getAngularVelocity(const Vector3d& coord)
{
	return Vector3d(b.z.y-b.y.z,b.x.z-b.z.x,b.y.x-b.x.y);
}
inline
Tensor3x3 LinearFlow::getRateOfStrainTensor(const Vector3d& coord)
{
	return Tensor3x3(b.x.x,            0.5*(b.y.x+b.x.y),0.5*(b.z.x+b.x.z),
	                 0.5*(b.y.x+b.x.y),b.y.y,            0.5*(b.y.z+b.z.y),
	                 0.5*(b.z.x+b.x.z),0.5*(b.y.z+b.z.y),b.z.z            );
}

#endif // _LINEAR_FLOW_H_
