#include "BoundaryCondition.h"

Vector3d BoundaryCondition::difPosition(const Vector3d& r0,const Vector3d& r1){
	return r1-r0;
}

void BoundaryCondition::difPositionAndVelocity(const Vector3d& r0,const Vector3d& r1,
	const Vector3d& v0,const Vector3d& v1,Vector3d& r01,Vector3d&v01){
	r01=r1-r0;
	v01=v1-v0;
}

Vector3d BoundaryCondition::getPosition(const Vector3d& r){
	return r;
}
//

SimplePeriodic::SimplePeriodic(const Vector3d& mS,const Vector3d& s)
:min(mS),size(s){}

Vector3d SimplePeriodic::difPosition(const Vector3d& r0,const Vector3d& r1){
	Vector3d r=r1-r0;
	return Vector3d(r.x-floor(r.x/size.x+0.5)*size.x,
	                r.y-floor(r.y/size.y+0.5)*size.y,
                  r.z-floor(r.z/size.z+0.5)*size.z);
}

void SimplePeriodic::difPositionAndVelocity(const Vector3d& r0,const Vector3d& r1,
	const Vector3d& v0,const Vector3d& v1,Vector3d& r01,Vector3d&v01){
	Vector3d r=r1-r0;
	r01.set(r.x-floor(r.x/size.x+0.5)*size.x,
	        r.y-floor(r.y/size.y+0.5)*size.y,
          r.z-floor(r.z/size.z+0.5)*size.z);
	v01=v1-v0;
}


Vector3d SimplePeriodic::getPosition(const Vector3d& r){
	return Vector3d(r.x-floor((r.x-min.x)/size.x)*size.x, 
		              r.y-floor((r.y-min.y)/size.y)*size.y,
		              r.z-floor((r.z-min.z)/size.z)*size.z );
}

//
LeesEdwards::LeesEdwards(const Vector3d& minSize,const Vector3d& s_size,const double t)
	:size(s_size),min(minSize),dt(t),dx(0)
{
	vmax=size.y;
}

LeesEdwards::~LeesEdwards()
{
}

void LeesEdwards::update()
{
	dx=vmax*dt*(*time);//イメージセルとのずれdx
	dx-=(int)floor(dx/size.x)*size.x;
}

Vector3d LeesEdwards::difPosition(const Vector3d& r0,const Vector3d& r1)
{
	Vector3d r=r1-r0;
	double cory=floor(r.y/size.y+0.5);
	double x=r.x-dx*cory;
	return Vector3d(x    -floor(x/size.x+0.5)*size.x,
	                r.y-cory*size.y,
                  r.z-floor(r.z/size.z+0.5)*size.z);
}

void LeesEdwards::difPositionAndVelocity(const Vector3d& r0,const Vector3d& r1,
	                                       const Vector3d& v0,const Vector3d& v1,
	                                       Vector3d& r01,Vector3d& v01)
{
	Vector3d r=r1-r0;
	double cory=floor(r.y/size.y+0.5);
	double x=r.x-dx*cory;
	r01 = Vector3d(x    -floor(x/size.x+0.5)*size.x,
	               r.y-cory*size.y,
                 r.z-floor(r.z/size.z+0.5)*size.z);
	v01=v1-v0;
	v01.x-=vmax*cory;
}

Vector3d LeesEdwards::getPosition(const Vector3d& coord)
{	
	Vector3d r=coord;
	double cory=floor((coord.y - min.y)/size.y);
	r.x-=dx*cory;
	return Vector3d(r.x - floor((coord.x-min.x)/size.x)*size.x, 
					coord.y - cory*size.y,
					coord.z - floor((coord.z-min.z)/size.z)*size.z );
}
