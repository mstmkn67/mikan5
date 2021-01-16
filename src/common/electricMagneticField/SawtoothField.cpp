#include "SawtoothField.h"

SawtoothField::SawtoothField(double p,double a,double l,double st,double d)
:L(l),tau(st),dt(d){
	field1.x=-p/(a*L);
	field2.x= p/((1.0-a)*L);
	aL=a*L;
}

SawtoothField::~SawtoothField(){}

void SawtoothField::initial(){}

void SawtoothField::update(){}

Vector3d SawtoothField::getField(const Vector3d& r){
	double tt=(*time)*dt;
	if(int(tt/tau)%2==1){
		return Vector3d(0.0,0.0,0.0);
	}
	double x=r.x-L*floor(r.x/L);
	if(x<aL){
		return field1;
	}
	return field2;
}
