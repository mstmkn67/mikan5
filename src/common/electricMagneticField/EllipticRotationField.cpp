#include "EllipticRotationField.h"

EllipticRotationField::EllipticRotationField(
double x,const Vector3d& _Bx,double y,const Vector3d& _By,double t)
:ax(x),Bx(_Bx),ay(y),By(_By),dt(t){
}

EllipticRotationField::~EllipticRotationField(){}

void EllipticRotationField::initial(){}

void EllipticRotationField::update(){}

Vector3d EllipticRotationField::getField(const Vector3d& position){
	double tt=*time;
	return Bx*cos(tt*ax)+By*sin(tt*ay);
}
