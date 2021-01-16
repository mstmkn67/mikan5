#include "DynamicFlow.h"

 DynamicShearFlow::DynamicShearFlow(double gd,double _w):gw(gd*w),w(_w){

 }

DynamicShearFlow::~DynamicShearFlow(){

}

void DynamicShearFlow::initial(){

}

void DynamicShearFlow::update(){

}

Vector3d DynamicShearFlow::getVelocity(const Vector3d& coord){
    return Vector3d(gw*cos(w*(*time))*coord.y,0,0);
}

Vector3d DynamicShearFlow::getAngularVelocity(const Vector3d& coord){
    return Vector3d(0,0,-0.5*gw*cos(w*(*time)));
}

Tensor3x3 DynamicShearFlow::getRateOfStrainTensor(const Vector3d& coord){
    double c=0.5*gw*cos(w*(*time));
    return Tensor3x3(0,c,0, c,0,0, 0,0,0);
}