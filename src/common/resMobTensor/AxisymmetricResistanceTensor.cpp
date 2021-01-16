#include "AxisymmetricResistanceTensor.h"

////////////////////////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265358979
//////////////////////////////
//BaseTensorOFAxissymmetric struct//
//////////////////////////////
BaseTensorOfAxisymmetry::BaseTensorOfAxisymmetry(const Vector3d& _u):u(_u){
	double l=u.length();
	u/=l;
	update();
}
void BaseTensorOfAxisymmetry::update(){
	uu0=dyad(u,u);//uu0
	uu1=tensor_utility::delta()-uu0;//uu1
	uuu=diamond(tensor_utility::epsilon()*u,u)+dyad(u,tensor_utility::epsilon()*u);//uuu
	uuuu0=1.5*dyad((uu0-tensor_utility::delta()/3.),(uu0-tensor_utility::delta()/3.));//uuuu0
	uuuu1=0.5*(diamond(dyad(u,tensor_utility::delta()),u)+diamond(diamond(u,tensor_utility::delta()),u)
	          +dyad(dyad(u,tensor_utility::delta()),u)
	          +dyad(diamond(u,tensor_utility::delta()),u)-4*dyad(uu0,uu0));//uuuu1
	//uuuu2
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				for(int l=0;l<3;l++){
					uuuu2(i,j,k,l)=0.5*(tensor_utility::delta(i,k)*tensor_utility::delta(j,l)
					                   +tensor_utility::delta(j,k)*tensor_utility::delta(i,l)
					                   -tensor_utility::delta(i,j)*tensor_utility::delta(k,l));
				}
	uuuu2+=0.5*(dyad(uu0,tensor_utility::delta())+dyad(tensor_utility::delta(),uu0)
	           -dyad(u,diamond(u,tensor_utility::delta()))
	           -diamond(u,diamond(u,tensor_utility::delta()))
	           -dyad(dyad(u,tensor_utility::delta()),u)
	           -dyad(diamond(u,tensor_utility::delta()),u)+dyad(uu0,uu0));
}
//////////////////////////////
//AxissymmetricElement class//
//////////////////////////////
AxisymmetricResistanceTensor::AxisymmetricResistanceTensor(
	const Vector3d& director,double longAxis)
:base(director),a(longAxis){}
void AxisymmetricResistanceTensor::translate(const Vector3d& r){
	position+=r;
	ResistanceTensor::translate(r);
}
void AxisymmetricResistanceTensor::rotate(const Tensor3x3& Q){
	ResistanceTensor::rotate(Q);
	position=convert(Q,position);
	base.u=convert(Q,base.u);
	base.update();
}
void AxisymmetricResistanceTensor::calcTensors(){
	Tensor3x3 AA=6.*a*PI*(getXA()*base.uu0+getYA()*base.uu1);
	set(AA,
	    Tensor3x3(),
	    8.*a*a*a*PI*(getXC()*base.uu0+getYC()*base.uu1),
	    Tensor3x3x3(),
	    4.*PI*a*a*a*getYH()*base.uuu,
	    20.*PI/3.*a*a*a*(getXK()*base.uuuu0+getYK()*base.uuuu1+getZK()*base.uuuu2)
	);
}
///////////////////////////
//Point class
///////////////////////////
Point::Point(){calcTensors();}
double Point::getXA()const{return 0.0;}
double Point::getYA()const{return 0.0;}
double Point::getXC()const{return 0.0;}
double Point::getYC()const{return 0.0;}
double Point::getYH()const{return 0.0;}
double Point::getXK()const{return 0.0;}
double Point::getYK()const{return 0.0;}
double Point::getZK()const{return 0.0;}
///////////////////////////
//Line class
///////////////////////////
Line::Line(const Vector3d& u,double h)
:AxisymmetricResistanceTensor(u,0),halfLength(h){calcTensors();}
double Line::getXA()const{return 0.0;}
double Line::getYA()const{return 0.0;}
double Line::getXC()const{return 0.0;}
double Line::getYC()const{return 0.0;}
double Line::getYH()const{return 0.0;}
double Line::getXK()const{return 0.0;}
double Line::getYK()const{return 0.0;}
double Line::getZK()const{return 0.0;}
//////////////////////////////
//Sphere class              //
//////////////////////////////
Sphere::Sphere(double radius)
:AxisymmetricResistanceTensor(Vector3d(1,0,0),radius){calcTensors();}
double Sphere::getXA()const{return 1.0;}
double Sphere::getYA()const{return 1.0;}
double Sphere::getXC()const{return 1.0;}
double Sphere::getYC()const{return 1.0;}
double Sphere::getYH()const{return 0.0;}
double Sphere::getXK()const{return 1.0;}
double Sphere::getYK()const{return 1.0;}
double Sphere::getZK()const{return 1.0;}
//////////////////////////////
//Disk class                //
//////////////////////////////
Disk::Disk(const Vector3d& director,double radius)
:AxisymmetricResistanceTensor(director,radius){calcTensors();}
double Disk::getXA()const{return 8./(3.*PI);}
double Disk::getYA()const{return 16./(9.*PI);}
double Disk::getXC()const{return 4./(3.*PI);}
double Disk::getYC()const{return 4./(3.*PI);}
double Disk::getYH()const{return -4./(3.*PI);}
double Disk::getXK()const{return 8./(15.*PI);}
double Disk::getYK()const{return 4./(5.*PI);}
double Disk::getZK()const{return 16./(15.*PI);}
//////////////////////////////
//Needle class              //
//////////////////////////////
Needle::Needle(const Vector3d& director,double a,double _b)
:AxisymmetricResistanceTensor(director,a),b(_b){P=1./log(2.*a*a/b/b);calcTensors();}
double Needle::getXA()const{return 4.*P/(6-3.*P);}
double Needle::getYA()const{return 8.*P/(6.+3.*P);}
double Needle::getXC()const{return 0.0;}
double Needle::getYC()const{return 2.*P/(6.-3.*P);}
double Needle::getYH()const{return 2.*P/(6.-3.*P);}
double Needle::getXK()const{return 4.*P/(30.-45.*P);}
double Needle::getYK()const{return 2.*P/(10.-5.*P);}
double Needle::getZK()const{return 0.0;}
//////////////////////////////
//ProlateSpheroid class     //
//////////////////////////////
ProlateSpheroid::ProlateSpheroid(const Vector3d& director,double a,double _b)
:AxisymmetricResistanceTensor(director,a),b(_b)
{s=sqrt(a*a-b*b)/a;s2=s*s;s3=s2*s;s5=s3*s2;L=log((1.+s)/(1.-s));calcTensors();}
double ProlateSpheroid::getXA()const{
	return 8./3.*s3/(-2.*s+(1.+s2)*L);}
double ProlateSpheroid::getYA()const{
	return 16./3.*s3/(2.*s+(3.*s2-1.)*L);}
double ProlateSpheroid::getXC()const{
	return 4./3.*(s3*(1.-s2))/(2.*s-(1.-s2)*L);}
double ProlateSpheroid::getYC()const{
	return 4./3.*(s3*(2.-s2))/(-2.*s+(1.+s2)*L);}
double ProlateSpheroid::getYH()const{
	return 4./3.*s5/(-2.*s+(1.+s2)*L);}
double ProlateSpheroid::getXK()const{
	return 8./15.*s5/((3.-s2)*L-6.*s);}
double ProlateSpheroid::getYK()const{
	double numerator=s5*(2.*s*(1.-2.*s2)-(1.-s2)*L);
	double denominator=(2.*s*(2.*s2-3.)+3.*(1.-s2)*L)*(-2.*s+(1.+s2)*L);
	return 0.8*numerator/denominator;}
double ProlateSpheroid::getZK()const{
	return 3.2*s5*(1.-s2)/(3.*(1.-s2)*(1.-s2)*L-2.*s*(3.-5.*s2));}
//////////////////////////////
//OblateSpheroid class      //
//////////////////////////////
OblateSpheroid::OblateSpheroid(const Vector3d& director,double a,double _b)
:AxisymmetricResistanceTensor(director,a),b(_b)
{s=sqrt(a*a-b*b)/a;s2=s*s;s3=s*s2;s5=s3*s2;sq=sqrt(1.-s*s);Q=atan(s/sq);calcTensors();}
double OblateSpheroid::getXA()const{
	return 4./3.*s3/((2.*s2-1)*Q+s*sq);}
double OblateSpheroid::getYA()const{
	return 8./3.*s3/((2.*s2+1.)*Q-s*sq);}
double OblateSpheroid::getXC()const{
	return 2./3.*s3/(Q-s*sq);}
double OblateSpheroid::getYC()const{
	return 2./3.*s3*(2.-s2)/(s*sq-(1.-2.*s2)*Q);}
double OblateSpheroid::getYH()const{
	return -2./3.*s5/(s*sq-(1.-2.*s2)*Q);}
double OblateSpheroid::getXK()const{
	return 4./15.*s5/((3.-2.*s2)*Q-3.*s*sq);}
double OblateSpheroid::getYK()const{
	double numerator=s5*(s*(1.+s2)-sq*Q);
	double denominator=(3.*s-s3-3.*sq*Q)*(s*sq-(1.-2*s2)*Q);
	return 0.4*numerator/denominator;}
double OblateSpheroid::getZK()const{
	return 1.6*s5/(3.*Q-(2.*s3+3.*s)*sq);
}
#undef PI

