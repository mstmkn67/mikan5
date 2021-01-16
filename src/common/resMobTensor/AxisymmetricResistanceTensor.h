//軸対称の物体の抵抗中心回りの抵抗テンソルの解析解

#ifndef _AXISYMETRIC_RESISTANCE_TENSOR_H_
#define _AXISYMETRIC_RESISTANCE_TENSOR_H_

#include "ResistanceMobilityTensor.h"
#include <string>
using namespace std;

//base tensor of axisymmetric element
struct BaseTensorOfAxisymmetry{
public:
	BaseTensorOfAxisymmetry(const Vector3d& u);//u is unit vector
	virtual ~BaseTensorOfAxisymmetry(){};
	virtual void update();
	Vector3d u;//director
	Tensor3x3 uu0,uu1;
	Tensor3x3x3 uuu;
	Tensor3x3x3x3 uuuu0,uuuu1,uuuu2;
};
class AxisymmetricResistanceTensor:public ResistanceTensor{
public:
	AxisymmetricResistanceTensor(
		const Vector3d& director=Vector3d(1,0,0),double longAxis=1.0);
	virtual ~AxisymmetricResistanceTensor(){};
	BaseTensorOfAxisymmetry base;
	virtual Vector3d getPosition()const{return position;};
	virtual Vector3d getDirector()const{return base.u;};
	virtual void translate(const Vector3d& r);
	virtual void rotate(const Tensor3x3& Q);
	virtual void calcTensors();
protected:
	//HydrodynamicResistanceFunctions
	virtual double getXA()const=0;
	virtual double getYA()const=0;
	virtual double getXC()const=0;
	virtual double getYC()const=0;
	virtual double getYH()const=0;
	virtual double getXK()const=0;
	virtual double getYK()const=0;
	virtual double getZK()const=0;
	double a;//long axis of ellipsoid(radius in cases of sphere and disk)
	Vector3d position;
private:
};
//Point
class Point:public AxisymmetricResistanceTensor{
public:
	Point();
	virtual ~Point(){};
	virtual double getXA()const;
	virtual double getYA()const;
	virtual double getXC()const;
	virtual double getYC()const;
	virtual double getYH()const;
	virtual double getXK()const;
	virtual double getYK()const;
	virtual double getZK()const;
private:
	virtual MobilityTensor getMobilityTensor()const{};//計算できないから使用できないように
};
//Line
class Line:public AxisymmetricResistanceTensor{
public:
	Line(const Vector3d& u,double halfLength);//u is the director of line
	virtual ~Line(){};
	virtual double getXA()const;
	virtual double getYA()const;
	virtual double getXC()const;
	virtual double getYC()const;
	virtual double getYH()const;
	virtual double getXK()const;
	virtual double getYK()const;
	virtual double getZK()const;
private:
	double halfLength;
	virtual MobilityTensor getMobilityTensor()const{};//計算できないから使用できないように
};
//Sphere
class Sphere:public AxisymmetricResistanceTensor{
public:
	Sphere(double radius);
	virtual ~Sphere(){};
protected:
	//HydrodynamicResistanceFunctions
	virtual double getXA()const;
	virtual double getYA()const;
	virtual double getXC()const;
	virtual double getYC()const;
	virtual double getYH()const;
	virtual double getXK()const;
	virtual double getYK()const;
	virtual double getZK()const;
};
//Disk
class Disk:public AxisymmetricResistanceTensor{
public:
	//director is the normal unit vector of disk surface.
	Disk(const Vector3d& director,double radius);
	virtual ~Disk(){};
protected:
	//HydrodynamicResistanceFunctions
	virtual double getXA()const;
	virtual double getYA()const;
	virtual double getXC()const;
	virtual double getYC()const;
	virtual double getYH()const;
	virtual double getXK()const;
	virtual double getYK()const;
	virtual double getZK()const;
};
//Needle
class Needle:public AxisymmetricResistanceTensor{
public:
	//a is the long axis, b is the short axis of very thin prolate spheroid
	Needle(const Vector3d& director,double a,double b);
	virtual ~Needle(){};
protected:
	//HydrodynamicResistanceFunctions
	virtual double getXA()const;
	virtual double getYA()const;
	virtual double getXC()const;
	virtual double getYC()const;
	virtual double getYH()const;
	virtual double getXK()const;
	virtual double getYK()const;
	virtual double getZK()const;
private:
	double b;//short axis
	double P;
};

class ProlateSpheroid:public AxisymmetricResistanceTensor{
public:
	//a is the long axis, b is the short axis
	ProlateSpheroid(const Vector3d& director,double a,double b);
	virtual ~ProlateSpheroid(){};
protected:
	//HydrodynamicResistanceFunctions
	virtual double getXA()const;
	virtual double getYA()const;
	virtual double getXC()const;
	virtual double getYC()const;
	virtual double getYH()const;
	virtual double getXK()const;
	virtual double getYK()const;
	virtual double getZK()const;
private:
	double b;//short axis
	double s,s2,s3,s5;
	double L;
};

class OblateSpheroid:public AxisymmetricResistanceTensor{
public:
	//a is the long axis, b is the short axis
	OblateSpheroid(const Vector3d& director,double a,double b);
	virtual ~OblateSpheroid(){};
protected:
	//HydrodynamicResistanceFunctions
	virtual double getXA()const;
	virtual double getYA()const;
	virtual double getXC()const;
	virtual double getYC()const;
	virtual double getYH()const;
	virtual double getXK()const;
	virtual double getYK()const;
	virtual double getZK()const;
private:
	double b;//short axis
	double s,s2,s3,s5;
	double sq;//sqrt(1-s*s)
	double Q;
};

#endif //_AXISYMETRIC_RESISTANCE_TENSOR_H_
