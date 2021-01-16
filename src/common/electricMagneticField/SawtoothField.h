//SawtoothField.h

#ifndef _SAWTOOTH_FIELD_H_
#define _SAWTOOTH_FIELD_H_

#include "ElectricMagneticField.h"
using namespace std;

class SawtoothField:public ElectricMagneticField{
public:
	SawtoothField(double potential_max,double asymmetric_parameter,
                double length,double switching_time,double dt);
	virtual ~SawtoothField();
	virtual void initial();
	virtual void update();
	virtual Vector3d getField(const Vector3d& position);
private:
	Vector3d field1,field2;
	double aL,L;
	double tau;
	double dt;
};

#endif // _SAWTOOTH_FIELD_H_
