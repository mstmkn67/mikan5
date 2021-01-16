//ElectricMagneticField.h

#ifndef _ELECTRIC_MAGNETIC_FIELD_H_
#define _ELECTRIC_MAGNETIC_FIELD_H_

#include "../Vector3d.h"
using namespace std;

class ElectricMagneticField{
public:
	ElectricMagneticField(){};
	virtual ~ElectricMagneticField(){};
	virtual void initial()=0;
	virtual void update()=0;
	virtual Vector3d getField(const Vector3d& position)=0;
	double* time;
protected:
};

#endif // _ELECTRIC_MAGNETIC_FIELD_H_
