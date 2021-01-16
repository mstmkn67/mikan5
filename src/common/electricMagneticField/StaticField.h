//StaticField.h

#ifndef _STATIC_FIELD_H_
#define _STATIC_FIELD_H_

#include "ElectricMagneticField.h"
using namespace std;

class StaticField:public ElectricMagneticField{
public:
	StaticField(const Vector3d& m):field(m){};
	virtual ~StaticField(){};
	virtual void initial(){};
	virtual void update(){};
	virtual Vector3d getField(const Vector3d& position){
		return field;
	}
private:
	Vector3d field;
};

#endif // _STATIC_FIELD_H_
