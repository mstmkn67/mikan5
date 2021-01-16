//ConstantField.h

#ifndef _CONSTANT_FIELD_H_
#define _CONSTANT_FIELD_H_

#include "ElectricMagneticField.h"
using namespace std;

class ConstantField:public ElectricMagneticField{
public:
	ConstantField(const Vector3d& m):field(m){};
	virtual ~ConstantField(){};
	virtual void initial(){};//‰Šú‰»
	virtual void update(){};//Às(ParticleSimulator 1step‚Éˆê‰ñ‚Ìupdate)
	virtual Vector3d getField(const Vector3d& position){
		return field;
	}
private:
	Vector3d field;
};

#endif // _CONSTANT_FIELD_H_
