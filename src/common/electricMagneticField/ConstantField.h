//ConstantField.h

#ifndef _CONSTANT_FIELD_H_
#define _CONSTANT_FIELD_H_

#include "ElectricMagneticField.h"
using namespace std;

class ConstantField:public ElectricMagneticField{
public:
	ConstantField(const Vector3d& m):field(m){};
	virtual ~ConstantField(){};
	virtual void initial(){};//初期化
	virtual void update(){};//実行(ParticleSimulator 1stepに一回のupdate)
	virtual Vector3d getField(const Vector3d& position){
		return field;
	}
private:
	Vector3d field;
};

#endif // _CONSTANT_FIELD_H_
