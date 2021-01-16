//VelocityField.h
// �����Ɋւ�����N���X

#ifndef _VELOCITY_FIELD_H_
#define _VELOCITY_FIELD_H_

#include "../Vector3d.h"
#include "../Tensor3x3.h"
using namespace std;

class VelocityField{
public:
	VelocityField(){};
	virtual ~VelocityField(){};
	virtual void initial()=0;//������
	virtual void update()=0;//���s(ParticleSimulator 1step�Ɉ���update)
	virtual Vector3d getVelocity(const Vector3d& coord)=0;//��̑��x
	virtual Vector3d getAngularVelocity(const Vector3d& coord)=0;//��̊p���x
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord)=0;//��̕ό`���x�e���\��
	//virtual double getShearStrain(){return 0.0;}//�Ђ��݂�Ԃ�(LeesEdwards)
	double* time;
protected:
};

#endif // _VELOCITY_FIELD_H_
