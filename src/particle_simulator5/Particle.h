//Particle

#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "../common/EulerAngle.h"
#include "../common/resMobTensor/MobilityTensorOfBrownianSimulation.h"
#include "../common/velocityField/VelocityField.h"
#include "../common/electricMagneticField/ElectricMagneticField.h"

#include <vector>
#include <string>
#include <algorithm>
using namespace std;

class Bead{
public:
	Vector3d o_pos;
	Vector3d position;
	Vector3d velocity;
	Vector3d force;
};

class Particle{
public:
	Particle(int index,
	         MobilityTensorOfBrownianSimulation* mob,
	         const Vector3d& position=Vector3d(0.0,0.0,0.0),
	         const EulerAngle& angle=EulerAngle(0.0,0.0,1.0));
	virtual ~Particle(){};
	virtual void updateEuler();//���Ԕ��W���s��
	virtual void updateRungeKutta1();
	virtual void updateRungeKutta2();
	virtual void updateRungeKutta3();
	virtual void updateRungeKutta4();
	virtual void forceAndTorqueReset();//�͂ƃg���N�̃��Z�b�g
	virtual Tensor3x3 getStress();//���݂̃����o�ŉ��͂��o��

	Vector3d position;       //���S�̈ʒu
	EulerAngle eulerAngle;   //�p�x(�I�C���[�p�A4����)
	Vector3d velocity,angularVelocity;//���q�̑��x�A�p���x
	Vector3d force,torque;//���q�ɓ�����
	Vector3d dRB,dTB;//�u���E���^���ɂ��ψ�
	int index;//index in udf files
	
	//vector<Vector3d>* original_bead_positions;
	vector<Bead*> bead;
	//virtual void beadsInitial();
	virtual void beadsPositionCalc();
	virtual void beadsForceAndTorqueCalc();
	
	static void setDt(double t){dt=t;};
	static void setVelocityField(VelocityField* vf){velocityField=vf;};
	static void setElectricField(ElectricMagneticField* ef){electricField=ef;};
	static void setBrownianFlag(bool f){brownianFlag=f;};
protected:
	virtual void calcDiff(Vector3d& dR,EulerAngle& dq);//updateEuler,updateRungeKutta�Ŏg��
private:
	MobilityTensorOfBrownianSimulation* mob;
	//
	static ElectricMagneticField* electricField;
	//static ElectricMagneticField* magneticField;
	//RungeKutta�p�ϐ�
	Vector3d R0,dR0,dR1,dR2,dR3;
	EulerAngle q0,dq0,dq1,dq2,dq3;
	static VelocityField* velocityField;//���x��
	static double dt;                //���ԍ���
	static bool brownianFlag;//Brown�^�������邩�H
	//���̃N���X�̃R�s�[���R�s�[�R���X�g���N�^�������̂��ߋ֎~
	Particle(const Particle& ob);
	Particle& operator=(const Particle& ob);
};

#endif // _PARTICLE_H_
