//ParticleSimulator
//

#ifndef _PARTICLE_SIMULATOR_H_
#define _PARTICLE_SIMULATOR_H_

#include <vector>
#include <map>
using namespace std;

#include "Particle.h"
#include "Force.h"
#include "../common/velocityField/VelocityField.h"
#include "../common/electricMagneticField/ElectricMagneticField.h"
#include "../common/Tensor3x3x3x3.h"
#include "../common/Timer.h"
#include "udfmanager.h"

class ParticleSimulator{
public:
	ParticleSimulator(UDFManager* inUdf,UDFManager* outUdf);
	virtual ~ParticleSimulator();
	virtual void run();
protected:
	virtual void run_Euler();//�I�C���[�@
	virtual void run_4thOrderRungeKutta();//4���̃����Q�N�b�^
	virtual void assign_particle();//���q�̍쐬
	virtual void assign_tensor(int id);//�e���\���̊��蓖��
	virtual void assign_force();//�͂̊��蓖�� 
	virtual void assign_electric_torque();//�g���N�����蓖��
	virtual void assign_velocityField();//���x��̊��蓖��
	virtual void assign_electricField();//�d��̊��蓖��
	virtual void report();
private:
	map<int,MobilityTensorOfBrownianSimulation*> tensor;//type_id and tensor
	map<int,vector<Vector3d> > bead_positions;//type_id and bead_positions in particle fixed frame
	map<int,vector<Particle*> > particle;//type_id, and particle
	vector<Force*> force;
	double volume;//volume of system
	
	VelocityField* velocityField;
	ElectricMagneticField* electricField;
	UDFManager* inudf;//input_udf
	UDFManager* outudf;
	Timer timer;
	double time;
	
	//���̃N���X�̃R�s�[�̋֎~������
	ParticleSimulator(const ParticleSimulator& ob);
	ParticleSimulator& operator=(const ParticleSimulator& ob);
};

#endif // _PARTICLE_SIMULATOR_H_
