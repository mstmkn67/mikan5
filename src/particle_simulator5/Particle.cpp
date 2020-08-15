#include "Particle.h"

VelocityField* Particle::velocityField=0;
double Particle::dt=0.01;
bool Particle::brownianFlag=false;
ElectricMagneticField* Particle::electricField=0;

Particle::Particle(int i,MobilityTensorOfBrownianSimulation* m,const Vector3d& p,const EulerAngle& a)
:index(i),mob(m),position(p),eulerAngle(a){
}

//EulerMethod
void Particle::updateEuler(){
	beadsForceAndTorqueCalc();
	if(brownianFlag)mob->getBrownianDisplacement(dRB,dTB);//�u���E���^���ɂ��ψ�
	Vector3d dR;EulerAngle dq;
	calcDiff(dR,dq);
	position+=dR;//���i�̎��Ԕ��W
	eulerAngle+=dq;//��]�̎��Ԕ��W
	eulerAngle.renormalize();
}

//RungeKuttaMethod
void Particle::updateRungeKutta1(){
	beadsForceAndTorqueCalc();
	if(brownianFlag)mob->getBrownianDisplacement(dRB,dTB);//�u���E���^���ɂ��ψ�
	R0=position;
	q0=eulerAngle;
	calcDiff(dR0,dq0);
	position=R0+dR0;//���i�̎��Ԕ��W
	eulerAngle=q0+dq0;//��]�̎��Ԕ��W
}

void Particle::updateRungeKutta2(){
	beadsForceAndTorqueCalc();
	calcDiff(dR1,dq1);
	position=R0+0.5*dR1;//���i�̎��Ԕ��W
	eulerAngle=q0+0.5*dq1;//��]�̎��Ԕ��W
}

void Particle::updateRungeKutta3(){
	beadsForceAndTorqueCalc();
	calcDiff(dR2,dq2);
	position=R0+0.5*dR2;//���i�̎��Ԕ��W
	eulerAngle=q0+0.5*dq2;//��]�̎��Ԕ��W
}

void Particle::updateRungeKutta4(){
	beadsForceAndTorqueCalc();
	calcDiff(dR3,dq3);
	position=R0+(dR0+2*dR1+2*dR2+dR3)/6.;//���i�̎��Ԕ��W
	eulerAngle=q0+(dq0+2*dq1+2*dq2+dq3)/6.;//��]�̎��Ԕ��W
	eulerAngle.renormalize();
}

void Particle::calcDiff(Vector3d& dR,EulerAngle& dq){
	Tensor3x3 &a=mob->a,&b=mob->b,&c=mob->c;
	Tensor3x3x3 g=mob->g.tilde(),h=mob->h.tilde();
	//�ψʂ̌��ς���
	Tensor3x3 m=eulerAngle.rotateMat();  //��]�s��
	Tensor3x3 mt=eulerAngle.rotateMat_();//��]�s��̓]�u
	Vector3d mforce=m*force;//�͂̍��W�ϊ�
	Vector3d mtorque=m*torque;//�g���N�̍��W�ϊ�
	Tensor3x3 mE=m*velocityField->getRateOfStrainTensor(position)*mt;//���x�ό`�e���\���̍��W�ϊ�
	Vector3d melectric=m*electricField->getField(position);//�d��̍��W�ϊ�
	velocity=velocityField->getVelocity(position)+mt*(a*mforce+b*mtorque+dot2(g,mE));//���i���x�̌v�Z
	dR=velocity*dt+mt*dRB;//���i�ϕ�
	Vector3d angularVelocity_P=m*velocityField->getAngularVelocity(position)
	                           +b.transpose()*mforce+c*mtorque+dot2(h,mE);
	angularVelocity=mt*angularVelocity_P;//��]���x�̌v�Z
	Vector3d dH=angularVelocity_P*dt+dTB;
	dq=eulerAngle.Q_convert(dH);//�I�C���[�p(4����)�̕ϕ�
}

Tensor3x3 Particle::getStress(){///check
	Tensor3x3x3& g=mob->g,&h=mob->h;Tensor3x3x3x3& k=mob->k;
	Vector3d ele=electricField->getField(position);
	//volume*(S+RF+0.5*�ÁET)
	Tensor3x3 sigma;
	sigma=-1.0*(dyad(position,force));//�r���A������
	sigma-=0.5*(Tensor3x3x3()*torque);//�g���N����
	//���͋ɂ���̍�
	Tensor3x3 m=eulerAngle.rotateMat();  //��]�s��
	Tensor3x3 mt=eulerAngle.rotateMat_();//��]�s��̓]�u
	Tensor3x3 E=velocityField->getRateOfStrainTensor(position);//���x�ό`�e���\��
	Tensor3x3 mEmt=m*E*mt;
	Tensor3x3 temp=g*m*force+h*m*torque+dot2(k,mEmt)+mob->get_kBTRh_();
	sigma+=-mt*temp*m;
	return sigma;
}

void Particle::forceAndTorqueReset(){
	force.set(0.0,0.0,0.0);
	torque.set(0.0,0.0,0.0);
	vector<Bead*>::iterator i=bead.begin();
	for(;i!=bead.end();i++){
		(*i)->force.set(0.0,0.0,0.0);
	}
}

//void Particle::beadsInitial(){
//	
//}

void Particle::beadsPositionCalc(){
	vector<Bead*>::iterator i=bead.begin();
	//vector<Vector3d>::iterator a=original_bead_positions->begin();
	for(;i!=bead.end();i++){
		(*i)->position=eulerAngle.rotate_((*i)->o_pos);
		(*i)->velocity=velocity+(angularVelocity^(*i)->position);
		(*i)->position+=position;
	}
}

void Particle::beadsForceAndTorqueCalc(){
	vector<Bead*>::iterator i=bead.begin();
	for(;i!=bead.end();i++){
		force+=(*i)->force;
		torque+=(((*i)->position-position)^(*i)->force);
	}
}
