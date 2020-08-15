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
	virtual void updateEuler();//時間発展を行う
	virtual void updateRungeKutta1();
	virtual void updateRungeKutta2();
	virtual void updateRungeKutta3();
	virtual void updateRungeKutta4();
	virtual void forceAndTorqueReset();//力とトルクのリセット
	virtual Tensor3x3 getStress();//現在のメンバで応力を出力

	Vector3d position;       //中心の位置
	EulerAngle eulerAngle;   //角度(オイラー角、4元数)
	Vector3d velocity,angularVelocity;//粒子の速度、角速度
	Vector3d force,torque;//粒子に働く力
	Vector3d dRB,dTB;//ブラウン運動による変位
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
	virtual void calcDiff(Vector3d& dR,EulerAngle& dq);//updateEuler,updateRungeKuttaで使う
private:
	MobilityTensorOfBrownianSimulation* mob;
	//
	static ElectricMagneticField* electricField;
	//static ElectricMagneticField* magneticField;
	//RungeKutta用変数
	Vector3d R0,dR0,dR1,dR2,dR3;
	EulerAngle q0,dq0,dq1,dq2,dq3;
	static VelocityField* velocityField;//速度場
	static double dt;                //時間差分
	static bool brownianFlag;//Brown運動させるか？
	//このクラスのコピーをコピーコンストラクタ未実装のため禁止
	Particle(const Particle& ob);
	Particle& operator=(const Particle& ob);
};

#endif // _PARTICLE_H_
