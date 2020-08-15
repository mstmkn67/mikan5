
#ifndef _BOUNDARY_CONDITION_H_
#define _BOUNDARY_CONDITION_H_

#include "../common/Vector3d.h"

class BoundaryCondition{
public:
	BoundaryCondition(){};
	virtual ~BoundaryCondition(){};
	virtual void update(){};
	virtual Vector3d difPosition(const Vector3d& r0,const Vector3d& r1);
	virtual void difPositionAndVelocity(const Vector3d& r0,const Vector3d& r1,
	                                    const Vector3d& v0,const Vector3d& v1,
	                                    Vector3d& r01,Vector3d&v01);
	virtual Vector3d getPosition(const Vector3d& r);
	int* time;//����
private:
};

class SimplePeriodic:public BoundaryCondition{
public:
	SimplePeriodic(const Vector3d& minSize,const Vector3d& size);
	virtual ~SimplePeriodic(){};
	virtual void update(){};
	virtual Vector3d difPosition(const Vector3d& r0,const Vector3d& r1);
	virtual void difPositionAndVelocity(const Vector3d& r0,const Vector3d& r1,
	                                    const Vector3d& v0,const Vector3d& v1,
	                                    Vector3d& r01,Vector3d&v01);
	virtual Vector3d getPosition(const Vector3d& r);
private:
	Vector3d min,size;
};

class LeesEdwards:public BoundaryCondition{
public:
	LeesEdwards(const Vector3d& minSize,const Vector3d& size,double dt);
	virtual ~LeesEdwards();
	virtual void update();//����̌v�Z���s��
	virtual Vector3d difPosition(const Vector3d& r0,const Vector3d& r1);
	virtual void difPositionAndVelocity(const Vector3d& r0,const Vector3d& r1,
	                                    const Vector3d& v0,const Vector3d& v1,
	                                    Vector3d& r01,Vector3d&v01);
	virtual Vector3d getPosition(const Vector3d& r);

	double dx;//�V�~�����[�V�����Z���ƃC���[�W�Z���Ƃ̂���
private:
	double dt;//��������
	double vmax;//�Z���̍ő�A�ŏ��̍��ɂ����邸��̑��x
	Vector3d min,size;
};

#endif // _BOUNDARY_CONDITION_H_
