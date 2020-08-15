
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
	int* time;//時間
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
	virtual void update();//ずれの計算を行う
	virtual Vector3d difPosition(const Vector3d& r0,const Vector3d& r1);
	virtual void difPositionAndVelocity(const Vector3d& r0,const Vector3d& r1,
	                                    const Vector3d& v0,const Vector3d& v1,
	                                    Vector3d& r01,Vector3d&v01);
	virtual Vector3d getPosition(const Vector3d& r);

	double dx;//シミュレーションセルとイメージセルとのずれ
private:
	double dt;//微小時間
	double vmax;//セルの最大、最小の差におけるずりの速度
	Vector3d min,size;
};

#endif // _BOUNDARY_CONDITION_H_
