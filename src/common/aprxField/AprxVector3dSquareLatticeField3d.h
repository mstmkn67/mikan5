//AprxVector3dSquareLatticeField3d
//���`��ԕt��Vector3d�^�����i�q3������N���X

#ifndef _APRX_VECTOR3D_SQUARE_LATTICE_FIELD_3D_H_
#define _APRX_VECTOR3D_SQUARE_LATTICE_FIELD_3D_H_

#include "AprxSquareLatticeField3d.h"

class AprxVector3dSquareLatticeField3d:public AprxSquareLatticeField3d<Vector3d>{
public:
	AprxVector3dSquareLatticeField3d();
	AprxVector3dSquareLatticeField3d(
		Array3d<Vector3d>* field,
	  const int* lMin,const Vector3d& sMin,//�n�̍ŏ��̊i�q�ƌn�̍ŏ��l
	  const int* lMax,const Vector3d& sMax);//�n�̍ő�̊i�q�ƌn�̍ő�l
	virtual ~AprxVector3dSquareLatticeField3d();
	virtual double differential(const Vector3d& r,int i,int j);//�ʒur�ɂ�������i������j�����Ŕ���
	virtual Vector3d rot(const Vector3d& r);       //�ʒur�ɂ�������`��Ԃ���rot��Ԃ�
	virtual double div(const Vector3d& r);       //�ʒur�ɂ�������`��Ԃ���div��Ԃ�
protected:
private:
	double element(const Vector3d& r,int i)const;//r��i������Ԃ�
};

#endif // _APRX_VECTOR3D_SQUARE_LATTICE_FIELD_3D_H_
