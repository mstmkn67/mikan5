//AprxDoubleSquareLatticeField3d
//���`��ԕt���{���x�^�����i�q3������N���X

#ifndef _APRX_DOUBLE_SQUARE_LATTICE_FIELD_3D_H_
#define _APRX_DOUBLE_SQUARE_LATTICE_FIELD_3D_H_

#include "AprxSquareLatticeField3d.h"

class AprxDoubleSquareLatticeField3d:public AprxSquareLatticeField3d<double>{
public:
	AprxDoubleSquareLatticeField3d();
	AprxDoubleSquareLatticeField3d(
		Array3d<double>* field,
	  const int* lMin,const Vector3d& sMin,//�n�̍ŏ��̊i�q�ƌn�̍ŏ��l
	  const int* lMax,const Vector3d& sMax);//�n�̍ő�̊i�q�ƌn�̍ő�l
	virtual ~AprxDoubleSquareLatticeField3d();
	virtual double differential(const Vector3d& r,int i);//�ʒur�ɂ�������i(=x,y,z)�����Ŕ���
	virtual Vector3d grad(const Vector3d& r);//�ʒur�ɂ�������grad��Ԃ�
private:
};

#endif // _APRX_DOUBLE_SQUARE_LATTICE_FIELD_3D_H_
