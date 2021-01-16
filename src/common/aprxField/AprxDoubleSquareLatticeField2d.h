//AprxDoubleSquareLatticeField2d
//���`��ԕt���{���x�^�����i�q2����(xy����)��N���X

#ifndef _APRX_DOUBLE_SQUARE_LATTICE_FIELD_2D_H_
#define _APRX_DOUBLE_SQUARE_LATTICE_FIELD_2D_H_

#include "AprxSquareLatticeField2d.h"

class AprxDoubleSquareLatticeField2d:public AprxSquareLatticeField2d<double>{
public:
	AprxDoubleSquareLatticeField2d();
	AprxDoubleSquareLatticeField2d(
		Array2d<double>* field,
		const int* lMin,const Vector2d& sMin,//�n�̍ŏ��̊i�q�ƌn�̍ŏ��l
	  const int* lMax,const Vector2d& sMax);//�n�̍ő�̊i�q�ƌn�̍ő�l
	virtual ~AprxDoubleSquareLatticeField2d();
	virtual double differential(const Vector2d& r,int i);//�ʒur�ɂ�������i(=x,y,z)�����Ŕ���
	virtual Vector2d grad(const Vector2d& r);//�ʒur�ɂ�������grad��Ԃ�
private:
};

#endif // _APRX_DOUBLE_SQUARE_LATTICE_FIELD_2D_H_
