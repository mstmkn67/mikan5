//AprxVector2dSquareLatticeField2d
//���`��ԕt��Vector2d�^�����i�q2����(xy����)��N���X

#ifndef _APRX_VECTOR2D_SQUARE_LATTICE_FIELD_2D_H_
#define _APRX_VECTOR2D_SQUARE_LATTICE_FIELD_2D_H_

#include "AprxSquareLatticeField2d.h"

class AprxVector2dSquareLatticeField2d:public AprxSquareLatticeField2d<Vector2d>{
public:
	AprxVector2dSquareLatticeField2d();
	AprxVector2dSquareLatticeField2d(
		Array2d<Vector2d>* field,
	  const int* lMin,const Vector2d& sMin,//�n�̍ŏ��̊i�q�ƌn�̍ŏ��l
	  const int* lMax,const Vector2d& sMax);//�n�̍ő�̊i�q�ƌn�̍ő�l
	virtual ~AprxVector2dSquareLatticeField2d();
	virtual double differential(const Vector2d& r,int i,int j);//�ʒur�ɂ�������i������j�����Ŕ���
	virtual double rot(const Vector2d& r);       //�ʒur�ɂ�������`��Ԃ���rot��Ԃ�
	virtual double div(const Vector2d& r);       //�ʒur�ɂ�������`��Ԃ���div��Ԃ�
protected:
private:
	double element(const Vector2d& r,int i)const;//r��i������Ԃ�
};

#endif // _APRX_VECTOR2D_SQUARE_LATTICE_FIELD_2D_H_
