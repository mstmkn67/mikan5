//AprxSquareLatticeField2d.h
//�����i�q�ō���Ă����̐��`��ԃN���X
//2����

#ifndef _APRX_SQUARE_LATTICE_FIELD_2D_H_
#define _APRX_SQUARE_LATTICE_FIELD_2D_H_

#include "../Array.h"
#include "../Vector2d.h"
#include <iostream>
using namespace std;

template<class FieldType>
class AprxSquareLatticeField2d{
public:
	AprxSquareLatticeField2d();
	AprxSquareLatticeField2d(Array2d<FieldType>* field,
	                         const int* lMin,const Vector2d& sMin,//�n�̍ŏ��̊i�q�ƌn�̍ŏ��l
	              					 const int* lMax,const Vector2d& sMax);//�n�̍ő�̊i�q�ƌn�̍ő�l)
	virtual ~AprxSquareLatticeField2d();
	
	virtual void setMin(const int* lMin,const Vector2d& sMin);//�n�̍ŏ��̊i�q�ƌn�̍ŏ��l
	virtual void setMax(const int* lMax,const Vector2d& sMax);//�n�̍ő�̊i�q�ƌn�̍ő�l
	
	virtual Vector2d getSystemMin()const;     //�n�͈̔͂̍ŏ��l��Ԃ�
	virtual Vector2d getSystemMax()const;     //�n�͈̔͂̍ő�l��Ԃ�
	virtual void getSystemLatticeMin(int* s); //�n�̊i�q�̍ŏ��l��Ԃ�
	virtual void getSystemLatticeMax(int* s); //�n�̊i�q�̍ő�l��Ԃ�
	
	virtual void setField(int x,int y,const FieldType& v);//�i�q�_�ւ̒l�̐ݒ�
	virtual void setField(const int* l,const FieldType& v);//�i�q�_�ւ̒l�̐ݒ�
	
	virtual Vector2d getDr()const;//�i�q�Ԋu���o��
	virtual void calcDr();//�i�q�Ԋu���v�Z
	
	virtual FieldType getField(int x,int y);//�i�q��index(x,y)�̒l��Ԃ�
	virtual FieldType getField(const int* l);
	virtual FieldType* operator[](int index);//a[][]�̌`�ŏ��������邱�Ƃ��o����悤�ɂ���
	virtual FieldType operator()(const Vector2d& r);//�ʒur�ɂ�������`��Ԃ������Ԃ�
	
	virtual Array2d<FieldType>* getField();
	virtual void setField(Array2d<FieldType>* field);
protected:
	Array2d<FieldType>* field;//��
	Vector2d systemSizeMin,systemSizeMax;  //��͈̔͂̍ŏ��l�A�傫��
	int systemLatticeMin[2],systemLatticeMax[2];  //��͈̔͂̊i�q�̍ŏ��l�A�傫��
	Vector2d dr;//�i�q�̊Ԋu
	//�ʒur�ɂ�����_���͂ނ����ł����Ƃ��������i�q�_�̍��W�̓_��Ԃ�
	virtual void getNeighborLatticeIndex(const Vector2d& r,int* s);
	virtual void getNeighborCoordAndField(const Vector2d& r,FieldType* nearField,Vector2d& g);
	virtual Vector2d getCoord(const int* m);//�i�q�_m�ɑ���������W��Ԃ�
private:
	//�ȉ����g���Ȃ�����
	AprxSquareLatticeField2d& operator=(const AprxSquareLatticeField2d<FieldType>&);//������Z�q
	AprxSquareLatticeField2d(const AprxSquareLatticeField2d<FieldType>&);           //�R�s�[�R���X�g���N�^
};

template<class FieldType>
AprxSquareLatticeField2d<FieldType>::AprxSquareLatticeField2d(){}

template<class FieldType>
AprxSquareLatticeField2d<FieldType>::AprxSquareLatticeField2d(
                Array2d<FieldType>* f,
	              const int* lMin,const Vector2d& sMin,
	              const int* lMax,const Vector2d& sMax)
{
	field=f;
	setMin(lMin,sMin);
	setMax(lMax,sMax);
	calcDr();
}

template<class FieldType>
AprxSquareLatticeField2d<FieldType>::~AprxSquareLatticeField2d(){}

template<class FieldType>
void AprxSquareLatticeField2d<FieldType>::setMin(const int* lMin,const Vector2d& sMin)
{
	systemLatticeMin[0]=lMin[0];systemLatticeMin[1]=lMin[1];
	systemSizeMin=sMin;
}

template<class FieldType>
void AprxSquareLatticeField2d<FieldType>::setMax(const int* lMax,const Vector2d& sMax)
{
	systemLatticeMax[0]=lMax[0];systemLatticeMax[1]=lMax[1];
	systemSizeMax=sMax;
}

template <class FieldType>
Vector2d AprxSquareLatticeField2d<FieldType>::getSystemMin()const
{
	return systemSizeMin;
}

template <class FieldType>
Vector2d AprxSquareLatticeField2d<FieldType>::getSystemMax()const
{
	return systemSizeMax;
}

template <class FieldType>
void AprxSquareLatticeField2d<FieldType>::getSystemLatticeMin(int* s)
{
	s[0]=systemLatticeMin[0];
	s[1]=systemLatticeMin[1];
}

template <class FieldType>
void AprxSquareLatticeField2d<FieldType>::getSystemLatticeMax(int* s)
{
	s[0]=systemLatticeMax[0];
	s[1]=systemLatticeMax[1];
}

template <class FieldType>
void AprxSquareLatticeField2d<FieldType>::setField(int x,int y,const FieldType& v)
{
	(*field)[x][y]=v;
}

template <class FieldType>
void AprxSquareLatticeField2d<FieldType>::setField(const int* l,const FieldType& v)
{
	(*field)[l[0]][l[1]]=v;
}

template <class FieldType>
inline Vector2d AprxSquareLatticeField2d<FieldType>::getDr()const
{
	return dr;
}

template <class FieldType>
FieldType AprxSquareLatticeField2d<FieldType>::getField(int x,int y)
{
	return (*field)[x][y];
}

template <class FieldType>
FieldType AprxSquareLatticeField2d<FieldType>::getField(const int* l)
{
	return (*field)[l[0]][l[0]];
}

template <class FieldType>
FieldType* AprxSquareLatticeField2d<FieldType>::operator[](int index)
{
	return (*field)[index];
}

template <class FieldType>
void AprxSquareLatticeField2d<FieldType>::getNeighborLatticeIndex(const Vector2d& r,int* l)
{
	l[0]=int((r.x-systemSizeMin.x)/dr.x+systemLatticeMin[0]);
	l[1]=int((r.y-systemSizeMin.y)/dr.y+systemLatticeMin[1]);
}

template <class FieldType>
Vector2d AprxSquareLatticeField2d<FieldType>::getCoord(const int* i)
{
	return Vector2d(systemSizeMin.x+dr.x*(i[0]-systemLatticeMin[0]),
	                systemSizeMin.y+dr.y*(i[1]-systemLatticeMin[1]));
}

template <class FieldType>
void AprxSquareLatticeField2d<FieldType>
::getNeighborCoordAndField(const Vector2d& r,FieldType* nearField,Vector2d& g)
{
	int nearPoint[2];
	getNeighborLatticeIndex(r,nearPoint);
	Vector2d coord=getCoord(nearPoint);
	//cout << coord << endl;
	nearField[0]=(*field)[nearPoint[0]  ][nearPoint[1]  ];
	//cout << nearPoint[0]+1 << " " << nearPoint[1]+1 << endl;
	nearField[1]=(*field)[nearPoint[0]+1][nearPoint[1]  ];
	nearField[2]=(*field)[nearPoint[0]+1][nearPoint[1]+1];
	nearField[3]=(*field)[nearPoint[0]  ][nearPoint[1]+1];
	g.x=2*(r.x-coord.x)/dr.x-1;
	g.y=2*(r.y-coord.y)/dr.y-1;
}

template <class FieldType>
FieldType AprxSquareLatticeField2d<FieldType>::operator()(const Vector2d& r)
{
	Vector2d g;
	FieldType nearField[4];
	getNeighborCoordAndField(r,nearField,g);
	double Phi[4];
	Phi[0]=(1.0-g.x)*(1.0-g.y)/4.0;Phi[1]=(1.0+g.x)*(1.0-g.y)/4.0;
	Phi[2]=(1.0+g.x)*(1.0+g.y)/4.0;Phi[3]=(1.0-g.x)*(1.0+g.y)/4.0;
	FieldType returnValue;
	returnValue=0;
	//OnePointAdvaice
	//C++�̕��@�ɂ��ƁAVector2d a;a=1.0;�́Aa=Vector2d(1.0,0.0);�Ɖ��߂����
	//���Ȃ݂ɁA����0.0�́A�f�t�H���g�l�ł���B
	for(int i=0;i<4;i++)returnValue+=Phi[i]*nearField[i];
	return returnValue;
}

template <class FieldType>
void AprxSquareLatticeField2d<FieldType>::calcDr()
{
	Vector2d s=systemSizeMax-systemSizeMin;
	int lx=systemLatticeMax[0]-systemLatticeMin[0];
	int ly=systemLatticeMax[1]-systemLatticeMin[1];
	if(lx==0)dr.x=0;
	else dr.x=s.x/lx;
	if(ly==0)dr.y=0;
	else dr.y=s.y/ly;
}

template<class FieldType>
Array2d<FieldType>* AprxSquareLatticeField2d<FieldType>::getField()
{
	return field;
}

template<class FieldType>
void AprxSquareLatticeField2d<FieldType>::setField(Array2d<FieldType>* f)
{
	field=f;
}

#endif // _APRX_SQUARE_LATTICE_FIELD_2D_H_
