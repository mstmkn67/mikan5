//AprxSquareLatticeField2d.h
//正方格子で作られている場の線形補間クラス
//2次元

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
	                         const int* lMin,const Vector2d& sMin,//系の最小の格子と系の最小値
	              					 const int* lMax,const Vector2d& sMax);//系の最大の格子と系の最大値)
	virtual ~AprxSquareLatticeField2d();
	
	virtual void setMin(const int* lMin,const Vector2d& sMin);//系の最小の格子と系の最小値
	virtual void setMax(const int* lMax,const Vector2d& sMax);//系の最大の格子と系の最大値
	
	virtual Vector2d getSystemMin()const;     //系の範囲の最小値を返す
	virtual Vector2d getSystemMax()const;     //系の範囲の最大値を返す
	virtual void getSystemLatticeMin(int* s); //系の格子の最小値を返す
	virtual void getSystemLatticeMax(int* s); //系の格子の最大値を返す
	
	virtual void setField(int x,int y,const FieldType& v);//格子点への値の設定
	virtual void setField(const int* l,const FieldType& v);//格子点への値の設定
	
	virtual Vector2d getDr()const;//格子間隔を出力
	virtual void calcDr();//格子間隔を計算
	
	virtual FieldType getField(int x,int y);//格子のindex(x,y)の値を返す
	virtual FieldType getField(const int* l);
	virtual FieldType* operator[](int index);//a[][]の形で書きかえることが出来るようにする
	virtual FieldType operator()(const Vector2d& r);//位置rにおける線形補間した場を返す
	
	virtual Array2d<FieldType>* getField();
	virtual void setField(Array2d<FieldType>* field);
protected:
	Array2d<FieldType>* field;//場
	Vector2d systemSizeMin,systemSizeMax;  //場の範囲の最小値、大きさ
	int systemLatticeMin[2],systemLatticeMax[2];  //場の範囲の格子の最小値、大きさ
	Vector2d dr;//格子の間隔
	//位置rにおける点を囲むうちでもっとも小さい格子点の座標の点を返す
	virtual void getNeighborLatticeIndex(const Vector2d& r,int* s);
	virtual void getNeighborCoordAndField(const Vector2d& r,FieldType* nearField,Vector2d& g);
	virtual Vector2d getCoord(const int* m);//格子点mに相当する座標を返す
private:
	//以下を使えなくする
	AprxSquareLatticeField2d& operator=(const AprxSquareLatticeField2d<FieldType>&);//代入演算子
	AprxSquareLatticeField2d(const AprxSquareLatticeField2d<FieldType>&);           //コピーコンストラクタ
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
	//C++の文法によると、Vector2d a;a=1.0;は、a=Vector2d(1.0,0.0);と解釈される
	//ちなみに、後ろの0.0は、デフォルト値である。
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
