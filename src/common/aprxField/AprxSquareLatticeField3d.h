//AprxSquareLatticeField3d.h
//正方格子で作られている場の線形補間クラス
//3次元

#ifndef _APRX_SQUARE_LATTICE_FIELD_3D_H_
#define _APRX_SQUARE_LATTICE_FIELD_3D_H_

#include "../Array.h"
#include "../Vector3d.h"
#include <iostream>
using namespace std;

template<class FieldType>
class AprxSquareLatticeField3d{
public:
	AprxSquareLatticeField3d();
	AprxSquareLatticeField3d(Array3d<FieldType>* field,
	              const int* lMin,const Vector3d& sMin,//系の最小の格子と系の最小値
	              const int* lMax,const Vector3d& sMax);//系の最大の格子と系の最大値)
	virtual ~AprxSquareLatticeField3d();

	virtual void setMin(const int* lMin,const Vector3d& sMin);//系の最小の格子と系の最小値
	virtual void setMax(const int* lMax,const Vector3d& sMax);//系の最大の格子と系の最大値
	
	virtual Vector3d getSystemMin()const;     //系の範囲の最小値を返す
	virtual Vector3d getSystemMax()const;     //系の範囲の最大値を返す
	virtual void getSystemLatticeMin(int* s); //系の格子の最小値を返す
	virtual void getSystemLatticeMax(int* s); //系の格子の最大値を返す
	
	virtual void setField(int x,int y,int z,const FieldType& v);//場の設定
	virtual void setField(const int* l,const FieldType& v);   //場の設定
	
	virtual Vector3d getDr()const;//格子間隔を出力
	virtual void calcDr();//格子間隔を計算
	
	virtual FieldType getField(int x,int y,int z);//格子のindex(x,y,z)の値を返す
	virtual FieldType getField(const int* l);
	virtual FieldType** operator[](int index);//a[][][]の形で書きかえることが出来るようにする
	virtual FieldType operator()(const Vector3d& r);//位置rにおける線形補間した場を返す
	
	virtual Array3d<FieldType>* getField();
	virtual void setField(Array3d<FieldType>* field);
protected:
	Array3d<FieldType>* field;//場
	Vector3d systemSizeMin,systemSizeMax;  //場の範囲の最小値、大きさ
	int systemLatticeMin[3],systemLatticeMax[3];  //場の範囲の格子の最小値、大きさ
	Vector3d dr;//格子の間隔
	//位置rにおける点を囲むうちでもっとも小さい格子点の座標の点を返す
	virtual void getNeighborLatticeIndex(const Vector3d& r,int* s);
	virtual void getNeighborCoordAndField(const Vector3d& r,FieldType* nearField,Vector3d& g);
	virtual Vector3d getCoord(const int* m)const;//格子点mに相当する座標を返す
private:
	//以下を使えなくする
	AprxSquareLatticeField3d& operator=(const AprxSquareLatticeField3d<FieldType>&);//代入演算子
	AprxSquareLatticeField3d(const AprxSquareLatticeField3d<FieldType>&);           //コピーコンストラクタ
};

template<class FieldType>
AprxSquareLatticeField3d<FieldType>::AprxSquareLatticeField3d(){}

template<class FieldType>
AprxSquareLatticeField3d<FieldType>::AprxSquareLatticeField3d(
                Array3d<FieldType>* f,
	              const int* lMin,const Vector3d& sMin,
	              const int* lMax,const Vector3d& sMax)
{
	field=f;
	setMin(lMin,sMin);
	setMax(lMax,sMax);
	calcDr();
}

template<class FieldType>
AprxSquareLatticeField3d<FieldType>::~AprxSquareLatticeField3d(){}

template<class FieldType>
void AprxSquareLatticeField3d<FieldType>::setMin(const int* lMin,const Vector3d& sMin)
{
	systemLatticeMin[0]=lMin[0];systemLatticeMin[1]=lMin[1];systemLatticeMin[2]=lMin[2];
	systemSizeMin=sMin;
}

template<class FieldType>
void AprxSquareLatticeField3d<FieldType>::setMax(const int* lMax,const Vector3d& sMax)
{
	systemLatticeMax[0]=lMax[0];systemLatticeMax[1]=lMax[1];systemLatticeMax[2]=lMax[2];
	systemSizeMax=sMax;
}

template <class FieldType>
Vector3d AprxSquareLatticeField3d<FieldType>::getSystemMin()const
{
	return systemSizeMin;
}

template <class FieldType>
Vector3d AprxSquareLatticeField3d<FieldType>::getSystemMax()const
{
	return systemSizeMax;
}

template <class FieldType>
void AprxSquareLatticeField3d<FieldType>::getSystemLatticeMin(int* l)
{
	l[0]=systemLatticeMin[0];
	l[1]=systemLatticeMin[1];
	l[2]=systemLatticeMin[2];
}

template <class FieldType>
void AprxSquareLatticeField3d<FieldType>::getSystemLatticeMax(int* l)
{
	l[0]=systemLatticeMax[0];
	l[1]=systemLatticeMax[1];
	l[2]=systemLatticeMax[2];
}

template <class FieldType>
void AprxSquareLatticeField3d<FieldType>::setField(int x,int y,int z,const FieldType& v)
{
	(*field)[x][y][z]=v;
}

template <class FieldType>
void AprxSquareLatticeField3d<FieldType>::setField(const int* l,const FieldType& v)
{
	(*field)[l[0]][l[1]][l[2]]=v;
}

template <class FieldType>
inline Vector3d AprxSquareLatticeField3d<FieldType>::getDr()const
{
	return dr;
}

template <class FieldType>
FieldType AprxSquareLatticeField3d<FieldType>::getField(int x,int y,int z)
{
	return (*field)[x][y][z];
}

template <class FieldType>
FieldType AprxSquareLatticeField3d<FieldType>::getField(const int* l)
{
	return (*field)[l[0]][l[1]][l[2]];
}

template <class FieldType>
FieldType** AprxSquareLatticeField3d<FieldType>::operator[](int index)
{
	return (*field)[index];
}

template <class FieldType>
void AprxSquareLatticeField3d<FieldType>::getNeighborLatticeIndex(const Vector3d& r,int* l)
{
	l[0]=int((r.x-systemSizeMin.x)/dr.x+systemLatticeMin[0]);
	l[1]=int((r.y-systemSizeMin.y)/dr.y+systemLatticeMin[1]);
	l[2]=int((r.z-systemSizeMin.z)/dr.z+systemLatticeMin[2]);
}

template <class FieldType>
Vector3d AprxSquareLatticeField3d<FieldType>::getCoord(const int* i)const
{
	return Vector3d(systemSizeMin.x+dr.x*(i[0]-systemLatticeMin[0]),
	                systemSizeMin.y+dr.y*(i[1]-systemLatticeMin[1]),
	                systemSizeMin.z+dr.z*(i[2]-systemLatticeMin[2]));
}

template <class FieldType>
void AprxSquareLatticeField3d<FieldType>
::getNeighborCoordAndField(const Vector3d& r,FieldType* nearField,Vector3d& g)
{
	int nearPoint[3];
	getNeighborLatticeIndex(r,nearPoint);
	Vector3d coord=getCoord(nearPoint);
	nearField[0]=(*field)[nearPoint[0]  ][nearPoint[1]  ][nearPoint[2]  ];
	nearField[1]=(*field)[nearPoint[0]+1][nearPoint[1]  ][nearPoint[2]  ];
	nearField[2]=(*field)[nearPoint[0]+1][nearPoint[1]+1][nearPoint[2]  ];
	nearField[3]=(*field)[nearPoint[0]  ][nearPoint[1]+1][nearPoint[2]  ];
	nearField[4]=(*field)[nearPoint[0]  ][nearPoint[1]  ][nearPoint[2]+1];
	nearField[5]=(*field)[nearPoint[0]+1][nearPoint[1]  ][nearPoint[2]+1];
	nearField[6]=(*field)[nearPoint[0]+1][nearPoint[1]+1][nearPoint[2]+1];
	nearField[7]=(*field)[nearPoint[0]  ][nearPoint[1]+1][nearPoint[2]+1];
	g.x=2*(r.x-coord.x)/dr.x-1;
	g.y=2*(r.y-coord.y)/dr.y-1;
	g.z=2*(r.z-coord.z)/dr.z-1;
}

template <class FieldType>
FieldType AprxSquareLatticeField3d<FieldType>::operator()(const Vector3d& r)
{
	Vector3d g;
	FieldType nearField[8];
	getNeighborCoordAndField(r,nearField,g);
	double Phi[8];
	Phi[0]=(1.0-g.x  )*(1.0-g.y )*(1.0-g.z)/8.0;Phi[1]=(1.0+g.x  )*(1.0-g.y )*(1.0-g.z)/8.0;
	Phi[2]=(1.0+g.x  )*(1.0+g.y )*(1.0-g.z)/8.0;Phi[3]=(1.0-g.x  )*(1.0+g.y )*(1.0-g.z)/8.0;
	Phi[4]=(1.0-g.x  )*(1.0-g.y )*(1.0+g.z)/8.0;Phi[5]=(1.0+g.x  )*(1.0-g.y )*(1.0+g.z)/8.0;
	Phi[6]=(1.0+g.x  )*(1.0+g.y )*(1.0+g.z)/8.0;Phi[7]=(1.0-g.x  )*(1.0+g.y )*(1.0+g.z)/8.0;
	FieldType returnValue;
	returnValue=0;
	//OnePointAdvaice
	//C++の文法によると、Vector3d a;a=1.0;は、a=Vector3d(1.0,0.0,0.0);と解釈される
	//ちなみに、後ろ二つの0.0は、デフォルト値である。
	for(int i=0;i<8;i++)returnValue+=Phi[i]*nearField[i];
	return returnValue;
}

template <class FieldType>
void AprxSquareLatticeField3d<FieldType>::calcDr()
{
	Vector3d s=systemSizeMax-systemSizeMin;
	int lx=systemLatticeMax[0]-systemLatticeMin[0];
	int ly=systemLatticeMax[1]-systemLatticeMin[1];
	int lz=systemLatticeMax[2]-systemLatticeMin[2];
	if(lx==0)dr.x=0;
	else dr.x=s.x/lx;
	if(ly==0)dr.y=0;
	else dr.y=s.y/ly;
	if(lz==0)dr.z=0;
	else dr.z=s.z/lz;
}

template<class FieldType>
Array3d<FieldType>* AprxSquareLatticeField3d<FieldType>::getField()
{
	return field;
}

template<class FieldType>
void AprxSquareLatticeField3d<FieldType>::setField(Array3d<FieldType>* f)
{
	field=f;
}


#endif // _APRX_SQUARE_LATTICE_FIELD_3D_H_
