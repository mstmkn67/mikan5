//AprxVector2dSquareLatticeField2d
//線形補間付きVector2d型正方格子2次元(xy平面)場クラス

#ifndef _APRX_VECTOR2D_SQUARE_LATTICE_FIELD_2D_H_
#define _APRX_VECTOR2D_SQUARE_LATTICE_FIELD_2D_H_

#include "AprxSquareLatticeField2d.h"

class AprxVector2dSquareLatticeField2d:public AprxSquareLatticeField2d<Vector2d>{
public:
	AprxVector2dSquareLatticeField2d();
	AprxVector2dSquareLatticeField2d(
		Array2d<Vector2d>* field,
	  const int* lMin,const Vector2d& sMin,//系の最小の格子と系の最小値
	  const int* lMax,const Vector2d& sMax);//系の最大の格子と系の最大値
	virtual ~AprxVector2dSquareLatticeField2d();
	virtual double differential(const Vector2d& r,int i,int j);//位置rにおける場のi成分をj成分で微分
	virtual double rot(const Vector2d& r);       //位置rにおける線形補間したrotを返す
	virtual double div(const Vector2d& r);       //位置rにおける線形補間したdivを返す
protected:
private:
	double element(const Vector2d& r,int i)const;//rのi成分を返す
};

#endif // _APRX_VECTOR2D_SQUARE_LATTICE_FIELD_2D_H_
