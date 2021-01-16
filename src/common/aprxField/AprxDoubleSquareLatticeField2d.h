//AprxDoubleSquareLatticeField2d
//線形補間付き倍精度型正方格子2次元(xy平面)場クラス

#ifndef _APRX_DOUBLE_SQUARE_LATTICE_FIELD_2D_H_
#define _APRX_DOUBLE_SQUARE_LATTICE_FIELD_2D_H_

#include "AprxSquareLatticeField2d.h"

class AprxDoubleSquareLatticeField2d:public AprxSquareLatticeField2d<double>{
public:
	AprxDoubleSquareLatticeField2d();
	AprxDoubleSquareLatticeField2d(
		Array2d<double>* field,
		const int* lMin,const Vector2d& sMin,//系の最小の格子と系の最小値
	  const int* lMax,const Vector2d& sMax);//系の最大の格子と系の最大値
	virtual ~AprxDoubleSquareLatticeField2d();
	virtual double differential(const Vector2d& r,int i);//位置rにおける場をi(=x,y,z)成分で微分
	virtual Vector2d grad(const Vector2d& r);//位置rにおける場のgradを返す
private:
};

#endif // _APRX_DOUBLE_SQUARE_LATTICE_FIELD_2D_H_
