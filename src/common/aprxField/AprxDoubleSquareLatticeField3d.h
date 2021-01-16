//AprxDoubleSquareLatticeField3d
//線形補間付き倍精度型正方格子3次元場クラス

#ifndef _APRX_DOUBLE_SQUARE_LATTICE_FIELD_3D_H_
#define _APRX_DOUBLE_SQUARE_LATTICE_FIELD_3D_H_

#include "AprxSquareLatticeField3d.h"

class AprxDoubleSquareLatticeField3d:public AprxSquareLatticeField3d<double>{
public:
	AprxDoubleSquareLatticeField3d();
	AprxDoubleSquareLatticeField3d(
		Array3d<double>* field,
	  const int* lMin,const Vector3d& sMin,//系の最小の格子と系の最小値
	  const int* lMax,const Vector3d& sMax);//系の最大の格子と系の最大値
	virtual ~AprxDoubleSquareLatticeField3d();
	virtual double differential(const Vector3d& r,int i);//位置rにおける場をi(=x,y,z)成分で微分
	virtual Vector3d grad(const Vector3d& r);//位置rにおける場のgradを返す
private:
};

#endif // _APRX_DOUBLE_SQUARE_LATTICE_FIELD_3D_H_
