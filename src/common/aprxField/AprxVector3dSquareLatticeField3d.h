//AprxVector3dSquareLatticeField3d
//線形補間付きVector3d型正方格子3次元場クラス

#ifndef _APRX_VECTOR3D_SQUARE_LATTICE_FIELD_3D_H_
#define _APRX_VECTOR3D_SQUARE_LATTICE_FIELD_3D_H_

#include "AprxSquareLatticeField3d.h"

class AprxVector3dSquareLatticeField3d:public AprxSquareLatticeField3d<Vector3d>{
public:
	AprxVector3dSquareLatticeField3d();
	AprxVector3dSquareLatticeField3d(
		Array3d<Vector3d>* field,
	  const int* lMin,const Vector3d& sMin,//系の最小の格子と系の最小値
	  const int* lMax,const Vector3d& sMax);//系の最大の格子と系の最大値
	virtual ~AprxVector3dSquareLatticeField3d();
	virtual double differential(const Vector3d& r,int i,int j);//位置rにおける場のi成分をj成分で微分
	virtual Vector3d rot(const Vector3d& r);       //位置rにおける線形補間したrotを返す
	virtual double div(const Vector3d& r);       //位置rにおける線形補間したdivを返す
protected:
private:
	double element(const Vector3d& r,int i)const;//rのi成分を返す
};

#endif // _APRX_VECTOR3D_SQUARE_LATTICE_FIELD_3D_H_
