#include "AprxVector3dSquareLatticeField3d.h"

AprxVector3dSquareLatticeField3d::AprxVector3dSquareLatticeField3d()
	:AprxSquareLatticeField3d<Vector3d>(){}

AprxVector3dSquareLatticeField3d::AprxVector3dSquareLatticeField3d(
		Array3d<Vector3d>* field,
    const int* lMin,const Vector3d& sMin,
    const int* lMax,const Vector3d& sMax)
	:AprxSquareLatticeField3d<Vector3d>(field,lMin,sMin,lMax,sMax){}

AprxVector3dSquareLatticeField3d::~AprxVector3dSquareLatticeField3d(){}

double AprxVector3dSquareLatticeField3d::differential(const Vector3d& r,int i,int j)
{
	static int k[8][3]={
		{ 1, 1, 1},{-1, 1, 1},{-1,-1, 1},{ 1,-1, 1},
		{ 1, 1,-1},{-1, 1,-1},{-1,-1,-1},{ 1,-1,-1}
	};
	Vector3d nearField[8],g;
	getNeighborCoordAndField(r,nearField,g);
	double PhiPrime[8];
	if(j==0){
		for(int l=0;l<8;l++){
			PhiPrime[l]=-0.25*k[l][j]*(1-k[l][1]*g.y)*(1-k[l][2]*g.z)/dr.x;
		}
	}else if(j==1){
		for(int l=0;l<8;l++){
			PhiPrime[l]=-0.25*k[l][j]*(1-k[l][0]*g.x)*(1-k[l][2]*g.z)/dr.y;
		}
	}else{
		for(int l=0;l<8;l++){
			PhiPrime[l]=-0.25*k[l][j]*(1-k[l][0]*g.x)*(1-k[l][1]*g.y)/dr.z;
		}
	}
	double returnValue=0;
	for(int m=0;m<8;m++){
		returnValue+=PhiPrime[m]*element(nearField[m],i);
	}
	return returnValue;
}

Vector3d AprxVector3dSquareLatticeField3d::rot(const Vector3d& r)
{
	return Vector3d(differential(r,2,1)-differential(r,1,2),
		              differential(r,0,2)-differential(r,2,0),
                  differential(r,1,0)-differential(r,0,1));
}

double AprxVector3dSquareLatticeField3d::div(const Vector3d& r)
{
	return differential(r,0,0)+differential(r,1,1)+differential(r,2,2);
}

double AprxVector3dSquareLatticeField3d::element(const Vector3d& r,int i)const
{
	if(i==0)return r.x;
	if(i==1)return r.y;
	return r.z;
}

