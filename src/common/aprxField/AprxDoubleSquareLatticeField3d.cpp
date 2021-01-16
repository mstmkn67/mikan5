#include "AprxDoubleSquareLatticeField3d.h"

AprxDoubleSquareLatticeField3d::AprxDoubleSquareLatticeField3d()
	:AprxSquareLatticeField3d<double>(){}

AprxDoubleSquareLatticeField3d::AprxDoubleSquareLatticeField3d(
		Array3d<double>* field,
    const int* lMin,const Vector3d& sMin,
    const int* lMax,const Vector3d& sMax)
	:AprxSquareLatticeField3d<double>(field,lMin,sMin,lMax,sMax){}

AprxDoubleSquareLatticeField3d::~AprxDoubleSquareLatticeField3d(){}

double AprxDoubleSquareLatticeField3d::differential(const Vector3d& r,int i)
{
	static int k[8][3]={
		{ 1, 1, 1},{-1, 1, 1},{-1,-1, 1},{ 1,-1, 1},
		{ 1, 1,-1},{-1, 1,-1},{-1,-1,-1},{ 1,-1,-1}
	};
	double nearField[8];
	Vector3d g;
	getNeighborCoordAndField(r,nearField,g);
	double PhiPrime[8];
	if(i==0){
		for(int l=0;l<8;l++){
			PhiPrime[l]=-0.25*k[l][i]*(1.0-k[l][1]*g.y)*(1.0-k[l][2]*g.z)/dr.x;
		}
	}else if(i==1){
		for(int l=0;l<8;l++){
			PhiPrime[l]=-0.25*k[l][i]*(1.0-k[l][0]*g.x)*(1.0-k[l][2]*g.z)/dr.y;
		}
	}else{
		for(int l=0;l<8;l++){
			PhiPrime[l]=-0.25*k[l][i]*(1.0-k[l][0]*g.x)*(1.0-k[l][1]*g.y)/dr.z;
		}
	}
	double returnValue=0;
	for(int m=0;m<8;m++){
		returnValue+=PhiPrime[m]*nearField[m];
	}
	return returnValue;
}

Vector3d AprxDoubleSquareLatticeField3d::grad(const Vector3d& r)
{
	return Vector3d(differential(r,0),differential(r,1),differential(r,2));
}
