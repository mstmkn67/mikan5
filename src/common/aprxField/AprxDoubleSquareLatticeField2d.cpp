#include "AprxDoubleSquareLatticeField2d.h"

AprxDoubleSquareLatticeField2d::AprxDoubleSquareLatticeField2d()
	:AprxSquareLatticeField2d<double>(){}

AprxDoubleSquareLatticeField2d::AprxDoubleSquareLatticeField2d(
		Array2d<double>* field,
    const int* lMin,const Vector2d& sMin,
    const int* lMax,const Vector2d& sMax)
	:AprxSquareLatticeField2d<double>(field,lMin,sMin,lMax,sMax){}

AprxDoubleSquareLatticeField2d::~AprxDoubleSquareLatticeField2d(){}

double AprxDoubleSquareLatticeField2d::differential(const Vector2d& r,int i)
{
	static int k[4][2]={
		{ 1, 1},{-1, 1},{-1,-1},{ 1,-1}
	};
	double nearField[4];
	Vector2d g;
	getNeighborCoordAndField(r,nearField,g);
	double PhiPrime[4];
	if(i==0){//x‚Å”÷•ª‚Ì‚Æ‚«
		for(int l=0;l<4;l++){
			PhiPrime[l]=-0.5*k[l][i]*(1.0-k[l][1]*g.y)/dr.x;
		}
	}else if(i==1){//y‚Å”÷•ª‚Ì‚Æ‚«
		for(int l=0;l<4;l++){
			PhiPrime[l]=-0.5*k[l][i]*(1.0-k[l][0]*g.x)/dr.y;
		}
	}else{//z‚Å”÷•ª‚Ì‚Æ‚«
		for(int l=0;l<4;l++){
			PhiPrime[l]=0;
		}
	}
	double returnValue=0;
	for(int m=0;m<4;m++){
		returnValue+=PhiPrime[m]*nearField[m];
	}
	return returnValue;
}

Vector2d AprxDoubleSquareLatticeField2d::grad(const Vector2d& r)
{
	return Vector2d(differential(r,0),differential(r,1));
}
