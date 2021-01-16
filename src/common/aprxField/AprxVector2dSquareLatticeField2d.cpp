#include "AprxVector2dSquareLatticeField2d.h"

AprxVector2dSquareLatticeField2d::AprxVector2dSquareLatticeField2d()
	:AprxSquareLatticeField2d<Vector2d>()
{
}

AprxVector2dSquareLatticeField2d::AprxVector2dSquareLatticeField2d(
		Array2d<Vector2d>* field,
	  const int* lMin,const Vector2d& sMin,
	  const int* lMax,const Vector2d& sMax)
	:AprxSquareLatticeField2d<Vector2d>(field,lMin,sMin,lMax,sMax){}

AprxVector2dSquareLatticeField2d::~AprxVector2dSquareLatticeField2d(){}

double AprxVector2dSquareLatticeField2d::differential(const Vector2d& r,int i,int j)
{
	static int k[4][2]={
		{ 1, 1},{-1, 1},{-1,-1},{ 1,-1}
	};
	Vector2d nearField[4];
	Vector2d g;
	getNeighborCoordAndField(r,nearField,g);
	double PhiPrime[4];
	if(j==0){//x‚Å”÷•ª‚Ì‚Æ‚«
		for(int l=0;l<4;l++){
			PhiPrime[l]=-0.5*k[l][j]*(1.0-k[l][1]*g.y)/dr.x;
		}
	}else if(j==1){//y‚Å”÷•ª‚Ì‚Æ‚«
		for(int l=0;l<4;l++){
			PhiPrime[l]=-0.5*k[l][j]*(1.0-k[l][0]*g.x)/dr.y;
		}
	}else{//z‚Å”÷•ª‚Ì‚Æ‚«
		for(int l=0;l<4;l++){
			PhiPrime[l]=0;
		}
	}
	double returnValue=0;
	for(int m=0;m<4;m++){
		returnValue+=PhiPrime[m]*element(nearField[m],i);
	}
	return returnValue;
}

double AprxVector2dSquareLatticeField2d::rot(const Vector2d& r)
{
	return differential(r,1,0)-differential(r,0,1);
}

double AprxVector2dSquareLatticeField2d::div(const Vector2d& r)
{
	return differential(r,0,0)+differential(r,1,1);
}

double AprxVector2dSquareLatticeField2d::element(const Vector2d& r,int i)const
{
	if(i==0)return r.x;
	return r.y;
}

