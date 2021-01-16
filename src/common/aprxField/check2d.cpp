#include "AprxSquareLatticeField2d.h"
#include "AprxDoubleSquareLatticeField2d.h"
#include "AprxVector2dSquareLatticeField2d.h"

int main()
{
	Array2d<double> a2(5,5);
	for(int i=0;i<a2.size()[0];i++)
		for(int j=0;j<a2.size()[1];j++)
			a2[i][j]=i;
	Array2d<Vector2d> v2(5,5);
	for(int i=0;i<v2.size()[0];i++)
		for(int j=0;j<v2.size()[1];j++)
			v2[i][j].set(i,1);
	int lmin[2]={1,1},lmax[2]={3,3};
	Vector2d smin(1,1),smax(5,5);
	////////////////////////////
	//AprxSquareLatticeField2d//
	////////////////////////////
	AprxSquareLatticeField2d<double> asf(&a2,lmin,smin,lmax,smax);
	cout << asf(Vector2d(1,3)) << endl;
	//////////////////////////////////
	//AprxDoubleSquareLatticeField2d//
	//////////////////////////////////
	AprxDoubleSquareLatticeField2d adsf(&a2,lmin,smin,lmax,smax);
	cout << adsf(Vector2d(1,3)) << endl;
	cout << adsf.grad(Vector2d(1,3))<< endl;
	////////////////////////////////////
	//AprxVector2dSquareLatticeField2d//
	////////////////////////////////////
	AprxVector2dSquareLatticeField2d avsf(&v2,lmin,smin,lmax,smax);
	cout << avsf(Vector2d(1,3));
	cout << avsf.rot(Vector2d(1,3)) << endl;
	cout << avsf.div(Vector2d(1,3)) << endl;
	return 0;
}
