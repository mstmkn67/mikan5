#include "AprxSquareLatticeField3d.h"
#include "AprxDoubleSquareLatticeField3d.h"
#include "AprxVector3dSquareLatticeField3d.h"

int main()
{
	Array3d<double> a3(5,5,5);
	for(int i=0;i<a3.size()[0];i++)
		for(int j=0;j<a3.size()[1];j++)
			for(int k=0;k<a3.size()[2];k++)
				a3[i][j][k]=i;
	Array3d<Vector3d> v3(5,5,5);
	for(int i=0;i<v3.size()[0];i++)
		for(int j=0;j<v3.size()[1];j++)
			for(int k=0;k<v3.size()[2];k++)
				v3[i][j][k].set(i,0,0);
	int lmin[3]={1,1,1},lmax[3]={3,3,3};
	Vector3d smin(1,1,1),smax(5,5,5);
	////////////////////////////
	//AprxSquareLatticeField3d//
	////////////////////////////
	AprxSquareLatticeField3d<double> asf(&a3,lmin,smin,lmax,smax);
	cout << asf(Vector3d(1,3,3)) << endl;
	//////////////////////////////////
	//AprxDoubleSquareLatticeField3d//
	//////////////////////////////////
	AprxDoubleSquareLatticeField3d adsf(&a3,lmin,smin,lmax,smax);
	cout << adsf(Vector3d(1,3,3)) << endl;
	cout << adsf.grad(Vector3d(1,3,3)) << endl;
	////////////////////////////////////
	//AprxVector3dSquareLatticeField3d//
	////////////////////////////////////
	AprxVector3dSquareLatticeField3d avsf(&v3,lmin,smin,lmax,smax);
	cout << avsf(Vector3d(1,3,3)) << endl;
	cout << avsf.rot(Vector3d(1,3,3)) << endl;
	cout << avsf.div(Vector3d(1,3,3)) << endl;
	return 0;
}
