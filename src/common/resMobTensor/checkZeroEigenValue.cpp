#include "ResistanceMobilityTensor.h"

int main(){
	double eta=6.0*3.1415;
	double eta2=8.0*3.1415;
	ResistanceTensor res;
	res.A.set(eta,0.0,0.0,
	          0.0,eta,0.0,
	          0.0,0.0,eta);
	res.C.set(eta2,0.0,0.0,
	          0.0,eta2,0.0,
	          0.0,0.0,0.0);
	MobilityTensor mob(res);
	cout << mob.a << endl << endl;
	cout << mob.c << endl;
	return 0;
}
