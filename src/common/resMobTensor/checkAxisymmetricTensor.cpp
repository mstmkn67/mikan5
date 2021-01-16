#include "AxisymmetricResistanceTensor.h"
#include <vector>
using namespace std;

int main(){
	Vector3d ux(1,0,0);
	vector<AxisymmetricResistanceTensor*> tensor;
	tensor.push_back(new Point);
	tensor.push_back(new Line(ux,1));
	tensor.push_back(new Sphere(1));
	tensor.push_back(new Disk(ux,1));
	tensor.push_back(new Needle(ux,2,1));
	tensor.push_back(new ProlateSpheroid(ux,2,1));
	tensor.push_back(new OblateSpheroid(ux,2,1));
	vector<AxisymmetricResistanceTensor*>::iterator i;
	for(i=tensor.begin();i!=tensor.end();i++)delete *i;
	
	cout << "check of sphere" << endl;
	Sphere* s=new Sphere(1);//put a sphere
	s->translate(Vector3d( 0.0,0.0,1.0));//the sphere moves, but hydrodynamic center is fixed 
	cout << "center of resistance: " << s->getCenterOfResistance() << endl;
	cout << "center of mobility  : " << MobilityTensor(*s).getCenterOfMobility() << endl;
	Tensor3x3 A=s->get_A(),B=s->get_B();
	cout << tensor_utility::calcCenterOfResistance(A,B) << endl;
	MobilityTensor m(*s);
	A=m.get_c();B=m.get_b();
	cout << A << B << endl;
	cout << tensor_utility::calcCenterOfMobility(B,A) << endl;

	cout << endl;
	cout << "check of dumbbell" << endl;
	Sphere* s1=new Sphere(1);//put a sphere whose radius is 1 on origin
	Sphere* s2=new Sphere(2);//put a sphere whose radius is 2 on origin
	s1->translate(Vector3d(0.0,0.0,1.0));//the sphere moves, but hydrodynamic center is fixed 
	s2->translate(Vector3d(0.0,0.0,-1.0));
	ResistanceTensor* dumbbell=new ResistanceTensor(*s1+*s2);//superposition of the spheres
	Vector3d hcr=dumbbell->getCenterOfResistance();
	cout << "original coupling tensor between translation and rotation:" << endl;
	cout << dumbbell->get_B() << endl; 
	cout << "original center of resistance: " << hcr << endl;
	//dumbbell->translate(-hcr);
	dumbbell->selectCenterOfResistance();
	cout << "coupling tensor between translation and rotation after traslation:" << endl;
	cout << dumbbell->get_B() << endl; 
	cout << "center of resistance after traslation: " 
	     << dumbbell->getCenterOfResistance() << endl;
	delete s1,s2,dumbbell;
	
	cout << endl;
	cout << "check of propeller" << endl;
 	Disk* d1=new Disk(Vector3d(1,0,0),1);
	Disk* d2=new Disk(Vector3d(0,1,0),2);
	d1->translate(Vector3d(0,0, 1));
	d2->translate(Vector3d(0,0,-1));
	ResistanceTensor* propeller=new ResistanceTensor(*d1+*d2);
	hcr=propeller->getCenterOfResistance();
	propeller->selectCenterOfResistance();
	cout << "center of resistance: " << propeller->getCenterOfResistance() << endl;
	cout << "center of mobility  : " << MobilityTensor(*propeller).getCenterOfMobility() << endl;
	cout << propeller->select_B_diagonalizedSystem() << endl;
	cout << "tensor B is: " << endl;
	cout << propeller->get_B() << endl;
	cout << endl;
	propeller->selectCenterOfMobility();
	cout << "center of resistance: " << propeller->getCenterOfResistance() << endl;
	cout << "center of mobility  : " << MobilityTensor(*propeller).getCenterOfMobility() << endl;
	cout << propeller->select_b_diagonalizedSystem() << endl;
	cout << "tensor b is: " << endl;
	cout << MobilityTensor(*propeller).get_b() << endl;
	cout << "center of resistance: " << propeller->getCenterOfResistance() << endl;
	cout << "center of mobility  : " << MobilityTensor(*propeller).getCenterOfMobility() << endl;
	cout << endl;
	
	return 0;
}
