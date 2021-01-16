#include "AxisymmetricResistanceTensor.h"
#include "MobilityTensorOfBrownianSimulation.h"
#include <fstream>

int main(){
	Vector3d director(0,0,1);
	double kBT=1.0,dt=0.01;
	int max_ite=1000,report_ite=10;
	ofstream out("trajectory.dat");
	Vector3d::parenOff();
	EulerAngle::parenOff();
	
	//resistance tensor of axisymmetric object
	AxisymmetricResistanceTensor* res_tensor;
	
	//res_tensor=new Sphere(1);//sphere with radius=1
	res_tensor=new ProlateSpheroid(director,2,1);//ellipsoid long axis=2,short one=1
	//mobility tensor of Brownian dynamics
	MobilityTensorOfBrownianSimulation* mob
		=new MobilityTensorOfBrownianSimulation(kBT,dt,*res_tensor);
	
	Vector3d pos;EulerAngle angle;
	Vector3d dR,dPhi;
	for(int i=0;i<max_ite;i++){
		Tensor3x3 mt=angle.rotateMat_();//rotational matrix
		mob->getBrownianDisplacement(dR,dPhi);
		pos+=mt*dR;//translational displacement
		angle+=angle.Q_convert(dPhi);//rotational difference
		angle.renormalize();
		if(i%report_ite==0){
			out << pos << " " << angle.rotate(director) << endl;
		}
	}
	
	out.close();
	delete mob;
	delete res_tensor;
	return 0;
}
