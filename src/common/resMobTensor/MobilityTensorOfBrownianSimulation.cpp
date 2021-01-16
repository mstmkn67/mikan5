#include "MobilityTensorOfBrownianSimulation.h"

MobilityTensorOfBrownianSimulation::MobilityTensorOfBrownianSimulation(double kBT_,double dt_,const string c)
:kBT(kBT_),dt(dt_),center(c){
	lowerMob=new double[21];
}
MobilityTensorOfBrownianSimulation::MobilityTensorOfBrownianSimulation(
double kBT_,double dt_,const MobilityTensor& m,const string c)
:MobilityTensor(m),kBT(kBT_),dt(dt_),center(c){
	lowerMob=new double[21];
	calcProperties();
}

MobilityTensorOfBrownianSimulation::MobilityTensorOfBrownianSimulation(
double kBT_,double dt_,const ResistanceTensor& r,const string c)
:MobilityTensor(r),kBT(kBT_),dt(dt_),center(c){
	lowerMob=new double[21];
	calcProperties();
}

MobilityTensorOfBrownianSimulation::MobilityTensorOfBrownianSimulation(
	double kBT_,double dt_,
	const Tensor3x3& a,const Tensor3x3& b,const Tensor3x3& c,
	const Tensor3x3x3& g,const Tensor3x3x3& h,const Tensor3x3x3x3& k,const string ce)
:MobilityTensor(a,b,c,g,h,k),kBT(kBT_),dt(dt_),center(ce){
	lowerMob=new double[21];
	calcProperties();
}

MobilityTensorOfBrownianSimulation::~MobilityTensorOfBrownianSimulation(){
	delete[] lowerMob;
}

void MobilityTensorOfBrownianSimulation::set(
	const Tensor3x3& a,const Tensor3x3& b,const Tensor3x3& c,
	const Tensor3x3x3& g,const Tensor3x3x3& h,const Tensor3x3x3x3& k){
	MobilityTensor::set(a,b,c,g,h,k);
	calcProperties();
}

void MobilityTensorOfBrownianSimulation::calcProperties(){
	if(center=="mobility"){
		selectCenterOfMobility();
	}else if(center=="resistance"){
		selectCenterOfResistance();
	}
	Tensor3x3x3 h=get_h();
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			for(int l=0;l<3;l++)
				for(int m=0;m<3;m++)
					//kBTRh_(i,j)+=kBT*(tensor_utility::epsilon(i,m,l)*h(l,m,j)
					//                +tensor_utility::epsilon(j,l,m)*h(m,i,l));
					kBTRh_(i,j)+=kBT*(tensor_utility::epsilon(l,m,i)*h(j,m,l)
					                +tensor_utility::epsilon(j,l,m)*h(m,i,l));
	double t=2*kBT*dt;
	Tensor3x3 aa=t*a,bb=t*b,cc=t*c;
	lowerMob[lap::sindex(0,0)]=aa.x.x;
	lowerMob[lap::sindex(1,0)]=aa.y.x;lowerMob[lap::sindex(1,1)]=aa.y.y;
	lowerMob[lap::sindex(2,0)]=aa.z.x;lowerMob[lap::sindex(2,1)]=aa.z.y;lowerMob[lap::sindex(2,2)]=aa.z.z;
	lowerMob[lap::sindex(3,0)]=bb.x.x;lowerMob[lap::sindex(3,1)]=bb.x.y;lowerMob[lap::sindex(3,2)]=bb.x.z;
	lowerMob[lap::sindex(4,0)]=bb.y.x;lowerMob[lap::sindex(4,1)]=bb.y.y;lowerMob[lap::sindex(4,2)]=bb.y.z;
	lowerMob[lap::sindex(5,0)]=bb.z.x;lowerMob[lap::sindex(5,1)]=bb.z.y;lowerMob[lap::sindex(5,2)]=bb.z.z;
	lowerMob[lap::sindex(3,3)]=cc.x.x;
	lowerMob[lap::sindex(4,3)]=cc.y.x;lowerMob[lap::sindex(4,4)]=cc.y.y;
	lowerMob[lap::sindex(5,3)]=cc.z.x;lowerMob[lap::sindex(5,4)]=cc.z.y;lowerMob[lap::sindex(5,5)]=cc.z.z;
	lap::get_Cholesky_decomposition(6,lowerMob);
}

void MobilityTensorOfBrownianSimulation::getBrownianDisplacement(
Vector3d& dR,Vector3d& dPhi){
	double s[6];
	s[0]=grn();s[1]=grn();s[2]=grn();s[3]=grn();s[4]=grn();s[5]=grn();
	double x[6];
	for(int i=0;i<6;i++){
		x[i]=0.0;
		for(int j=i;j<6;j++){
			x[i]+=lowerMob[lap::sindex(i,j)]*s[j];
		}
	}
	dR.set(x[0],x[1],x[2]);dPhi.set(x[3],x[4],x[5]);
}

Tensor3x3 MobilityTensorOfBrownianSimulation::get_kBTRh_()const{
	return kBTRh_;
}
