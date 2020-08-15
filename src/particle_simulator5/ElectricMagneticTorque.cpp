#include "ElectricMagneticTorque.h"

DipoleTorque::DipoleTorque(const vector<Particle*>& p,const Vector3d& d,
	ElectricMagneticField* f):Force(p),dipole(d),emf(f){}

DipoleTorque::~DipoleTorque(){}

void DipoleTorque::update(){
	vector<Particle*>::iterator p=particles.begin();
	Vector3d B=emf->getField((*p)->position);
	for(;p!=particles.end();p++){
		EulerAngle& e=(*p)->eulerAngle;
		Tensor3x3 mt=e.rotateMat_();//‰ñ“]s—ñ‚Ì“]’u
		Vector3d d=mt*dipole;
		(*p)->torque+=(d^B);
		//cout << d << " " << B <<  "  " << (d^B) << endl;
	}
}

SusceptibilityTorque::SusceptibilityTorque(const vector<Particle*>& p,double a,const Tensor3x3& c,
	ElectricMagneticField* f):Force(p),alpha(a),chi(c),emf(f){}

SusceptibilityTorque::~SusceptibilityTorque(){}

void SusceptibilityTorque::update(){
	vector<Particle*>::iterator p=particles.begin();
	Vector3d B=emf->getField((*p)->position);
	for(;p!=particles.end();p++){
		EulerAngle& e=(*p)->eulerAngle;
		Tensor3x3 m=e.rotateMat();  //‰ñ“]s—ñ
		Tensor3x3 mt=e.rotateMat_();//‰ñ“]s—ñ‚Ì“]’u
		Tensor3x3 c_chi=mt*chi*m;
		(*p)->torque+=-alpha*(B^(c_chi*B));
	}
}
