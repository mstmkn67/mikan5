//ElectricMagneticTorque

#ifndef _ELECTRIC_MAGNETIC_TORQUE_H_
#define _ELECTRIC_MAGNETIC_TORQUE_H_

#include "Force.h"
#include "../common/Tensor3x3.h"
#include "../common/electricMagneticField/ElectricMagneticField.h"

class DipoleTorque:public Force{
public:
	DipoleTorque(const vector<Particle*>& particles,const Vector3d& dipole,
               ElectricMagneticField* f);
	virtual ~DipoleTorque();
	virtual void update();
	ElectricMagneticField* emf;
private:
	Vector3d dipole;
};

class SusceptibilityTorque:public Force{
public:
	SusceptibilityTorque(const vector<Particle*>& particles,double alpha,const Tensor3x3& chi,
                       ElectricMagneticField* f);
	virtual ~SusceptibilityTorque();
	virtual void update();
private:
	ElectricMagneticField* emf;
	double alpha;
	Tensor3x3 chi;
};

#endif // _ELECTRIC_MAGNETIC_TORQUE_H_
