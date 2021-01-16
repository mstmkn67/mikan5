#ifndef _MOBILITY_TENSOR_OF_BROWNIAN_SIMULATION_H_
#define _MOBILITY_TENSOR_OF_BROWNIAN_SIMULATION_H_

#include "ResistanceMobilityTensor.h"
//#include "../solver/CholeskyDecomposition.h"
#include "../solver/LapackFunctions.h"
#include "../EulerAngle.h"
#include "../random/random.h"

class MobilityTensorOfBrownianSimulation:public MobilityTensor{
public:
	MobilityTensorOfBrownianSimulation(double kBT=0.0,double dt=0.01,const string center="mobility");
	MobilityTensorOfBrownianSimulation(double kBT,double dt,const MobilityTensor& m,const string center="mobility");
	MobilityTensorOfBrownianSimulation(double kBT,double dt,const ResistanceTensor& r,const string center="mobility");
	MobilityTensorOfBrownianSimulation(double kBT,double dt,
		const Tensor3x3& a,const Tensor3x3& b,const Tensor3x3& c,
		const Tensor3x3x3& g,const Tensor3x3x3& h,const Tensor3x3x3x3& k,const string center="mobility");

	virtual ~MobilityTensorOfBrownianSimulation();
	virtual void set(const Tensor3x3& a,const Tensor3x3& b,const Tensor3x3& c,
	                 const Tensor3x3x3& g,const Tensor3x3x3& h,const Tensor3x3x3x3& k);
	virtual void getBrownianDisplacement(Vector3d& dR,Vector3d& dPhi);
	virtual Tensor3x3 get_kBTRh_()const;
protected:
	virtual void calcProperties();
private:
	double kBT;//Boltzman constant x temperature
	double dt; //time between steps
	Tensor3x3 kBTRh_;//values of calculating stress due to Brownian motion
	GaussianRandomNumber grn;
	double* lowerMob;//Cholesky decomposition of 6x6 matrix composed of a, b and c
	string center;
};

#endif // _MOBILITY_TENSOR_OF_BROWNIAN_SIMULATION_H_
