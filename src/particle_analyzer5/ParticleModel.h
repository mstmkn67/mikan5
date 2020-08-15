//ParticleModel

#ifndef _PARTICLE_MODEL_H_
#define _PARTICLE_MODEL_H_

#include "../common/resMobTensor/ResistanceMobilityTensor.h"
#include <vector>
using namespace std;

class ParticleModel{
public:
	ParticleModel(double eta,const string& center_="null",const string& diagonalizedTensor_="null")
	:viscosity(eta),center(center_),diagonalizedTensor(diagonalizedTensor_){};
	virtual ~ParticleModel(){};
	//diagonalizedTensor --> A,B,C,a,b,c
	virtual void setCondition(const string& center,const string& diagonalizedTensor);
	//output of translate tensor and rotational tensor 
	virtual void update(Vector3d& translate,Tensor3x3& rotate);
	virtual ResistanceTensor getResistanceTensor();
protected:
	ResistanceTensor res;
	string center;
	string diagonalizedTensor;
	double viscosity;
};

#endif // _PARTICLE_MODEL_H_
