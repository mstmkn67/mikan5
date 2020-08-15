
#ifndef _SUPERPOSITION_MODEL_H_
#define _SUPERPOSITION_MODEL_H_

#include "ParticleModel.h"
#include "../common/resMobTensor/AxisymmetricResistanceTensor.h"

//Superposition analyzer
//Calculating the hydrodynamic tensor using superpositon approximation

class SuperpositionModel:public ParticleModel{
public:
	SuperpositionModel(double viscosity,
	                   const vector<AxisymmetricResistanceTensor*>& r,
										 const string& center,
	                   const string& diagonalizedTensor)
							:ParticleModel(viscosity,center,diagonalizedTensor),ress(r){};
	virtual ~SuperpositionModel(){};
	virtual void setTensor(const vector<AxisymmetricResistanceTensor*>& r){ress=r;};
	virtual void update(Vector3d& translation,Tensor3x3& rotation);
protected:
private:
	vector<AxisymmetricResistanceTensor*> ress;
};

#endif // _SUPERPOSITION_MODEL_H_
