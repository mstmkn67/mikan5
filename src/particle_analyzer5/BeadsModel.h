//BeadsModel

#ifndef _BEADS_MODEL_H_
#define _BEADS_MODEL_H_

#include "ParticleModel.h"
#include "../common/solver/LapackFunctions.h"
#include "../common/GreenFunc/BeadsInteraction.h"
#include <vector>
#include <fstream>
using namespace std;


class BeadsModel:public ParticleModel{
public:
	BeadsModel(double viscosity,
	           vector<double>* radius_, vector<Vector3d>* pos_,
						 const string& center,
	           const string& diagonalizedTensor)
		:ParticleModel(viscosity,center,diagonalizedTensor),radius(radius_),
     pos(pos_),mu(0.053051647697298445256294587790838){};
	virtual ~BeadsModel(){};
	virtual void update(Vector3d& translate,Tensor3x3& rotate);
protected:
	virtual void calc_hydro_effect();
	virtual void reduce_tensors();
	vector<double>* radius;
	vector<Vector3d>* pos;
private:
	double mu;
	double* force_tensor;//(3*n)*(3*n+1)/2
};

#endif // _BEADS_MODEL_H_
