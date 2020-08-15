//Boundary Element Model

#ifndef _BOUNDARY_ELEMENT_MODEL_H_
#define _BOUNDARY_ELEMENT_MODEL_H_

#include "ParticleModel.h"
#include "../common/solver/LapackFunctions.h"
#include "../common/mesh/Mesh.h"
#include "../common/GreenFunc/StokesGreenFunc.h"
#include <vector>
#include <fstream>
using namespace std;

class BoundaryElementModel:public ParticleModel{
public:
	BoundaryElementModel(double viscosity,Mesh* mesh_,
                       const string& center,const string& diagonalizedTensor);
	virtual ~BoundaryElementModel();
	virtual void update(Vector3d& translate,Tensor3x3& rotate);
protected:
	virtual void calc_hydro_effect();
	virtual void reduce_tensors();
private:
	virtual Tensor3x3 calc_sub_hydro_tensor(const Vector3d& ri,const Vector3d& rj0,const Vector3d& rj1,
                                          const Vector3d& rj2,double s);
	Mesh* mesh;
	double* traction_tensor;//(3*n)*(3*n)
	double* electro_tensor;//n*3*3
};

#endif // _BOUNDARY_ELEMENT_ANALYZER_H_
