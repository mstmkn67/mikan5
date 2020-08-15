#include "SuperpositionModel.h"

void SuperpositionModel::update(Vector3d& translate,Tensor3x3& rotate){
	res.clear();
	vector<AxisymmetricResistanceTensor*>::iterator i=ress.begin();
	for(;i!=ress.end();i++){res=res+*(*i);}
	res.A*=viscosity;res.B*=viscosity;res.C*=viscosity;
	res.G*=viscosity;res.H*=viscosity;res.K*=viscosity;
	ParticleModel::update(translate,rotate);
	for(i=ress.begin();i!=ress.end();i++){
		(*i)->translate(translate);
		(*i)->rotate(rotate);
	}
}
