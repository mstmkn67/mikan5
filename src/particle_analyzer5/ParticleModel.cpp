#include "ParticleModel.h"

void ParticleModel::setCondition(const string& center,const string& d){
	diagonalizedTensor=d;
}

void ParticleModel::update(Vector3d& translate,Tensor3x3& rotate){
	//translate
	if(center=="resistance"){
		translate=res.selectCenterOfResistance();
	}else if(center=="mobility"){
		translate=res.selectCenterOfMobility();
	}
	//rotate
	if(diagonalizedTensor=="A"){
		rotate=res.select_A_diagonalizedSystem();
	}else if(diagonalizedTensor=="B"){
		rotate=res.select_B_diagonalizedSystem();
	}else if(diagonalizedTensor=="C"){
		rotate=res.select_C_diagonalizedSystem();
	}else if(diagonalizedTensor=="a"){
		rotate=res.select_a_diagonalizedSystem();
	}else if(diagonalizedTensor=="b"){
		rotate=res.select_b_diagonalizedSystem();
	}else if(diagonalizedTensor=="c"){
		rotate=res.select_c_diagonalizedSystem();
	}else{
		rotate.set(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0);
	}
	//cout << MobilityTensor(res).get_b() << endl;
}

ResistanceTensor ParticleModel::getResistanceTensor(){
	return res;
}
