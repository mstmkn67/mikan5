#include "BeadsModel.h"

void BeadsModel::update(Vector3d& translate,Tensor3x3& rotate){
	res.clear();
	int n=pos->size();
	force_tensor=new double[(3*n)*(3*n+1)/2];
	for(int i=0;i<(3*n)*(3*n+1)/2;i++){
		force_tensor[i]=0.0;
	}
	calc_hydro_effect();
	reduce_tensors();
	ParticleModel::update(translate,rotate);
	vector<Vector3d>::iterator i=pos->begin();
	for(;i!=pos->end();i++){
		(*i)=(*i)+translate;
		(*i)=rotate*(*i);
	}

	delete[] force_tensor;
}

void BeadsModel::calc_hydro_effect(){
	cout << "\t\t calculation of hydrodynamic effect...";
	int n=pos->size();
	for(int i=0;i<n;i++){
		for(int j=0;j<=i;j++){
			int index=lap::msindex(i,j,0,0);
			if(i==j){
				force_tensor[index]=mu;
				index=lap::msindex(i,j,1,0);
				force_tensor[index]=0.0;force_tensor[index+1]=mu;
				index=lap::msindex(i,j,2,0);
				force_tensor[index]=0.0;force_tensor[index+1]=0.0;force_tensor[index+2]=mu;
			}else{
				Tensor3x3 h=beads_interaction::getRPYTensor((*radius)[i],(*pos)[i],(*radius)[j],(*pos)[j]);
				//Tensor3x3 h=beads_interaction::getRPYTensor((*pos)[i],(*pos)[j]);
				force_tensor[index]=h.x.x;force_tensor[index+1]=h.x.y;force_tensor[index+2]=h.x.z;
				index=lap::msindex(i,j,1,0);
				force_tensor[index]=h.y.x;force_tensor[index+1]=h.y.y;force_tensor[index+2]=h.y.z;
				index=lap::msindex(i,j,2,0);
				force_tensor[index]=h.z.x;force_tensor[index+1]=h.z.y;force_tensor[index+2]=h.z.z;
			}
		}
	}
	lap::get_inverse_matrix_symmetric(3*n,force_tensor);
	cout << "  done" <<endl;
}

void BeadsModel::reduce_tensors(){
	double PI=3.14159265358979;
	cout << "\t\t reduction of tensors...";
	int n=pos->size();
	Tensor3x3 A,B,C;Tensor3x3x3 G,H;Tensor3x3x3x3 K;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			Tensor3x3 h(force_tensor[lap::msindex(i,j,0,0)],force_tensor[lap::msindex(i,j,0,1)],force_tensor[lap::msindex(i,j,0,2)],
	                force_tensor[lap::msindex(i,j,1,0)],force_tensor[lap::msindex(i,j,1,1)],force_tensor[lap::msindex(i,j,1,2)],
	                force_tensor[lap::msindex(i,j,2,0)],force_tensor[lap::msindex(i,j,2,1)],force_tensor[lap::msindex(i,j,2,2)]);
			Vector3d& ri=(*pos)[i];Vector3d& rj=(*pos)[j];
			A+=h;//A 
			B+=(rj^h);//B
			C+=-ri^h^rj;//C
			G+=0.5*(dyad(ri,h)+diamond(h,ri));//G
			H+=-0.5*(dyad(ri,(h^rj))+diamond((h^rj),ri));//H
			K+=0.25*(diamond(dyad(ri,h),rj)+dyad(dyad(ri,h),rj)+diamond(diamond(ri,h),rj)+diamond(ri,dyad(h,rj)));//K
		}
	}
	res.set(A,B,C,G,H,K);
	cout << "  done" << endl;
}

