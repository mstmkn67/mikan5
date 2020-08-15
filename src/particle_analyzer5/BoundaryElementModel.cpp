#include "BoundaryElementModel.h"

BoundaryElementModel::BoundaryElementModel(double viscosity,Mesh* mesh_,
	const string& center,const string& diagonalizedTensor)
:ParticleModel(viscosity,center,diagonalizedTensor),mesh(mesh_){}

BoundaryElementModel::~BoundaryElementModel(){}

void BoundaryElementModel::update(Vector3d& translate,Tensor3x3& rotate){
	res.clear();
	int n=mesh->face.size();
	traction_tensor=new double[3*n*3*n];
	mesh->calc_positions_areas_normals();
	calc_hydro_effect();
	reduce_tensors();
	ParticleModel::update(translate,rotate);
	vector<Vertex*>::iterator i=mesh->vertex.begin();
	for(;i!=mesh->vertex.end();i++){
		(*i)->position+=translate;
		(*i)->position=rotate*((*i)->position);
	}
	delete[] traction_tensor;
}

void BoundaryElementModel::calc_hydro_effect(){
	cout << "\t\t calculation of hydrodynamic effect...";
	int n=mesh->face.size();
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			Tensor3x3 H;
			Vector3d& ri=mesh->face[i]->position;
			Vector3d& rj0=mesh->face[j]->vertex[0]->position;
			Vector3d& rj1=mesh->face[j]->vertex[1]->position;
			Vector3d& rj2=mesh->face[j]->vertex[2]->position;
			if(i!=j){
				double sj=mesh->face[j]->area;
				H=calc_sub_hydro_tensor(ri,rj0,rj1,rj2,sj);
			}else{
				Vector3d u0=rj1-rj0,v0=ri-rj0;double s0=0.5*sqrt((u0^v0)*(u0^v0));
				Vector3d u1=rj2-rj1,v1=ri-rj1;double s1=0.5*sqrt((u1^v1)*(u1^v1));
				Vector3d u2=rj0-rj2,v2=ri-rj2;double s2=0.5*sqrt((u2^v2)*(u2^v2));
				H=(calc_sub_hydro_tensor(ri,ri,rj1,rj0,s0)+calc_sub_hydro_tensor(ri,ri,rj1,rj2,s1)
				  +calc_sub_hydro_tensor(ri,ri,rj2,rj0,s2));
			}
			for(int a=0;a<3;a++){
				for(int b=0;b<3;b++){
					traction_tensor[lap::mindex_(i,j,a,b,n)]=H.get(a,b);
				}
			}
		}
	}
	lap::get_inverse_matrix_general(3*n,traction_tensor);
	cout << "  done" <<endl;
}

void BoundaryElementModel::reduce_tensors(){
	cout << "\t\t reduction of tensors...";
	int n=mesh->face.size();
	for(int i=0;i<n;i++){
		double s=mesh->face[i]->area;
		Vector3d ri=mesh->face[i]->position;
		for(int j=0;j<n;j++){
			Vector3d rj=mesh->face[j]->position;
			//Tensor3x3 f(traction_tensor[(3*i+0)*3*n+3*j+0],traction_tensor[(3*i+0)*3*n+3*j+1],traction_tensor[(3*i+0)*3*n+3*j+2],
			//            traction_tensor[(3*i+1)*3*n+3*j+0],traction_tensor[(3*i+1)*3*n+3*j+1],traction_tensor[(3*i+1)*3*n+3*j+2],
			//            traction_tensor[(3*i+2)*3*n+3*j+0],traction_tensor[(3*i+2)*3*n+3*j+1],traction_tensor[(3*i+2)*3*n+3*j+2]);
			Tensor3x3 f(traction_tensor[lap::mindex_(i,j,0,0,n)],traction_tensor[lap::mindex_(i,j,0,1,n)],traction_tensor[lap::mindex_(i,j,0,2,n)],
			            traction_tensor[lap::mindex_(i,j,1,0,n)],traction_tensor[lap::mindex_(i,j,1,1,n)],traction_tensor[lap::mindex_(i,j,1,2,n)],
			            traction_tensor[lap::mindex_(i,j,2,0,n)],traction_tensor[lap::mindex_(i,j,2,1,n)],traction_tensor[lap::mindex_(i,j,2,2,n)]);
			res.A+=f*s;
			res.B+=(ri^f)*s;
			res.C+=-(ri^f^rj)*s;
			res.G+=s*0.5*(dyad(ri,f)+diamond(f,ri));
			res.H+=-s*0.5*(dyad(ri,(f^rj))+diamond((f^rj),ri));
			res.K+=s*0.25*(diamond(dyad(ri,f),rj)+dyad(dyad(ri,f),rj)
			          +diamond(diamond(ri,f),rj)+diamond(ri,dyad(f,rj)));
		}
	}
	cout << "  done" << endl;
}


Tensor3x3 BoundaryElementModel::calc_sub_hydro_tensor(
	const Vector3d& ri,const Vector3d& rj0,const Vector3d& rj1,const Vector3d& rj2,double s){
	return s/viscosity*stokes_green_func::get_free_G(ri,(rj0+rj1+rj2)/3.0);
}
