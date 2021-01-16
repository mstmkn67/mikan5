#include "UdfBemFlow.h"

////////////////////////////UdfBemFlow2d//////////////////////////////////
UdfBemFlow2d::UdfBemFlow2d(const string& udfFileName,double l,double t,double e)
:boundary_check(false){
	udf=new UDFManager(udfFileName);
	length_ratio=udf->d("input.unit_parameter.length")/l;
	cout << "\t\tlength ratio: " << length_ratio << endl;
	velocity_ratio=length_ratio*t/udf->d("input.unit_parameter.time");
	cout << "\t\tvelocity ratio: " << velocity_ratio << endl;
	string sp=udf->s("input.point_number_of_Gauss_Legendre_integration");
	int pp=1;
	if(sp=="1")pp=1;else if(sp=="2")pp=2;else if(sp=="3")pp=3;
	else if(sp=="4")pp=4;else if(sp=="5")pp=5;else if(sp=="6")pp=6;
	integral=new ParameterOfLegendreGaussFormula1d(pp);
}
UdfBemFlow2d::~UdfBemFlow2d(){
	delete integral;
	delete mesh;
	delete udf;
}

void UdfBemFlow2d::readData(){
	mesh=new StokesMesh2d;
	//assign vertex
	Location vloc("input.vertex[]");
	for(int i=0;i<udf->size(vloc);i++){
		vloc.next();
		Vector3d r(udf->d(vloc.sub("position.x")),udf->d(vloc.sub("position.y")));
		mesh->vertex.push_back(new Vertex2d(length_ratio*r));
	}
	cout << "number of vertex: " << mesh->vertex.size() << endl;
	//assign edge
	Location floc("input.edge[]");
	Location dloc("output.edge[]");
	for(int j=0;j<udf->size(floc);j++){
		floc.next();dloc.next();
		vector<int>& v=udf->iarray(floc.sub("vertex[]"));
		int id0=udf->getLocation("Vertex",v[0]).getIndex().get()[0];
		int id1=udf->getLocation("Vertex",v[1]).getIndex().get()[0];
		mesh->edge.push_back(new StokesEdge2d(mesh->vertex[id0],mesh->vertex[id1]));
		//read density
		double dx=velocity_ratio/length_ratio*udf->d(dloc.sub("density.x"));
		double dy=velocity_ratio/length_ratio*udf->d(dloc.sub("density.y"));
		((StokesEdge2d*)mesh->edge[j])->density.set(dx,dy);
	}
	cout << "number of edge: " << mesh->edge.size() << endl;
	mesh->calc_position_length_normal();
	if(udf->size("output.additional_information[]")!=0){
		Location loc("output.additional_information[]");
		for(int i=0;i<udf->size("output.additional_information[]");i++){
			loc.next();
			string item=udf->s(loc.sub("item"));
			string value=udf->s(loc.sub("value"));
			if(item=="min_x"){
				cout << item << ": " << value << endl;
				min.x=atof(value.c_str());
				boundary_check==true;
			}else if(item=="min_y"){
				cout << item << ": " << value << endl;
				min.y=atof(value.c_str());
				boundary_check==true;
			}else if(item=="size_x"){
				cout << item << ": " << value << endl;
				size.x=atof(value.c_str());
				boundary_check==true;
			}else if(item=="size_y"){
				cout << item << ": " << value << endl;
				size.y=atof(value.c_str());
				boundary_check==true;
			}
		}
		min*=length_ratio;
		size*=length_ratio;
	}
}

Vector3d UdfBemFlow2d::getVelocity(const Vector3d& coord){
	Vector2d r=reset_coord(coord);
	Vector2d velocity;
	vector<Edge2d*>::iterator i=mesh->edge.begin();
	for(;i!=mesh->edge.end();i++){
		Vector2d& ri0=(*i)->vertex[0]->position;
		Vector2d& ri1=(*i)->vertex[1]->position;
		double length_i=(*i)->length;
		Vector2d& density=((StokesEdge2d*)(*i))->density;
		Tensor3x3 s;
		for(int l=0;l<integral->point;l++){
			Vector2d r=0.5*((1.0-integral->xi[l])*ri0+(integral->xi[l]+1.0)*ri1);
			s+=integral->w[l]*stokes_green_func2d::get_free_G(Vector3d(coord.x,coord.y,0.0),Vector3d(r.x,r.y,0.0));
		}
		Vector3d temp=s*length_i*Vector3d(density.x,density.y,0.0);
		velocity+=Vector2d(temp.x,temp.y);
	}
	return Vector3d(velocity.x,velocity.y,0.0);
}

Vector3d UdfBemFlow2d::getAngularVelocity(const Vector3d& coord){
	Vector2d r=reset_coord(coord);
	vector<Edge2d*>::iterator i=mesh->edge.begin();
	Tensor3x3 temp;
	for(;i!=mesh->edge.end();i++){
		Vector2d& ri0=(*i)->vertex[0]->position;
		Vector2d& ri1=(*i)->vertex[1]->position;
		double length_i=(*i)->length;
		Vector2d& density=((StokesEdge2d*)(*i))->density;
		Tensor3x3x3 s;
		for(int l=0;l<integral->point;l++){
			Vector2d r=0.5*((1.0-integral->xi[l])*ri0+(integral->xi[l]+1.0)*ri1);
			s+=integral->w[l]*stokes_green_func2d::get_grad_free_G(Vector3d(coord.x,coord.y,0.0),Vector3d(r.x,r.y,0.0));
		}
	 temp+=s*length_i*Vector3d(density.x,density.y,0.0);
	}
	return Vector3d(0.0,0.0,0.5*(temp.y.x-temp.x.y));
}

Tensor3x3 UdfBemFlow2d::getRateOfStrainTensor(const Vector3d& coord){
	Vector2d r=reset_coord(coord);
	vector<Edge2d*>::iterator i=mesh->edge.begin();
	Tensor3x3 temp;
	for(;i!=mesh->edge.end();i++){
		Vector2d& ri0=(*i)->vertex[0]->position;
		Vector2d& ri1=(*i)->vertex[1]->position;
		double length_i=(*i)->length;
		Vector2d& density=((StokesEdge2d*)(*i))->density;
		Tensor3x3x3 s;
		for(int l=0;l<integral->point;l++){
			Vector2d r=0.5*((1.0-integral->xi[l])*ri0+(integral->xi[l]+1.0)*ri1);
			s+=integral->w[l]*stokes_green_func2d::get_grad_free_G(Vector3d(coord.x,coord.y,0.0),Vector3d(r.x,r.y,0.0));
		}
		temp+=s*length_i*Vector3d(density.x,density.y,0.0);
	}
	return Tensor3x3(temp.x.x,0.5*(temp.x.y+temp.y.x),0.0,
		               0.5*(temp.y.x+temp.x.y),temp.y.y,0.0,
	                 0.0,0.0,0.0);
}

Vector2d UdfBemFlow2d::reset_coord(const Vector2d& r){
	if(boundary_check==false)return r;
	return Vector2d(r.x-floor((r.x-min.x)/size.x)*size.x,
	                r.y-floor((r.y-min.y)/size.y)*size.y);
}

////////////////////////////UdfBemFlow3d//////////////////////////////////
UdfBemFlow3d::UdfBemFlow3d(const string& udfFileName,double l,double t,double e)
:boundary_check(false){
	udf=new UDFManager(udfFileName);
	length_ratio=udf->d("input.unit_parameter.length")/l;
	cout << "\t\tlength ratio: " << length_ratio << endl;
	velocity_ratio=length_ratio*t/udf->d("input.unit_parameter.time");
	cout << "\t\tvelocity ratio: " << velocity_ratio << endl;
	string sp=udf->s("input.point_number_of_Gauss_Legendre_integration");
	int pp=1;
	if(sp=="1")pp=1;else if(sp=="3")pp=3;else if(sp=="4")pp=4;else if(sp=="7")pp=7;
	integral=new ParameterOfLegendreGaussFormulaTiangle(pp);
}
UdfBemFlow3d::~UdfBemFlow3d(){
	delete integral;
	delete mesh;
	delete udf;
}

void UdfBemFlow3d::readData(){
	mesh=new StokesMesh;
	//assign vertex
	Location vloc("input.vertex[]");
	for(int i=0;i<udf->size(vloc);i++){
		vloc.next();
		Vector3d r(udf->d(vloc.sub("position.x")),udf->d(vloc.sub("position.y")),
		           udf->d(vloc.sub("position.z")));
		mesh->vertex.push_back(new Vertex(length_ratio*r));
	}
	cout << "number of vertex: " << mesh->vertex.size() << endl;
	//assign face
	Location floc("input.face[]");
	Location dloc("output.face[]");
	for(int j=0;j<udf->size(floc);j++){
		floc.next();dloc.next();
		vector<int>& v=udf->iarray(floc.sub("vertex[]"));
		int id0=udf->getLocation("Vertex",v[0]).getIndex().get()[0];
		int id1=udf->getLocation("Vertex",v[1]).getIndex().get()[0];
		int id2=udf->getLocation("Vertex",v[2]).getIndex().get()[0];
		mesh->face.push_back(new StokesFace(mesh->vertex[id0],mesh->vertex[id1],mesh->vertex[id2]));
		//read density
		double dx=velocity_ratio/length_ratio*udf->d(dloc.sub("density.x"));
		double dy=velocity_ratio/length_ratio*udf->d(dloc.sub("density.y"));
		double dz=velocity_ratio/length_ratio*udf->d(dloc.sub("density.z"));
		((StokesFace*)mesh->face[j])->density.set(dx,dy,dz);
	}
	cout << "number of face: " << mesh->face.size() << endl;
	mesh->calc_positions_areas_normals();
	if(udf->size("output.additional_information[]")!=0){
		Location loc("output.additional_information[]");
		for(int i=0;i<udf->size("output.additional_information[]");i++){
			loc.next();
			string item=udf->s(loc.sub("item"));
			string value=udf->s(loc.sub("value"));
			if(item=="min_x"){
				cout << item << ": " << value << endl;
				min.x=atof(value.c_str());
				boundary_check==true;
			}else if(item=="min_y"){
				cout << item << ": " << value << endl;
				min.y=atof(value.c_str());
				boundary_check==true;
			}else if(item=="min_z"){
				cout << item << ": " << value << endl;
				min.z=atof(value.c_str());
				boundary_check==true;
			}else if(item=="size_x"){
				cout << item << ": " << value << endl;
				size.x=atof(value.c_str());
				boundary_check==true;
			}else if(item=="size_y"){
				cout << item << ": " << value << endl;
				size.y=atof(value.c_str());
				boundary_check==true;
			}else if(item=="size_z"){
				cout << item << ": " << value << endl;
				size.z=atof(value.c_str());
				boundary_check==true;
			}
		}
		min*=length_ratio;
		size*=length_ratio;
	}
}

Vector3d UdfBemFlow3d::getVelocity(const Vector3d& coord){
	Vector3d c=reset_coord(coord);
	Vector3d velocity;
	vector<Face*>::iterator i=mesh->face.begin();
	for(;i!=mesh->face.end();i++){
		Vector3d& ri0=(*i)->vertex[0]->position;
		Vector3d& ri1=(*i)->vertex[1]->position;
		Vector3d& ri2=(*i)->vertex[2]->position;
		double area_i=(*i)->area;
		Vector3d& density=((StokesFace*)(*i))->density;
		Tensor3x3 s;
		for(int l=0;l<integral->point;l++){
			Vector3d r=integral->xi[l].x*(ri1-ri0)+integral->xi[l].y*(ri2-ri0)+ri0;
			s+=integral->w[l]*stokes_green_func::get_free_G(c,r);
		}
		velocity+=s*area_i*density;
	}
	return velocity;
}

Vector3d UdfBemFlow3d::getAngularVelocity(const Vector3d& coord){
	Vector3d c=reset_coord(coord);
	vector<Face*>::iterator i=mesh->face.begin();
	Tensor3x3 temp;
	for(;i!=mesh->face.end();i++){
		Vector3d& ri0=(*i)->vertex[0]->position;
		Vector3d& ri1=(*i)->vertex[1]->position;
		Vector3d& ri2=(*i)->vertex[2]->position;
		double area_i=(*i)->area;
		Vector3d& density=((StokesFace*)(*i))->density;
		Tensor3x3x3 s;
		for(int l=0;l<integral->point;l++){
			Vector3d r=integral->xi[l].x*(ri1-ri0)+integral->xi[l].y*(ri2-ri0)+ri0;
			s+=integral->w[l]*stokes_green_func::get_grad_free_G(c,r);
		}
		temp+=s*area_i*density;
	}
	return Vector3d(0.5*(temp.z.y-temp.y.z),0.5*(temp.x.z-temp.z.x),0.5*(temp.y.x-temp.x.y));
}

Tensor3x3 UdfBemFlow3d::getRateOfStrainTensor(const Vector3d& coord){
	Vector3d c=reset_coord(coord);
	vector<Face*>::iterator i=mesh->face.begin();
	Tensor3x3 temp;
	for(;i!=mesh->face.end();i++){
		Vector3d& ri0=(*i)->vertex[0]->position;
		Vector3d& ri1=(*i)->vertex[1]->position;
		Vector3d& ri2=(*i)->vertex[2]->position;
		double area_i=(*i)->area;
		Vector3d& density=((StokesFace*)(*i))->density;
		Tensor3x3x3 s;
		for(int l=0;l<integral->point;l++){
			Vector3d r=integral->xi[l].x*(ri1-ri0)+integral->xi[l].y*(ri2-ri0)+ri0;
			s+=integral->w[l]*stokes_green_func::get_grad_free_G(c,r);
		}
		temp+=s*area_i*density;
	}
	return Tensor3x3(temp.x.x,0.5*(temp.x.y+temp.y.x),0.5*(temp.x.z+temp.z.x),
		               0.5*(temp.y.x+temp.x.y),temp.y.y,0.5*(temp.y.z+temp.z.y),
	                 0.5*(temp.z.x+temp.x.z),0.5*(temp.z.y+temp.y.z),temp.z.z);
}

Vector3d UdfBemFlow3d::reset_coord(const Vector3d& r){
	if(boundary_check==false)return r;
	return Vector3d(r.x-floor((r.x-min.x)/size.x)*size.x,
	                r.y-floor((r.y-min.y)/size.y)*size.y,
	                r.z-floor((r.z-min.z)/size.z)*size.z);
}

