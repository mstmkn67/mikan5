#include "ParticleAnalyzer.h"

#include "../common/Timer.h"

ParticleAnalyzer::ParticleAnalyzer(UDFManager* in,UDFManager* out)
:in_udf(in),out_udf(out)
{}

ParticleAnalyzer::~ParticleAnalyzer(){}

void ParticleAnalyzer::update(){
	loc_type.seek("particle_type[]");
	for(int l=0;l<in_udf->size(loc_type);l++){
		loc_type.next();
		cout << loc_type.str() << endl;
		string type=in_udf->s(loc_type.sub("model.model_type"));
		cout << "model_type: " << type << endl;
		if(type=="beads"){
			analyze_beads();
		}else if(type=="BEM"){
			analyze_BEM();
		}else if(type=="superposition"){
			analyze_superposition();
		}else{
			cout << "We don't know " << type << "!!!" << endl;
			cout << "We ignore ParticleType[" << l << "]." << endl;
		}
	}
}


void ParticleAnalyzer::analyze_beads(){
	vector<double> radii;
	vector<Vector3d> positions;
	Location bd(loc_type.sub("model.beads.bead[]"));
	for(int b=0;b<in_udf->size(bd);b++){
		bd.next();
		double a=in_udf->d(bd.sub("radius"));
		double x=in_udf->d(bd.sub("position.x"));
		double y=in_udf->d(bd.sub("position.y"));
		double z=in_udf->d(bd.sub("position.z"));
		radii.push_back(a);
		positions.push_back(Vector3d(x,y,z));//Set bead position
	}
	cout << "\tnumber of beads: " << positions.size() << endl;
	//string c=in_udf->s(loc_type.sub("center"));
	string c("resistance");
	//cout << "\tcenter of tensor: " << c << endl;
	string d=in_udf->s(loc_type.sub("tensor_diagonalization"));
	cout << "\ttensor_diagonalization: " << d << endl;
	double viscosity=in_udf->d("fluid.viscosity");
	BeadsModel* hba=new BeadsModel(viscosity,&radii,&positions,c,d);
	Vector3d translate;Tensor3x3 rotate;
	Timer timer;
	hba->update(translate,rotate);
	double calcTime=timer.get();
	//Bead
	vector<Vector3d>::iterator i=positions.begin();
	Location lb=loc_type.sub("model.beads.bead[]");
	for(;i!=positions.end();i++){
		lb.next();
		out_udf->put(lb.sub("position.x"),i->x);
		out_udf->put(lb.sub("position.y"),i->y);
		out_udf->put(lb.sub("position.z"),i->z);
	}
	ResistanceTensor rt=hba->getResistanceTensor();
	output_common(calcTime,translate,rotate,rt);
	out_udf->write();
	delete hba;
	cout << endl;
}

void ParticleAnalyzer::analyze_BEM(){
	Mesh mesh;
	Location bd(loc_type.sub("model.BEM.vertex[]"));
	mesh.vertex.reserve(in_udf->size(bd));
	map<int,int> vertex_id_index;
	for(int i=0;i<in_udf->size(bd);i++){
		bd.next();
		vertex_id_index[in_udf->i(bd.sub("id"))]=i;
		Vector3d r(in_udf->d(bd.sub("position.x")),
		           in_udf->d(bd.sub("position.y")),
		           in_udf->d(bd.sub("position.z")));
		mesh.vertex.push_back(new Vertex(r));
	}
	cout << "\tnumber of verticies: " << mesh.vertex.size() << endl;
	Location loc_face(loc_type.sub("model.BEM.face[]"));
	mesh.face.reserve(in_udf->size(loc_face));
	for(int j=0;j<in_udf->size(loc_face);j++){
		loc_face.next();
		Location loc(loc_face.sub("vertex[]"));
		loc.next();Vertex* v0=mesh.vertex[vertex_id_index[in_udf->i(loc)]];
		loc.next();Vertex* v1=mesh.vertex[vertex_id_index[in_udf->i(loc)]];
		loc.next();Vertex* v2=mesh.vertex[vertex_id_index[in_udf->i(loc)]];
		mesh.face.push_back(new Face(v0,v1,v2));
	}
	cout << "\tnumber of faces: " << mesh.face.size() << endl;
	string c=in_udf->s(loc_type.sub("center"));
	cout << "\tcenter of tensor: " << c << endl;
	string d=in_udf->s(loc_type.sub("tensor_diagonalization"));
	cout << "\ttensor_diagonalization: " << d << endl;
	double viscosity=in_udf->d("fluid.viscosity");
	BoundaryElementModel* hba=new BoundaryElementModel(viscosity,&mesh,c,d);
	Vector3d translate;Tensor3x3 rotate;
	Timer timer;
	hba->update(translate,rotate);
	double calcTime=timer.get();
	vector<Vertex*>::iterator i=mesh.vertex.begin();
	Location lb=loc_type.sub("model.BEM.vertex[]");
	for(;i!=mesh.vertex.end();i++){
		lb.next();
		Vector3d r=(*i)->position;
		out_udf->put(lb.sub("position.x"),r.x);
		out_udf->put(lb.sub("position.y"),r.y);
		out_udf->put(lb.sub("position.z"),r.z);
	}
	ResistanceTensor rt=hba->getResistanceTensor();
	output_common(calcTime,translate,rotate,rt);
	out_udf->write();
	//delete objects
	for(i=mesh.vertex.begin();i!=mesh.vertex.end();i++){
		delete *i;
	}
	vector<Face*>::iterator j=mesh.face.begin();
	for(;j!=mesh.face.end();j++){
		delete *j;
	}
	mesh.vertex.clear();mesh.face.clear();
	delete hba;
	cout << endl;
}

void ParticleAnalyzer::analyze_superposition(){
	//superposition model input from UDF file
	vector<AxisymmetricResistanceTensor*> tensor;
	Location el(loc_type.sub("model.superposition.element[]"));
	for(int e=0;e<in_udf->size(el);e++){
		el.next();
		string elementName=in_udf->s(el.sub("element"));
		if(elementName=="sphere"){
			double a=in_udf->d(el.sub("sphere.radius"));
			tensor.push_back(new Sphere(a));
		}else if(elementName=="prolate_spheroid"){
			double ux=in_udf->d(el.sub("prolate_spheroid.director.x"));
			double uy=in_udf->d(el.sub("prolate_spheroid.director.y"));
			double uz=in_udf->d(el.sub("prolate_spheroid.director.z"));
			double la=in_udf->d(el.sub("prolate_spheroid.long_axis"));
			double sa=in_udf->d(el.sub("prolate_spheroid.short_axis"));
			tensor.push_back(new ProlateSpheroid(Vector3d(ux,uy,uz),la,sa));
		}else if(elementName=="oblate_spheroid"){
			double ux=in_udf->d(el.sub("oblate_spheroid.director.x"));
			double uy=in_udf->d(el.sub("oblate_spheroid.director.y"));
			double uz=in_udf->d(el.sub("oblate_spheroid.director.z"));
			double la=in_udf->d(el.sub("oblate_spheroid.long_axis"));
			double sa=in_udf->d(el.sub("oblate_spheroid.short_axis"));
			tensor.push_back(new OblateSpheroid(Vector3d(ux,uy,uz),la,sa));
		}else if(elementName=="disk"){
			double ux=in_udf->d(el.sub("disk.director.x"));
			double uy=in_udf->d(el.sub("disk.director.y"));
			double uz=in_udf->d(el.sub("disk.director.z"));
			double a=in_udf->d(el.sub("disk.radius"));
			tensor.push_back(new Disk(Vector3d(ux,uy,uz),a));
		}else if(elementName=="needle"){
			double ux=in_udf->d(el.sub("needle.director.x"));
			double uy=in_udf->d(el.sub("needle.director.y"));
			double uz=in_udf->d(el.sub("needle.director.z"));
			double la=in_udf->d(el.sub("needle.long_axis"));
			double sa=in_udf->d(el.sub("needle.short_axis"));
			tensor.push_back(new Needle(Vector3d(ux,uy,uz),la,sa));
		}else if(elementName=="point"){
			tensor.push_back(new Point);
		}else if(elementName=="line"){
			double ux=in_udf->d(el.sub("line.director.x"));
			double uy=in_udf->d(el.sub("line.director.y"));
			double uz=in_udf->d(el.sub("line.director.z"));
			double h=in_udf->d(el.sub("line.half_length"));
			tensor.push_back(new Line(Vector3d(ux,uy,uz),h));
		}
		Vector3d r(in_udf->d(el.sub("position.x")),in_udf->d(el.sub("position.y")),in_udf->d(el.sub("position.z")));
		(*(tensor.end()-1))->translate(r);
	}
	//calculation of superposition model
	string c=in_udf->s(loc_type.sub("center"));
	cout << "\tcenter of tensor: " << c << endl;
	string t=in_udf->s(loc_type.sub("tensor_diagonalization"));
	cout << "\ttensor_diagonalization: " << t << endl;
	double viscosity=in_udf->d("fluid.viscosity");
	SuperpositionModel* sa=new SuperpositionModel(viscosity,tensor,c,t);
	Vector3d translate;Tensor3x3 rotate;
	Timer timer;
	sa->update(translate,rotate);
	double calcTime=timer.get();
	//output to analysis data to UDF file
	el.seek(loc_type.sub("model.superposition.element[]").str());
	for(int e=0;e<in_udf->size(el);e++){
		el.next();
		Vector3d r(in_udf->d(el.sub("position.x")),in_udf->d(el.sub("position.y")),in_udf->d(el.sub("position.z")));
		string elementName=out_udf->s(el.sub("element"));
		r+=translate;r=convert(rotate,r);
		out_udf->put(el.sub("position.x"),r.x);out_udf->put(el.sub("position.y"),r.y);out_udf->put(el.sub("position.z"),r.z);
		if(elementName=="prolate_spheroid" || elementName=="oblate_spheroid" ||
		   elementName=="disk" || elementName=="needle" || elementName=="line"){
			Location loc_dir(el.sub(elementName+".director"));
			Vector3d u(in_udf->d(loc_dir.sub("x")),in_udf->d(loc_dir.sub("y")),in_udf->d(loc_dir.sub("z")));
			u=convert(rotate,u);
			out_udf->put(loc_dir.sub("x"),u.x);out_udf->put(loc_dir.sub("y"),u.y);out_udf->put(loc_dir.sub("z"),u.z);
		}
	}
	ResistanceTensor rt=sa->getResistanceTensor();
	output_common(calcTime,translate,rotate,rt);
	out_udf->write();
	delete sa;
	vector<AxisymmetricResistanceTensor*>::iterator n=tensor.begin();
	for(;n!=tensor.end();n++)delete *n;
	cout << endl;
}

void ParticleAnalyzer::output_common(
	double time,Vector3d& translate,Tensor3x3& rotate,ResistanceTensor& rt){
	//�ϊ������e���\���A���s�ړ������x�N�g�����o��
	Location infoTensor(loc_type.sub("analysis_information.convert_tensor"));
	out_udf->put(infoTensor.sub("x.x"),rotate.x.x);
	out_udf->put(infoTensor.sub("x.y"),rotate.x.y);
	out_udf->put(infoTensor.sub("x.z"),rotate.x.z);
	out_udf->put(infoTensor.sub("y.x"),rotate.y.x);
	out_udf->put(infoTensor.sub("y.y"),rotate.y.y);
	out_udf->put(infoTensor.sub("y.z"),rotate.y.z);
	out_udf->put(infoTensor.sub("z.x"),rotate.z.x);
	out_udf->put(infoTensor.sub("z.y"),rotate.z.y);
	out_udf->put(infoTensor.sub("z.z"),rotate.z.z);
	Location infoVector(loc_type.sub("analysis_information.shift_vector"));
	out_udf->put(infoVector.sub("x"),translate.x);
	out_udf->put(infoVector.sub("y"),translate.y);
	out_udf->put(infoVector.sub("z"),translate.z);
	//�v�Z���Ԃ̏o��
	out_udf->put(loc_type.sub("analysis_information.cpu_time"),time);
	//��R�e���\���̏o��
	Tensor3x3 A=rt.get_A();
	Location tt=loc_type.sub("resistance_tensor.A");
	out_udf->put(tt.sub("xx"),A(0,0));
	out_udf->put(tt.sub("yy"),A(1,1));
	out_udf->put(tt.sub("zz"),A(2,2));
	out_udf->put(tt.sub("yx"),A(1,0));
	out_udf->put(tt.sub("zx"),A(2,0));
	out_udf->put(tt.sub("zy"),A(2,1));
	Tensor3x3 B=rt.get_B();
	Location tr=loc_type.sub("resistance_tensor.B");
	out_udf->put(tr.sub("xx"),B(0,0));
	out_udf->put(tr.sub("yy"),B(1,1));
	out_udf->put(tr.sub("zz"),B(2,2));
	out_udf->put(tr.sub("yx"),B(1,0));
	out_udf->put(tr.sub("zx"),B(2,0));
	out_udf->put(tr.sub("zy"),B(2,1));
	Tensor3x3 C=rt.get_C();
	Location rr=loc_type.sub("resistance_tensor.C");
	out_udf->put(rr.sub("xx"),C(0,0));
	out_udf->put(rr.sub("yy"),C(1,1));
	out_udf->put(rr.sub("zz"),C(2,2));
	out_udf->put(rr.sub("yx"),C(1,0));
	out_udf->put(rr.sub("zx"),C(2,0));
	out_udf->put(rr.sub("zy"),C(2,1));
	Tensor3x3x3 g=rt.get_G();
	Location te=loc_type.sub("resistance_tensor.G");
	out_udf->put(te.sub("xx.x"),g(0,0,0));
	out_udf->put(te.sub("xx.y"),g(0,0,1));
	out_udf->put(te.sub("xx.z"),g(0,0,2));//G.xx
	out_udf->put(te.sub("yy.x"),g(1,1,0));
	out_udf->put(te.sub("yy.y"),g(1,1,1));
	out_udf->put(te.sub("yy.z"),g(1,1,2));//G.yy
	out_udf->put(te.sub("zz.x"),g(2,2,0));
	out_udf->put(te.sub("zz.y"),g(2,2,1));
	out_udf->put(te.sub("zz.z"),g(2,2,2));//G.zz
	out_udf->put(te.sub("yx.x"),g(1,0,0));
	out_udf->put(te.sub("yx.y"),g(1,0,1));
	out_udf->put(te.sub("yx.z"),g(1,0,2));//G.yx
	out_udf->put(te.sub("zx.x"),g(2,0,0));
	out_udf->put(te.sub("zx.y"),g(2,0,1));
	out_udf->put(te.sub("zx.z"),g(2,0,2));//G.zx
	out_udf->put(te.sub("zy.x"),g(2,1,0));
	out_udf->put(te.sub("zy.y"),g(2,1,1));
	out_udf->put(te.sub("zy.z"),g(2,1,2));//G.zy
	Tensor3x3x3 h=rt.get_H();
	Location re=loc_type.sub("resistance_tensor.H");
	out_udf->put(re.sub("xx.x"),h(0,0,0));
	out_udf->put(re.sub("xx.y"),h(0,0,1));
	out_udf->put(re.sub("xx.z"),h(0,0,2));//H.xx
	out_udf->put(re.sub("yy.x"),h(1,1,0));
	out_udf->put(re.sub("yy.y"),h(1,1,1));
	out_udf->put(re.sub("yy.z"),h(1,1,2));//H.yy
	out_udf->put(re.sub("zz.x"),h(2,2,0));
	out_udf->put(re.sub("zz.y"),h(2,2,1));
	out_udf->put(re.sub("zz.z"),h(2,2,2));//H.zz
	out_udf->put(re.sub("yx.x"),h(1,0,0));
	out_udf->put(re.sub("yx.y"),h(1,0,1));
	out_udf->put(re.sub("yx.z"),h(1,0,2));//H.yx
	out_udf->put(re.sub("zx.x"),h(2,0,0));
	out_udf->put(re.sub("zx.y"),h(2,0,1));
	out_udf->put(re.sub("zx.z"),h(2,0,2));//H.zx
	out_udf->put(re.sub("zy.x"),h(2,1,0));
	out_udf->put(re.sub("zy.y"),h(2,1,1));
	out_udf->put(re.sub("zy.z"),h(2,1,2));//H.zy
	Tensor3x3x3x3 k=rt.get_K();
	//begin of traceless-like operation
	Tensor3x3 deltaK=dot2(tensor_utility::delta(),k),Kdelta=dot2(k,tensor_utility::delta());
	double deltaKdelta=dot2(dot2(tensor_utility::delta(),k),tensor_utility::delta());
	for(int i0=0;i0<3;i0++)
		for(int i1=0;i1<3;i1++)
			for(int i2=0;i2<3;i2++)
				for(int i3=0;i3<3;i3++)
					k(i0,i1,i2,i3)+=-1./3.*(tensor_utility::delta(i0,i1)*deltaK(i2,i3)+Kdelta(i0,i1)*tensor_utility::delta(i2,i3))
						            +1./9.*tensor_utility::delta(i0,i1)*deltaKdelta*tensor_utility::delta(i2,i3);
	//end of traceless-like operation
	Location ke=loc_type.sub("resistance_tensor.K");
	out_udf->put(ke.sub("xx.xx"),k(0,0,0,0));//K.xx
	out_udf->put(ke.sub("yy.xx"),k(1,1,0,0));
	out_udf->put(ke.sub("yy.yy"),k(1,1,1,1));//K.yy
	out_udf->put(ke.sub("zz.xx"),k(2,2,0,0));
	out_udf->put(ke.sub("zz.yy"),k(2,2,1,1));
	out_udf->put(ke.sub("zz.zz"),k(2,2,2,2));//K.zz
	out_udf->put(ke.sub("yx.xx"),k(1,0,0,0));
	out_udf->put(ke.sub("yx.yy"),k(1,0,1,1));
	out_udf->put(ke.sub("yx.zz"),k(1,0,2,2));
	out_udf->put(ke.sub("yx.yx"),k(1,0,1,0));//K.yx
	out_udf->put(ke.sub("zx.xx"),k(2,0,0,0));
	out_udf->put(ke.sub("zx.yy"),k(2,0,1,1));
	out_udf->put(ke.sub("zx.zz"),k(2,0,2,2));
	out_udf->put(ke.sub("zx.yx"),k(2,0,1,0));
	out_udf->put(ke.sub("zx.zx"),k(2,0,2,0));//K.xz
	out_udf->put(ke.sub("zy.xx"),k(2,1,0,0));
	out_udf->put(ke.sub("zy.yy"),k(2,1,1,1));
	out_udf->put(ke.sub("zy.zz"),k(2,1,2,2));
	out_udf->put(ke.sub("zy.yx"),k(2,1,1,0));
	out_udf->put(ke.sub("zy.zx"),k(2,1,2,0));
	out_udf->put(ke.sub("zy.zy"),k(2,1,2,1));//K.yz
	//dipole, susceptibility convert
	Vector3d d(in_udf->d(loc_type.sub("dipole.x")),
	           in_udf->d(loc_type.sub("dipole.y")),
	           in_udf->d(loc_type.sub("dipole.z")));
	d=rotate*d;
	out_udf->put(loc_type.sub("dipole.x"),d.x);
	out_udf->put(loc_type.sub("dipole.y"),d.y);
	out_udf->put(loc_type.sub("dipole.z"),d.z);
	//
	Tensor3x3 chi;
	chi.x.x=in_udf->d(loc_type.sub("susceptibility.xx"));
	chi.y.y=in_udf->d(loc_type.sub("susceptibility.yy"));
	chi.z.z=in_udf->d(loc_type.sub("susceptibility.zz"));
	chi.x.y=chi.y.x=in_udf->d(loc_type.sub("susceptibility.yx"));
	chi.x.z=chi.z.x=in_udf->d(loc_type.sub("susceptibility.zx"));
	chi.y.z=chi.z.y=in_udf->d(loc_type.sub("susceptibility.zy"));
	chi=rotate*chi*(rotate.transpose());
	out_udf->put(loc_type.sub("susceptibility.xx"),chi.x.x);
	out_udf->put(loc_type.sub("susceptibility.yy"),chi.y.y);
	out_udf->put(loc_type.sub("susceptibility.zz"),chi.z.z);
	out_udf->put(loc_type.sub("susceptibility.yx"),chi.y.x);
	out_udf->put(loc_type.sub("susceptibility.zx"),chi.z.x);
	out_udf->put(loc_type.sub("susceptibility.zy"),chi.z.y);
}
