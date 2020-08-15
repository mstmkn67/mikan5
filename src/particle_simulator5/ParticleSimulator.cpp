#include "ParticleSimulator.h"

//Force
#include "ConstantForce.h"
#include "ElectricMagneticTorque.h"
//VelocityField
#include "../common/velocityField/LinearFlow.h"
#include "../common/velocityField/PoiseuilleFlow.h"
#include "../common/velocityField/UdfFlow.h"
#include "../common/velocityField/DynamicFlow.h"
//ElectricMagneticField
#include "../common/electricMagneticField/StaticField.h"
#include "../common/electricMagneticField/EllipticRotationField.h"
#include "../common/electricMagneticField/SawtoothField.h"

ParticleSimulator::ParticleSimulator(UDFManager* inUdf,UDFManager* outUdf)
:inudf(inUdf),outudf(outUdf){
	timer.start();
	double lx=inudf->d("simulation.system_size.size.x");
	double ly=inudf->d("simulation.system_size.size.y");
	double lz=inudf->d("simulation.system_size.size.z");
	volume=lx*ly*lz;
	assign_particle();
	assign_force();
	assign_electricField();
	assign_electric_torque();
	assign_velocityField();
}

ParticleSimulator::~ParticleSimulator(){
	delete velocityField;
	vector<Force*>::iterator f=force.begin();
	for(;f!=force.end();f++)delete *f;
	map<int,vector<Particle*> >::iterator p=particle.begin();
	for(;p!=particle.end();p++){
		vector<Particle*>::iterator q=p->second.begin();
		for(;q!=p->second.end();q++)delete *q;
	}
	map<int,MobilityTensorOfBrownianSimulation*>::iterator t=tensor.begin();
	for(;t!=tensor.end();t++)delete t->second;
}

void ParticleSimulator::assign_tensor(int id){
	Location pt=inudf->getLocation(string("ParticleType"),id);
	pt.up();
	cout << "assign mobility tensor of " << pt << endl;
	Tensor3x3 A,B,C;Tensor3x3x3 G,H;Tensor3x3x3x3 K;
	Location fa=pt.sub("resistance_tensor.A");
	A.x.x=inudf->d(fa.sub("xx"));A.y.y=inudf->d(fa.sub("yy"));A.z.z=inudf->d(fa.sub("zz"));
	A.y.x=A.x.y=inudf->d(fa.sub("xy"));A.z.x=A.x.z=inudf->d(fa.sub("xz"));A.z.y=A.y.z=inudf->d(fa.sub("yz"));
	Location fb=pt.sub("resistance_tensor.B");
	B.x.x=inudf->d(fb.sub("xx"));B.y.y=inudf->d(fb.sub("yy"));B.z.z=inudf->d(fb.sub("zz"));
	B.y.x=B.x.y=inudf->d(fb.sub("xy"));B.z.x=B.x.z=inudf->d(fb.sub("xz"));B.z.y=B.y.z=inudf->d(fb.sub("yz"));
	Location fc=pt.sub("resistance_tensor.C");
	C.x.x=inudf->d(fc.sub("xx"));C.y.y=inudf->d(fc.sub("yy"));C.z.z=inudf->d(fc.sub("zz"));
	C.y.x=C.x.y=inudf->d(fc.sub("xy"));C.z.x=C.x.z=inudf->d(fc.sub("xz"));C.z.y=C.y.z=inudf->d(fc.sub("yz"));
	Location fg=pt.sub("resistance_tensor.G");
	G.x.x.x=inudf->d(fg.sub("xx.x"));
	G.x.x.y=inudf->d(fg.sub("xx.y"));
	G.x.x.z=inudf->d(fg.sub("xx.z"));//G.xx
	G.y.y.x=inudf->d(fg.sub("yy.x"));
	G.y.y.y=inudf->d(fg.sub("yy.y"));
	G.y.y.z=inudf->d(fg.sub("yy.z"));//G.yy
	G.z.z.x=inudf->d(fg.sub("zz.x"));
	G.z.z.y=inudf->d(fg.sub("zz.y"));
	G.z.z.z=inudf->d(fg.sub("zz.z"));//G.zz
	G.y.x.x=G.x.y.x=inudf->d(fg.sub("yx.x"));
	G.y.x.y=G.x.y.y=inudf->d(fg.sub("yx.y"));
	G.y.x.z=G.x.y.z=inudf->d(fg.sub("yx.z"));//G.yx
	G.z.x.x=G.x.z.x=inudf->d(fg.sub("zx.x"));
	G.z.x.y=G.x.z.y=inudf->d(fg.sub("zx.y"));
	G.z.x.z=G.x.z.z=inudf->d(fg.sub("zx.z"));//G.zx
	G.z.y.x=G.y.z.x=inudf->d(fg.sub("zy.x"));
	G.z.y.y=G.y.z.y=inudf->d(fg.sub("zy.y"));
	G.z.y.z=G.y.z.z=inudf->d(fg.sub("zy.z"));//G.zy
	Location fh=pt.sub("resistance_tensor.H");
	H.x.x.x=inudf->d(fh.sub("xx.x"));
	H.x.x.y=inudf->d(fh.sub("xx.y"));
	H.x.x.z=inudf->d(fh.sub("xx.z"));//G.xx
	H.y.y.x=inudf->d(fh.sub("yy.x"));
	H.y.y.y=inudf->d(fh.sub("yy.y"));
	H.y.y.z=inudf->d(fh.sub("yy.z"));//G.yy
	H.z.z.x=inudf->d(fh.sub("zz.x"));
	H.z.z.y=inudf->d(fh.sub("zz.y"));
	H.z.z.z=inudf->d(fh.sub("zz.z"));//G.zz
	H.y.x.x=H.x.y.x=inudf->d(fh.sub("yx.x"));
	H.y.x.y=H.x.y.y=inudf->d(fh.sub("yx.y"));
	H.y.x.z=H.x.y.z=inudf->d(fh.sub("yx.z"));//G.yx
	H.z.x.x=H.x.z.x=inudf->d(fh.sub("zx.x"));
	H.z.x.y=H.x.z.y=inudf->d(fh.sub("zx.y"));
	H.z.x.z=H.x.z.z=inudf->d(fh.sub("zx.z"));//G.zx
	H.z.y.x=H.y.z.x=inudf->d(fh.sub("zy.x"));
	H.z.y.y=H.y.z.y=inudf->d(fh.sub("zy.y"));
	H.z.y.z=H.y.z.z=inudf->d(fh.sub("zy.z"));//G.zy
	Location fk=pt.sub("resistance_tensor.K");
	K.x.x.x.x=inudf->d(fk.sub("xx.xx"));//K.xx
	K.y.y.x.x=K.x.x.y.y=inudf->d(fk.sub("yy.xx"));
	K.y.y.y.y=inudf->d(fk.sub("yy.yy"));//K.yy
	K.z.z.x.x=K.x.x.z.z=inudf->d(fk.sub("zz.xx"));
	K.z.z.y.y=K.y.y.z.z=inudf->d(fk.sub("zz.yy"));	
	K.z.z.z.z=inudf->d(fk.sub("zz.zz"));//K.zz
	K.y.x.x.x=K.x.x.y.x=K.x.y.x.x=K.x.x.x.y=inudf->d(fk.sub("yx.xx"));
	K.y.x.y.y=K.y.y.y.x=K.x.y.y.y=K.y.y.x.y=inudf->d(fk.sub("yx.yy"));
	K.y.x.z.z=K.z.z.y.x=K.x.y.z.z=K.z.z.x.y=inudf->d(fk.sub("yx.zz"));
	K.y.x.y.x=K.x.y.y.x=K.y.x.x.y=K.x.y.x.y=inudf->d(fk.sub("yx.yx"));//K.yx
	K.z.x.x.x=K.x.x.z.x=K.x.z.x.x=K.x.x.x.z=inudf->d(fk.sub("zx.xx"));
	K.z.x.y.y=K.y.y.z.x=K.x.z.y.y=K.y.y.x.z=inudf->d(fk.sub("zx.yy"));
	K.z.x.z.z=K.z.z.z.x=K.x.z.z.z=K.z.z.x.z=inudf->d(fk.sub("zx.zz"));
	K.z.x.y.x=K.x.z.y.x=K.z.x.x.y=K.x.z.x.y
		=K.y.x.z.x=K.x.y.z.x=K.y.x.x.z=K.x.y.x.z=inudf->d(fk.sub("zx.yx"));
	K.z.x.z.x=K.x.z.z.x=K.z.x.x.z=K.x.z.x.z=inudf->d(fk.sub("zx.zx"));//K.zx
	K.z.y.x.x=K.x.x.z.y=K.y.z.x.x=K.x.x.y.z=inudf->d(fk.sub("zy.xx"));
	K.z.y.y.y=K.y.y.z.y=K.y.z.y.y=K.y.y.y.z=inudf->d(fk.sub("zy.yy"));
	K.z.y.z.z=K.z.z.z.y=K.y.z.z.z=K.z.z.y.z=inudf->d(fk.sub("zy.zz"));
	K.z.y.y.x=K.y.z.y.x=K.z.y.x.y=K.y.z.x.y
		=K.y.x.z.y=K.x.y.z.y=K.y.x.y.z=K.x.y.y.z=inudf->d(fk.sub("zy.yx"));
	K.z.y.z.x=K.y.z.z.x=K.z.y.x.z=K.y.z.x.z
		=K.z.x.z.y=K.x.z.z.y=K.z.x.y.z=K.x.z.y.z=inudf->d(fk.sub("zy.zx"));
	K.z.y.z.y=K.y.z.z.y=K.z.y.y.z=K.y.z.y.z=inudf->d(fk.sub("zy.zy"));//K.zy
	double dt=inudf->doubleValue("simulation.time.dt");
	double kBT=inudf->d("simulation.Brownian_motion.true.kBT");
	tensor[id]=new MobilityTensorOfBrownianSimulation(kBT,dt,ResistanceTensor(A,B,C,G,H,K));
}

void ParticleSimulator::assign_particle(){
	double dt=inudf->doubleValue("simulation.time.dt");
	Particle::setDt(dt);
	string b(inudf->stringValue("simulation.Brownian_motion.Brownian_flag"));
	if(b=="true"){
		Particle::setBrownianFlag(true);
		//double kBT=inudf->d("simulation.Brownian_motion.true.kBT");
	}
	Location p("simulation.particle[]");
	for(int i=0;i<inudf->size(p);i++){
		p.next();
		int type_id=inudf->intValue(p.sub("particle_type"));
		Vector3d r(inudf->d(p.sub("initial_position.x")),
               inudf->d(p.sub("initial_position.y")),
               inudf->d(p.sub("initial_position.z")));
		EulerAngle e(inudf->d(p.sub("initial_quaternion.xi")),
                 inudf->d(p.sub("initial_quaternion.eta")),
		             inudf->d(p.sub("initial_quaternion.zeta")),
                 inudf->d(p.sub("initial_quaternion.chi")));
		double s=e.xi*e.xi+e.eta*e.eta+e.zeta*e.zeta+e.chi*e.chi;
		if(s<0.99 || s>1.01){
			cout << "[0,0,0,1] is set in simulation.particle[" << i << "].initial_quaternion " << endl;
			e.set(0.0,0.0,0.0,1.0);
		}
		if(tensor.find(type_id)==tensor.end()){
			assign_tensor(type_id);//when type_id's mobility tensor isn't assigned
		}
		Particle* pp=new Particle(i,tensor[type_id],r,e);
		particle[type_id].push_back(pp);
	}
}

void ParticleSimulator::assign_force(){
	Location fl("simulation.external_force[]");
	for(int i=0;i<inudf->size(fl);i++){
		fl.next();
		vector<Particle*> pt;
		Location idl=fl.sub("particle_type[]");
		int n=inudf->size(idl);
		for(int j=0;j<n;j++){
			idl.next();
			int k=inudf->i(idl);
			copy(particle[k].begin(),particle[k].end(),back_inserter(pt));
		}
		Vector3d f(inudf->d(fl.sub("force.x")),inudf->d(fl.sub("force.y")),inudf->d(fl.sub("force.z")));
		Vector3d t(inudf->d(fl.sub("torque.x")),inudf->d(fl.sub("torque.y")), inudf->d(fl.sub("torque.z")));
		Vector3d p(inudf->d(fl.sub("acting_point.x")),inudf->d(fl.sub("acting_point.y")),inudf->d(fl.sub("acting_point.z")));
		bool flag=true;
		if(inudf->s(fl.sub("coordination_system"))=="laboratory_system")flag=false;
		force.push_back(new ConstantForce(pt,f,t,p,flag));
	}
}

void ParticleSimulator::assign_electric_torque(){
	string name=inudf->stringValue("simulation.external_field.electric_field.field_type");
	if(name=="false")return;
	map<int,vector<Particle*> >::iterator p;
	for(p=particle.begin();p!=particle.end();p++){
		Location loc=inudf->getLocation("ParticleType",p->first);
		loc.up();
		double dx=inudf->d(loc.sub("dipole.x"));
		double dy=inudf->d(loc.sub("dipole.y"));
		double dz=inudf->d(loc.sub("dipole.z"));
		force.push_back(new DipoleTorque(p->second,Vector3d(dx,dy,dz),electricField));
		double alpha=inudf->d(loc.sub("susceptibility.electric_alpha"));
		double xx=inudf->d(loc.sub("susceptibility.xx"));
		double yy=inudf->d(loc.sub("susceptibility.yy"));
		double zz=inudf->d(loc.sub("susceptibility.zz"));
		double yx=inudf->d(loc.sub("susceptibility.yx"));
		double zx=inudf->d(loc.sub("susceptibility.zx"));
		double zy=inudf->d(loc.sub("susceptibility.zy"));
		force.push_back(new SusceptibilityTorque(p->second,alpha,Tensor3x3(xx,yx,zx, yx,yy,zy, zx,zy,zz),electricField));
	}
}

void ParticleSimulator::assign_velocityField(){
	Particle::setVelocityField(velocityField);
	string name=inudf->stringValue("simulation.external_field.velocity_field.velocity_field_type");
	cout << "velocity field: " ;
	if(name=="simple_shear_flow"){
		cout << name << endl;
		double gamma=inudf->doubleValue("simulation.external_field.velocity_field.simple_shear_flow.shear_rate");
		velocityField=new SimpleShearFlow(gamma);
	}else if(name=="dynamic_shear_flow"){
		double strain=inudf->doubleValue("simulation.external_field.velocity_field.dynamic_shear_flow.strain");
		double omega=inudf->doubleValue("simulation.external_field.velocity_field.dynamic_shear_flow.angular_frequency");
		velocityField=new DynamicShearFlow(strain,omega);
	}else if(name=="elongational_flow"){
		cout << name << endl;
		double s=inudf->d("simulation.external_field.velocity_field.elongational_flow.strain_rate");
		double k=inudf->d("simulation.external_field.velocity_field.elongational_flow.k");
		velocityField=new ElongationalFlow(k,s);
	}else if(name=="linear_flow"){
		cout << name << endl;
		Location loc("simulation.external_field.velocity_field.linear_flow");
		Vector3d a(inudf->d(loc.sub("A.x")),inudf->d(loc.sub("A.y")),inudf->d(loc.sub("A.z")));
		Tensor3x3 b(inudf->d(loc.sub("B.x.x")),inudf->d(loc.sub("B.x.y")),inudf->d(loc.sub("B.x.z")),
		            inudf->d(loc.sub("B.y.x")),inudf->d(loc.sub("B.y.y")),inudf->d(loc.sub("B.y.z")),
		            inudf->d(loc.sub("B.z.x")),inudf->d(loc.sub("B.z.y")),inudf->d(loc.sub("B.z.z")));
		velocityField=new LinearFlow(a,b);
	}else if(name=="plane_Poiseuille_flow"){
		cout << name << endl;
		Location loc("simulation.external_field.velocity_field.plane_Poiseuille_flow");
		double h=inudf->d(loc.sub("half_length_between_plates"));
		double v=inudf->d(loc.sub("max_speed"));
		velocityField=new PlanePoiseuilleFlow(h,v);
	}else if(name=="cylinder_Poiseuille_flow"){
		cout << name << endl;
		Location loc("simulation.external_field.velocity_field.cylinder_Poiseuille_flow");
		double r=inudf->d(loc.sub("radius"));
		double v=inudf->d(loc.sub("max_speed"));
	}else if(name=="static_udf_flow2d"){
		cout << name << endl;
		string udfName=inudf->s("simulation.external_field.velocity_field.static_udf_flow2d.file_name");
		cout << "\tudf file: " << udfName << endl;
		double l=inudf->d("unit_parameter.length");
		double t=inudf->d("unit_parameter.time");
		double e=inudf->d("unit_parameter.solvent_viscosity");
		velocityField=new UdfFlow2d(udfName,l,t,e);
	}else if(name=="static_udf_flow3d"){
		cout << name << endl;
		string udfName=inudf->s("simulation.external_field.velocity_field.static_udf_flow3d.file_name");
		cout << "\tudf file: " << udfName << endl;
		double l=inudf->d("unit_parameter.length");
		double t=inudf->d("unit_parameter.time");
		double e=inudf->d("unit_parameter.solvent_viscosity");
		velocityField=new UdfFlow3d(udfName,l,t,e);
	}else if(name=="dynamic_udf_flow2d"){
		cout << name << endl;
		string udfName=inudf->s("simulation.external_field.velocity_field.dynamic_udf_flow2d.file_name");
		cout << "\tudf file: " << udfName << endl;
		double l=inudf->d("unit_parameter.length");
		double t=inudf->d("unit_parameter.time");
		double e=inudf->d("unit_parameter.solvent_viscosity");
		double dt=inudf->d("simulation.time.dt");
		velocityField=new DynamicUdfFlow2d(udfName,l,t,e,dt);
	}else if(name=="dynamic_udf_flow3d"){
		cout << name << endl;
		string udfName=inudf->s("simulation.external_field.velocity_field.dynamic_udf_flow3d.file_name");
		cout << "\tudf file: " << udfName << endl;
		double l=inudf->d("unit_parameter.length");
		double t=inudf->d("unit_parameter.time");
		double e=inudf->d("unit_parameter.solvent_viscosity");
		double dt=inudf->d("simulation.time.dt");
		velocityField=new DynamicUdfFlow3d(udfName,l,t,e,dt);
	}else{
		cout << "non field" << endl;
		velocityField=new NonFlow;
	}
	Particle::setVelocityField(velocityField);
}

void ParticleSimulator::assign_electricField(){
	string name=inudf->stringValue("simulation.external_field.electric_field.field_type");
	cout << "electric field: ";
	if(name=="static_field"){
		cout << name << endl;
		Vector3d E(inudf->d("simulation.external_field.electric_field.constant_field.E.x"),
		           inudf->d("simulation.external_field.electric_field.constant_field.E.y"),
		           inudf->d("simulation.external_field.electric_field.constant_field.E.z"));
		electricField=new StaticField(E);
	}else if(name=="elliptic_rotation_field"){
		cout << name << endl;
		double w1=inudf->d("simulation.external_field.electric_field.elliptic_rotation_field.angular_velocity1");
		Vector3d E1(inudf->d("simulation.external_field.electric_field.elliptic_rotation_field.E1.x"),
		            inudf->d("simulation.external_field.electric_field.elliptic_rotation_field.E1.y"),
		            inudf->d("simulation.external_field.electric_field.elliptic_rotation_field.E1.z"));
		double w2=inudf->d("simulation.external_field.electric_field.elliptic_rotation_field.angular_velocity2");
		Vector3d E2(inudf->d("simulation.external_field.electric_field.elliptic_rotation_field.E2.x"),
		            inudf->d("simulation.external_field.electric_field.elliptic_rotation_field.E2.y"),
		            inudf->d("simulation.external_field.electric_field.elliptic_rotation_field.E2.z"));
		double dt=inudf->d("simulation.time.dt");
		electricField=new EllipticRotationField(w1,E1,w2,E2,dt);
	}else if(name=="sawtooth_field"){
		cout << name << endl;
		double umax=inudf->d("simulation.external_field.electric_field.sawtooth_field.umax");
		double alpha=inudf->d("simulation.external_field.electric_field.sawtooth_field.alpha");
		double length=inudf->d("simulation.external_field.electric_field.sawtooth_field.length");
		double tau=inudf->d("simulation.external_field.electric_field.sawtooth_field.swith_time");
		double dt=inudf->d("simulation.time.dt");
		electricField=new SawtoothField(umax,alpha,length,tau,dt);
	}else{
		cout << "non field" << endl;
		electricField=new StaticField(Vector3d(0.0,0.0,0.0));
	}
	Particle::setElectricField(electricField);
}

void ParticleSimulator::run(){
	string name=inudf->stringValue("simulation.integrator");
	if(name=="Euler")run_Euler();
	else if(name=="4th_order_Runge_Kutta")run_4thOrderRungeKutta();
	else{
		cout << "simulation.integrator is not selected." << endl;
	}
}

void ParticleSimulator::run_Euler(){
	cout << "simulation starts using Euler method" << endl;
	time=0.0;
	velocityField->initial();
	int total=inudf->i("simulation.time.simulation_steps");
	int record=inudf->i("simulation.time.report_steps");
	double dt=inudf->d("simulation.time.dt");
	report();
	velocityField->time=&time;
	electricField->time=&time;
	for(int itime=1;itime<=total;itime++){
		time=dt*itime;
		velocityField->update();//���x��X�V
		map<int,vector<Particle*> >::iterator p=particle.begin();
		for(;p!=particle.end();p++){
			vector<Particle*>::iterator r=p->second.begin();
			for(;r!=p->second.end();r++){
				(*r)->forceAndTorqueReset();//���q�ɂ�����͂̃��Z�b�g
				(*r)->beadsPositionCalc();
			}
		}
		vector<Force*>::iterator f=force.begin();//�͂̌v�Z
		for(;f!=force.end();f++)(*f)->update();
		for(p=particle.begin();p!=particle.end();p++){//���Ԕ��W
			vector<Particle*>::iterator r=p->second.begin();
			for(;r!=p->second.end();r++){
				(*r)->updateEuler();
			}
		}
		if(itime%record==0){
			cout << "simulation steps: "  << itime << endl;
			report();
		}
	}
}

void ParticleSimulator::run_4thOrderRungeKutta(){
	cout << "simulation starts using 4th order Runge Kutta method" << endl;
	time=0.0;
	velocityField->initial();
	int total=inudf->i("simulation.time.simulation_steps");
	int record=inudf->i("simulation.time.report_steps");
	double dt=inudf->d("simulation.time.dt");
	map<int,vector<Particle*> >::iterator p;
	vector<Force*>::iterator f;
	vector<Particle*>::iterator r;
	velocityField->time=&time;
	electricField->time=&time;
	report();
	for(int itime=1;itime<=total;itime++){
		time=dt*itime;
		velocityField->update();//���x��X�V
		/////RungeKutta1�i��/////
		for(p=particle.begin();p!=particle.end();p++){
			for(r=p->second.begin();r!=p->second.end();r++){
				(*r)->forceAndTorqueReset();
				(*r)->beadsPositionCalc();
			}
		}
		for(f=force.begin();f!=force.end();f++)(*f)->update();//�͂̌v�Z
		/////RungeKutta2�i��/////
		for(p=particle.begin();p!=particle.end();p++)
			for(r=p->second.begin();r!=p->second.end();r++)
				(*r)->updateRungeKutta1();
		for(p=particle.begin();p!=particle.end();p++){
			for(r=p->second.begin();r!=p->second.end();r++){
				(*r)->forceAndTorqueReset();
				(*r)->beadsPositionCalc();
			}
		}
		for(f=force.begin();f!=force.end();f++)(*f)->update();//�͂̌v�Z
		for(p=particle.begin();p!=particle.end();p++)
			for(r=p->second.begin();r!=p->second.end();r++)
				(*r)->updateRungeKutta2();
		/////RungeKutta3�i��/////
		for(p=particle.begin();p!=particle.end();p++){
			for(r=p->second.begin();r!=p->second.end();r++){
				(*r)->forceAndTorqueReset();
				(*r)->beadsPositionCalc();
			}
		}
		for(f=force.begin();f!=force.end();f++)(*f)->update();//�͂̌v�Z
		for(p=particle.begin();p!=particle.end();p++)
			for(r=p->second.begin();r!=p->second.end();r++)
				(*r)->updateRungeKutta3();
		/////RungeKutta4�i��/////
		for(p=particle.begin();p!=particle.end();p++){
			for(r=p->second.begin();r!=p->second.end();r++){
				(*r)->forceAndTorqueReset();
				(*r)->beadsPositionCalc();
			}
		}
		for(f=force.begin();f!=force.end();f++)(*f)->update();//�͂̌v�Z
		for(p=particle.begin();p!=particle.end();p++)
			for(r=p->second.begin();r!=p->second.end();r++)
				(*r)->updateRungeKutta4();
		if(itime%record==0){
			cout << "simulation steps: "  << itime << endl;
			report();
		}
	}
}	

void ParticleSimulator::report(){
	outudf->newRecord();
	Tensor3x3 stress;
	Vector3d induced_dipole;
	map<int,vector<Particle*> >::iterator p=particle.begin();
	for(;p!=particle.end();p++){
		vector<Particle*>::iterator r=p->second.begin();
		for(;r!=p->second.end();r++){
			int index=(*r)->index;
			Location loc("simulation_result.particle[]",INDEX(index));
			Vector3d pos=(*r)->position;EulerAngle angle=(*r)->eulerAngle;
			Vector3d v=(*r)->velocity;Vector3d o=(*r)->angularVelocity;
			outudf->put(loc.sub("position.x"),pos.x);outudf->put(loc.sub("position.y"),pos.y);
			outudf->put(loc.sub("position.z"),pos.z);
			outudf->put(loc.sub("quaternion.xi"),angle.xi);outudf->put(loc.sub("quaternion.eta"),angle.eta);
			outudf->put(loc.sub("quaternion.zeta"),angle.zeta);outudf->put(loc.sub("quaternion.chi"),angle.chi);
			outudf->put(loc.sub("velocity.x"),v.x);outudf->put(loc.sub("velocity.y"),v.y);
			outudf->put(loc.sub("velocity.z"),v.z);
			outudf->put(loc.sub("angular_velocity.x"),o.x);outudf->put(loc.sub("angular_velocity.y"),o.y);
			outudf->put(loc.sub("angular_velocity.z"),o.z);
			stress+=(*r)->getStress();
		}
	}
	stress/=volume;//induced_dipole/=volume;
	Location loc("simulation_result.stress");
	outudf->put(loc.sub("x.x"),stress.x.x);outudf->put(loc.sub("x.y"),stress.x.y);outudf->put(loc.sub("x.z"),stress.x.z);
	outudf->put(loc.sub("y.x"),stress.y.x);outudf->put(loc.sub("y.y"),stress.y.y);outudf->put(loc.sub("y.z"),stress.y.z);
	outudf->put(loc.sub("z.x"),stress.z.x);outudf->put(loc.sub("z.y"),stress.z.y);outudf->put(loc.sub("z.z"),stress.z.z);
	//loc.seek("simulation_result.induced_current_dipole");
	//outudf->put(loc.sub("x"),induced_dipole.x);outudf->put(loc.sub("y"),induced_dipole.y);outudf->put(loc.sub("z"),induced_dipole.z);
	double t=timer.get();
	outudf->put("simulation_result.cpu_time",t);
	outudf->write();
}

