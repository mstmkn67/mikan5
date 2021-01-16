#include "UdfFlow.h"

////////////////////////////UdfFlow2d//////////////////////////////////
UdfFlow2d::UdfFlow2d(const string& udfFileName,double l,double t,double e){
	udf=new UDFManager(udfFileName);
	avsl=new AprxVector2dSquareLatticeField2d;
	length_ratio=udf->d("unit_parameter.length")/l;
	cout << "\t\tlength ratio: " << length_ratio << endl;
	//angular_velocity_ratio=udf->d("unit_parameter.time")/t;
	velocity_ratio=length_ratio*t/udf->d("unit_parameter.time");
	cout << "\t\tvelocity ratio: " << velocity_ratio << endl;
	//�z��̊m��
	int latticeNumber[2]={udf->i("system_conditions.lattice_number.x"),
		                    udf->i("system_conditions.lattice_number.y")};
	field.resize(latticeNumber[0],latticeNumber[1]);
	//���`��ԃN���X�̏���
	avsl->setField(&field);
	int minLattice[2]={udf->i("system_conditions.lattice1.x"),
		                 udf->i("system_conditions.lattice1.y")};
	int maxLattice[2]={udf->i("system_conditions.lattice2.x"),
		                 udf->i("system_conditions.lattice2.y")};
	Vector2d minLatticeCoord(udf->d("system_conditions.position1.x")*length_ratio,
		                       udf->d("system_conditions.position1.y")*length_ratio);
	Vector2d maxLatticeCoord(udf->d("system_conditions.position2.x")*length_ratio,
		                       udf->d("system_conditions.position2.y")*length_ratio);
	avsl->setMin(minLattice,minLatticeCoord);
	avsl->setMax(maxLattice,maxLatticeCoord);
	avsl->calcDr();
	//power=udf->d("system_conditions.scale_factor_of_velocity");
	min=minLatticeCoord;
	size=maxLatticeCoord-minLatticeCoord;
}

UdfFlow2d::~UdfFlow2d(){
	delete avsl;
	//delete udf; //�G���[���o���H
}

void UdfFlow2d::readData(){
	int* size=field.size();
	//�f�[�^�̓ǂݍ���
	for(int i=0;i<size[0];i++){
		for(int j=0;j<size[1];j++){
			field[i][j]=Vector2d(udf->d(Location("Field.velocity[].x",INDEX(j*size[0]+i)))*velocity_ratio,
				                   udf->d(Location("Field.velocity[].y",INDEX(j*size[0]+i)))*velocity_ratio);
		}
	}
}


Vector3d UdfFlow2d::getVelocity(const Vector3d& coord){
	Vector2d r=reset_coord(coord);
	Vector2d t=(*avsl)(Vector2d(r.x,r.y));
	return Vector3d(t.x,t.y,0);
}

Vector3d UdfFlow2d::getAngularVelocity(const Vector3d& coord){
	Vector2d r=reset_coord(coord);
	return Vector3d(0,0,0.5*avsl->rot(Vector2d(r.x,r.y)));
}

Tensor3x3 UdfFlow2d::getRateOfStrainTensor(const Vector3d& coord){
	Vector2d r=reset_coord(coord);
	double dvxdx=avsl->differential(r,0,0);
	double dvxdy=avsl->differential(r,0,1);
	double dvydx=avsl->differential(r,1,0);
	double dvydy=avsl->differential(r,1,1);
	return Tensor3x3(dvxdx,0.5*(dvydx+dvxdy),0.0,
		            0.5*(dvydx+dvxdy),dvydy,0.0,
					 0.0,0.0,0.0);
}

Vector2d UdfFlow2d::reset_coord(const Vector2d& r){
	return Vector2d(r.x-floor((r.x-min.x)/size.x)*size.x,
	                r.y-floor((r.y-min.y)/size.y)*size.y);
}

////////////////////////////UdfFlow3d//////////////////////////////////
UdfFlow3d::UdfFlow3d(const string& udfFileName,double l,double t,double e){
	udf=new UDFManager(udfFileName);
	avsl=new AprxVector3dSquareLatticeField3d;
	length_ratio=udf->d("unit_parameter.length")/l;
	cout << "\t\tlength ratio: " << length_ratio << endl;
	//angular_velocity_ratio=udf->d("unit_parameter.time")/t;
	velocity_ratio=length_ratio*t/udf->d("unit_parameter.time");
	cout << "\t\tvelocity ratio: " << velocity_ratio << endl;
	//�z��̊m��
	int latticeNumber[3]={udf->i("system_conditions.lattice_number.x"),
	 	                    udf->i("system_conditions.lattice_number.y"),
		                    udf->i("system_conditions.lattice_number.z")};
	field.resize(latticeNumber[0],latticeNumber[1],latticeNumber[2]);
	//���`��ԃN���X�̏���
	avsl->setField(&field);
	int minLattice[3]={udf->i("system_conditions.lattice1.x"),
		                 udf->i("system_conditions.lattice1.y"),
	                   udf->i("system_conditions.lattice1.z")};
	int maxLattice[3]={udf->i("system_conditions.lattice2.x"),
		                 udf->i("system_conditions.lattice2.y"),
                     udf->i("system_conditions.lattice2.z")};
		Vector3d minLatticeCoord(udf->d("system_conditions.position1.x")*length_ratio,
		                         udf->d("system_conditions.position1.y")*length_ratio,
		                         udf->d("system_conditions.position1.z")*length_ratio);
	Vector3d maxLatticeCoord(udf->d("system_conditions.position2.x")*length_ratio,
	                         udf->d("system_conditions.position2.y")*length_ratio,
	                         udf->d("system_conditions.position2.z")*length_ratio);
	avsl->setMin(minLattice,minLatticeCoord);
	avsl->setMax(maxLattice,maxLatticeCoord);
	avsl->calcDr();
	//power=udf->d("system_conditions.scale_factor_of_velocity");
	min=minLatticeCoord;
	size=maxLatticeCoord-minLatticeCoord;
}

UdfFlow3d::~UdfFlow3d(){
	delete avsl;
	//delete udf; //�G���[���o���H
}

void UdfFlow3d::readData(){
	//�f�[�^�̓ǂݍ���
	int* size=field.size();
	for(int i=0;i<size[0];i++){
		for(int j=0;j<size[1];j++){
			for(int k=0;k<size[2];k++){
				int index=k*size[0]*size[1]+j*size[0]+i;
				field[i][j][k]=velocity_ratio*Vector3d(udf->d(Location("Field.velocity[].x",INDEX(index)))*velocity_ratio,
				                              udf->d(Location("Field.velocity[].y",INDEX(index)))*velocity_ratio,
				                              udf->d(Location("Field.velocity[].z",INDEX(index)))*velocity_ratio);
			}
		}
	}
}

Vector3d UdfFlow3d::getVelocity(const Vector3d& coord){
	Vector3d r=reset_coord(coord);
	return (*avsl)(r);
}

Vector3d UdfFlow3d::getAngularVelocity(const Vector3d& coord){
	Vector3d r=reset_coord(coord);
	return 0.5*avsl->rot(r);
}

Tensor3x3 UdfFlow3d::getRateOfStrainTensor(const Vector3d& coord){
	Vector3d r=reset_coord(coord);
	double dvxdx=avsl->differential(r,0,0);
	double dvxdy=avsl->differential(r,0,1);
	double dvxdz=avsl->differential(r,0,2);//
	double dvydx=avsl->differential(r,1,0);
	double dvydy=avsl->differential(r,1,1);
	double dvydz=avsl->differential(r,1,2);//
	double dvzdx=avsl->differential(r,2,0);
	double dvzdy=avsl->differential(r,2,1);
	double dvzdz=avsl->differential(r,2,2);//
	return Tensor3x3(dvxdx,0.5*(dvydx+dvxdy),0.5*(dvzdx+dvxdz),
	                 0.5*(dvydx+dvxdy),dvydy,0.5*(dvzdy+dvydz),
	                 0.5*(dvzdx+dvxdz),0.5*(dvzdy+dvydz),dvzdz);
}

Vector3d UdfFlow3d::reset_coord(const Vector3d& r){
	return Vector3d(r.x-floor((r.x-min.x)/size.x)*size.x,
	                r.y-floor((r.y-min.y)/size.y)*size.y,
	                r.z-floor((r.z-min.z)/size.z)*size.z);
}


///////////////////////////DynamicUdfFlow2d//////////////////////////////////
DynamicUdfFlow2d::DynamicUdfFlow2d(const string& udfFileName,double l,double t,double e,double d)
:dt(d){
	avsl1=new AprxVector2dSquareLatticeField2d;
	avsl2=new AprxVector2dSquareLatticeField2d;
	udf=new UDFManager(udfFileName);
	length_ratio=udf->d("unit_parameter.length")/l;
	cout << "\t\tlength ratio: " << length_ratio << endl;
	//angular_velocity_ratio=udf->d("unit_parameter.time")/t;
	velocity_ratio=length_ratio*t/udf->d("unit_parameter.time");
	cout << "\t\tvelocity ratio: " << velocity_ratio << endl;
	//�z��̊m��
	int latticeNumber[2]={udf->i("system_conditions.lattice_number.x"),
		                    udf->i("system_conditions.lattice_number.y")};
	field1.resize(latticeNumber[0],latticeNumber[1]);field2.resize(latticeNumber[0],latticeNumber[1]);
	//���`��ԃN���X�̏���
	avsl1->setField(&field1);avsl2->setField(&field2);
	int minLattice[2]={udf->i("system_conditions.lattice1.x"),
		                 udf->i("system_conditions.lattice1.y")};
	int maxLattice[2]={udf->i("system_conditions.lattice2.x"),
		                 udf->i("system_conditions.lattice2.y")};
	Vector2d minLatticeCoord(udf->d("system_conditions.position1.x")*length_ratio,
		                       udf->d("system_conditions.position1.y")*length_ratio);
	Vector2d maxLatticeCoord(udf->d("system_conditions.position2.x")*length_ratio,
		                       udf->d("system_conditions.position2.y")*length_ratio);
	avsl1->setMin(minLattice,minLatticeCoord);avsl1->setMax(maxLattice,maxLatticeCoord);avsl1->calcDr();
	avsl2->setMin(minLattice,minLatticeCoord);avsl2->setMax(maxLattice,maxLatticeCoord);avsl2->calcDr();
	//���R�[�h������̎��Ԃ�ǂݎ��
	timeStepsPerRecord=int(udf->d("system_conditions.time_per_record")/dt);
	cout << "\t\ttime steps in particle simulator per one record velocity field: "
	     << timeStepsPerRecord << endl;//timeStepsPerRecord=1001;
	//power=udf->d("system_conditions.scale_factor_of_velocity");
	min=minLatticeCoord;
	size=maxLatticeCoord-minLatticeCoord;
}

DynamicUdfFlow2d::~DynamicUdfFlow2d(){
	delete udf;
}

void DynamicUdfFlow2d::initial(){
	if(!udf->jump(0)){
		cout << "Can't jump to initial record" << endl;
		exit(1);
	}
	cout << "\t\tvelocity record: " << 0 << endl;
	int* size=field1.size();
	for(int i=0;i<size[0];i++){
		for(int j=0;j<size[1];j++){
			int index=j*size[0]+i;
			field1[i][j]=velocity_ratio*Vector2d(udf->d(Location("Field.velocity[].x",INDEX(index))),
			                            udf->d(Location("Field.velocity[].y",INDEX(index))));
		}
	}
	udf->nextRecord();
	cout << "\t\tvelocity record: " << 1 << endl;
	{for(int i=0;i<size[0];i++){
		for(int j=0;j<size[1];j++){
			int index=j*size[0]+i;
			field2[i][j]=velocity_ratio*Vector2d(udf->d(Location("Field.velocity[].x",INDEX(index))),
			                            udf->d(Location("Field.velocity[].y",INDEX(index))));
		}
	}}
	region=0;
}

void DynamicUdfFlow2d::readData(int record,Array2d<Vector2d>* field){
	udf->jump(record);
	//�f�[�^�̓ǂݍ���
	int* size=(*field).size();
	for(int i=0;i<size[0];i++){
		for(int j=0;j<size[1];j++){
			int index=j*size[0]+i;
			(*field)[i][j]=velocity_ratio*Vector2d(udf->d(Location("Field.velocity[].x",INDEX(index))),
			                              udf->d(Location("Field.velocity[].y",INDEX(index))));
		}
	}
}

void DynamicUdfFlow2d::update(){
	int newRegion=*time/timeStepsPerRecord;
	if(region==newRegion)return;//���R�[�h�̓ǂݑւ��̕K�v�͂Ȃ�
	if(newRegion-region!=1){
		cout << "\t\tvelocity record: " << newRegion << " and " << newRegion+1 << endl;
		if(newRegion%2==0){
			readData(newRegion,&field1);//�f�[�^�̓ǂݍ���
			readData(newRegion+1,&field2);//�f�[�^�̓ǂݍ���
		}else{
			readData(newRegion+1,&field1);//�f�[�^�̓ǂݍ���
			readData(newRegion,&field2);//�f�[�^�̓ǂݍ���
		}
	}else{
		cout << "\t\tvelocity record: " << newRegion+1 << endl;
		if(newRegion%2==0){
			readData(newRegion+1,&field2);//�f�[�^�̓ǂݍ���
		}else{
			readData(newRegion+1,&field1);//�f�[�^�̓ǂݍ���
		}
	}
	region=newRegion;
}

Vector3d DynamicUdfFlow2d::getVelocity(const Vector3d& coord){
	Vector2d r=reset_coord(coord);
	int itime=(*time)/dt;
	double ta=itime%timeStepsPerRecord;double tb=timeStepsPerRecord-ta;
	double a=ta/timeStepsPerRecord;double b=tb/timeStepsPerRecord;
	if(region%2==0){double n=a;a=b;b=n;}
	Vector2d t1=(*avsl1)(Vector2d(r.x,r.y));Vector2d t2=(*avsl2)(Vector2d(r.x,r.y)); 
	Vector2d t=a*t1+b*t2;
	return Vector3d(t.x,t.y,0);
}

Vector3d DynamicUdfFlow2d::getAngularVelocity(const Vector3d& coord){
	Vector2d r=reset_coord(coord);
	int itime=(*time)/dt;
	double ta=itime%timeStepsPerRecord;double tb=timeStepsPerRecord-ta;
	double a=ta/timeStepsPerRecord;double b=tb/timeStepsPerRecord;
	if(region%2==0){double n=a;a=b;b=n;}
	return Vector3d(0,0,0.5*(a*avsl1->rot(Vector2d(r.x,r.y))+b*avsl2->rot(Vector2d(r.x,r.y))));
}

Tensor3x3 DynamicUdfFlow2d::getRateOfStrainTensor(const Vector3d& coord){
	Vector2d r=reset_coord(coord);
	//cout <<  timeStepsPerRecord << endl;
	int itime=(*time)/dt;
	double ta=itime%timeStepsPerRecord;double tb=timeStepsPerRecord-ta;
	double a=ta/timeStepsPerRecord;double b=tb/timeStepsPerRecord;
	if(region%2==0){double n=a;a=b;b=n;}
	double dvxdx=a*avsl1->differential(r,0,0)+b*avsl2->differential(r,0,0);
	double dvxdy=a*avsl1->differential(r,0,1)+b*avsl2->differential(r,0,1);
	double dvydx=a*avsl1->differential(r,1,0)+b*avsl2->differential(r,1,0);
	double dvydy=a*avsl1->differential(r,1,1)+b*avsl1->differential(r,1,1);
	return Tensor3x3(dvxdx,0.5*(dvydx+dvxdy),0.0,
		            0.5*(dvydx+dvxdy),dvydy,0.0,
					 0.0,0.0,0.0);
}

Vector2d DynamicUdfFlow2d::reset_coord(const Vector2d& r){
	return Vector2d(r.x-floor((r.x-min.x)/size.x)*size.x,
	                r.y-floor((r.y-min.y)/size.y)*size.y);
}

//////////////////////////DynamicUdfFlow3d//////////////////////////////////
DynamicUdfFlow3d::DynamicUdfFlow3d(const string& udfFileName,double l,double t,double e,double d)
:dt(d){
	avsl1=new AprxVector3dSquareLatticeField3d;
	avsl2=new AprxVector3dSquareLatticeField3d;
	udf=new UDFManager(udfFileName);
	length_ratio=udf->d("unit_parameter.length")/l;
	cout << "\t\tlength ratio: " << length_ratio << endl;
	//angular_velocity_ratio=udf->d("unit_parameter.time")/t;
	velocity_ratio=length_ratio*t/udf->d("unit_parameter.time");
	cout << "\t\tvelocity ratio: " << velocity_ratio << endl;
	//�z��̊m��
	int latticeNumber[3]={udf->i("system_conditions.lattice_number.x"),
	 	                    udf->i("system_conditions.lattice_number.y"),
		                    udf->i("system_conditions.lattice_number.z")};
	field1.resize(latticeNumber[0],latticeNumber[1],latticeNumber[2]);
	field2.resize(latticeNumber[0],latticeNumber[1],latticeNumber[2]);
	//���`��ԃN���X�̏���
	avsl1->setField(&field1);avsl2->setField(&field2);
	int minLattice[3]={udf->i("system_conditions.lattice1.x"),
		                 udf->i("system_conditions.lattice1.y"),
	                   udf->i("system_conditions.lattice1.z")};
	int maxLattice[3]={udf->i("system_conditions.lattice2.x"),
		                 udf->i("system_conditions.lattice2.y"),
                     udf->i("system_conditions.lattice2.z")};
   	Vector3d minLatticeCoord(length_ratio*udf->d("system_conditions.position1.x"),
		                         length_ratio*udf->d("system_conditions.position1.y"),
		                         length_ratio*udf->d("system_conditions.position1.z"));
	Vector3d maxLatticeCoord(length_ratio*udf->d("system_conditions.position2.x"),
	                         length_ratio*udf->d("system_conditions.position2.y"),
	                         length_ratio*udf->d("system_conditions.position2.z"));
	avsl1->setMin(minLattice,minLatticeCoord);avsl1->setMax(maxLattice,maxLatticeCoord);avsl1->calcDr();
	avsl2->setMin(minLattice,minLatticeCoord);avsl2->setMax(maxLattice,maxLatticeCoord);avsl2->calcDr();
	//���R�[�h������̎��Ԃ�ǂݎ��
	timeStepsPerRecord=int(udf->d("system_conditions.time_per_record")/dt);
	cout << "\t\ttime steps in particle simulator per one record velocity field: "
	     << timeStepsPerRecord << endl;
	//power=udf->d("system_conditions.scale_factor_of_velocity");
	min=minLatticeCoord;
	size=maxLatticeCoord-minLatticeCoord;
}

DynamicUdfFlow3d::~DynamicUdfFlow3d(){
	delete udf;
}

void DynamicUdfFlow3d::initial(){
	if(!udf->jump(0)){
		cout << "Can't jump to initial record" << endl;
		exit(1);
	}
	int* size=field1.size();
	for(int i=0;i<size[0];i++){
		for(int j=0;j<size[1];j++){
			for(int k=0;k<size[2];k++){
				int index=k*size[0]*size[1]+j*size[0]+i;
				field1[i][j][k]=velocity_ratio*Vector3d(udf->d(Location("Field.velocity[].x",INDEX(index))),
					                     udf->d(Location("Field.velocity[].y",INDEX(index))),
									     udf->d(Location("Field.velocity[].z",INDEX(index))));
			}
		}
	}
	udf->nextRecord();
	{for(int i=0;i<size[0];i++){
		for(int j=0;j<size[1];j++){
			for(int k=0;k<size[2];k++){
				int index=k*size[0]*size[1]+j*size[0]+i;
				field2[i][j][k]=velocity_ratio*Vector3d(udf->d(Location("Field.velocity[].x",INDEX(index))),
					                     udf->d(Location("Field.velocity[].y",INDEX(index))),
									     udf->d(Location("Field.velocity[].z",INDEX(index))));
			}
		}
	}}
	region=0;
}

void DynamicUdfFlow3d::readData(int record,Array3d<Vector3d>* field){
	udf->jump(record);
	//�f�[�^�̓ǂݍ���
	int* size=(*field).size();
	for(int i=0;i<size[0];i++){
		for(int j=0;j<size[1];j++){
			for(int k=0;k<size[2];k++){
				int index=k*size[0]*size[1]+j*size[0]+i;
				(*field)[i][j][k]=velocity_ratio*Vector3d(udf->d(Location("Field.velocity[].x",INDEX(index))),
					                       udf->d(Location("Field.velocity[].y",INDEX(index))),
									       udf->d(Location("Field.velocity[].z",INDEX(index))));
			}
		}
	}
}

void DynamicUdfFlow3d::update(){
	int newRegion=(*time)/timeStepsPerRecord;
	if(region==newRegion)return;//���R�[�h�̓ǂݑւ��̕K�v�͂Ȃ�
	if(newRegion-region!=1){
		cout << "\t\tvelocity record: " << newRegion << " and " << newRegion+1 << endl;
		if(newRegion%2==0){
			readData(newRegion,&field1);//�f�[�^�̓ǂݍ���
			readData(newRegion+1,&field2);//�f�[�^�̓ǂݍ���
		}else{
			readData(newRegion+1,&field1);//�f�[�^�̓ǂݍ���
			readData(newRegion,&field2);//�f�[�^�̓ǂݍ���
		}
	}else{
		cout << "\t\tvelocity record: " << newRegion+1 << endl;
		if(newRegion%2==0){
			readData(newRegion+1,&field2);//�f�[�^�̓ǂݍ���
		}else{
			readData(newRegion+1,&field1);//�f�[�^�̓ǂݍ���
		}
	}
	region=newRegion;
}


Vector3d DynamicUdfFlow3d::getVelocity(const Vector3d& coord){
	Vector3d r=reset_coord(coord);
	int itime=(*time)/dt;
	double ta=itime%timeStepsPerRecord;double tb=timeStepsPerRecord-ta;
	double a=ta/timeStepsPerRecord;double b=tb/timeStepsPerRecord;
	if(region%2==0){double n=a;a=b;b=n;}
	return a*(*avsl1)(r)+b*(*avsl2)(r);
}

Vector3d DynamicUdfFlow3d::getAngularVelocity(const Vector3d& coord){
	Vector3d r=reset_coord(coord);
	int itime=(*time)/dt;
	double ta=itime%timeStepsPerRecord;double tb=timeStepsPerRecord-ta;
	double a=ta/timeStepsPerRecord;double b=tb/timeStepsPerRecord;
	if(region%2==0){double n=a;a=b;b=n;}
	return 0.5*(a*avsl1->rot(r)+b*avsl2->rot(r));
}

Tensor3x3 DynamicUdfFlow3d::getRateOfStrainTensor(const Vector3d& coord){
	Vector3d r=reset_coord(coord);
	int itime=(*time)/dt;
	double ta=itime%timeStepsPerRecord;double tb=timeStepsPerRecord-ta;
	double a=ta/timeStepsPerRecord;double b=tb/timeStepsPerRecord;
	if(region%2==0){double n=a;a=b;b=n;}
	double dvxdx=0.5*(a*avsl1->differential(r,0,0)+b*avsl2->differential(r,0,0));
	double dvxdy=0.5*(a*avsl1->differential(r,0,1)+b*avsl2->differential(r,0,1));
	double dvxdz=0.5*(a*avsl1->differential(r,0,2)+b*avsl2->differential(r,0,2));//
	double dvydx=0.5*(a*avsl1->differential(r,1,0)+b*avsl2->differential(r,1,0));
	double dvydy=0.5*(a*avsl1->differential(r,1,1)+b*avsl2->differential(r,1,1));
	double dvydz=0.5*(a*avsl1->differential(r,1,0)+b*avsl2->differential(r,1,2));//
	double dvzdx=0.5*(a*avsl1->differential(r,2,0)+b*avsl2->differential(r,2,0));
	double dvzdy=0.5*(a*avsl1->differential(r,2,1)+b*avsl2->differential(r,2,1));
	double dvzdz=0.5*(a*avsl1->differential(r,2,2)+b*avsl2->differential(r,2,2));
	return Tensor3x3(dvxdx,0.5*(dvydx+dvxdy),0.5*(dvzdx+dvxdz),
	                 0.5*(dvydx+dvxdy),dvydy,0.5*(dvzdy+dvydz),
	                 0.5*(dvzdx+dvxdz),0.5*(dvzdy+dvydz),dvzdz);
}

Vector3d DynamicUdfFlow3d::reset_coord(const Vector3d& r){
	return Vector3d(r.x-floor((r.x-min.x)/size.x)*size.x,
	                r.y-floor((r.y-min.y)/size.y)*size.y,
	                r.z-floor((r.z-min.z)/size.z)*size.z);
}

//�ȉ��͎��ԕ����ɐ��`��Ԃ���Ă��Ȃ��o�[�W����
////////////////////////////DynamicUdfFlow2d//////////////////////////////////
//DynamicUdfFlow2d::DynamicUdfFlow2d(const string& udfFileName):UdfFlow2d(udfFileName){
//	//���R�[�h������̎��Ԃ�ǂݎ��
//	timePerRecord=udf->i("SystemConditions.TimeStepsPerRecord");
//}
//
//void DynamicUdfFlow2d::initial(){
//	if(!udf->jump(0)){
//		cout << "Can't jump to initial record" << endl;
//		exit(1);
//	}
//	readData();//�f�[�^�̓ǂݍ���
//}
//
//void DynamicUdfFlow2d::update(){
//	if((*time)%timePerRecord!=0)return;
//	if(!udf->nextRecord()){
//		cout << "Can't read the next record" << endl;
//		exit(1);
//	}
//	readData();//�f�[�^�̓ǂݍ���
//}

////////////////////////////DynamicUdfFlow2d//////////////////////////////////
//DynamicUdfFlow3d::DynamicUdfFlow3d(const string& udfFileName):UdfFlow3d(udfFileName){
//	//���R�[�h������̎��Ԃ�ǂݎ��
//	timePerRecord=udf->i("SystemConditions.TimeStepsPerRecord");
//}
//
//void DynamicUdfFlow3d::initial(){
//	if(!udf->jump(0)){
//		cout << "Can't jump to initial record" << endl;
//		exit(1);
//	}
//	readData();//�f�[�^�̓ǂݍ���
//}
//
//void DynamicUdfFlow3d::update(){
//	if((*time)%timePerRecord!=0)return;
//	if(!udf->nextRecord()){
//		cout << "Can't read the next record" << endl;
//		exit(1);
//	}
//	readData();//�f�[�^�̓ǂݍ���
//}

