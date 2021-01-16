#include "ResistanceMobilityTensor.h"
/////ResistanceTensor/////
ResistanceTensor::ResistanceTensor(const ResistanceTensor& r){
	A=r.A;B=r.B;C=r.C;G=r.G;H=r.H;K=r.K;
}
ResistanceTensor::ResistanceTensor(const MobilityTensor& m){
	ResistanceTensor r=m.getResistanceTensor();
	A=r.A;B=r.B;C=r.C;G=r.G;H=r.H;K=r.K;
}
ResistanceTensor::ResistanceTensor(
	const Tensor3x3& A_,const Tensor3x3& B_,const Tensor3x3& C_,
	const Tensor3x3x3& G_,const Tensor3x3x3& H_,const Tensor3x3x3x3& K_){
	A=A_;B=B_;C=C_;G=G_;H=H_;K=K_;
}
void ResistanceTensor::set(const Tensor3x3& A_,const Tensor3x3& B_,const Tensor3x3& C_,
	                         const Tensor3x3x3& G_,const Tensor3x3x3& H_,const Tensor3x3x3x3& K_){
	A=A_;B=B_;C=C_;G=G_;H=H_;K=K_;
}
MobilityTensor ResistanceTensor::getMobilityTensor()const{
	Tensor3x3 a,b,c;Tensor3x3x3 g,h;Tensor3x3x3x3 k;
	tensor_utility::convertResistance2Mobility(A,B,C,G,H,K,a,b,c,g,h,k);
	return MobilityTensor(a,b,c,g,h,k);
}

ResistanceTensor ResistanceTensor::getTranslateTensor(const Vector3d& r)const{
	Tensor3x3 C_=tensor_utility::translation_C(A,B,C,-r);//'-' means moving the center 
	Tensor3x3x3x3 K_=tensor_utility::translation_K(A,G,K,-r);
	Tensor3x3x3 H_=tensor_utility::translation_H(A,B,G,H,-r);
	Tensor3x3 B_=tensor_utility::translation_B(A,B,-r);
	Tensor3x3x3 G_=tensor_utility::translation_G(A,G,-r);
	return ResistanceTensor(A,B_,C_,G_,H_,K_);
}

void ResistanceTensor::translate(const Vector3d& r){
	C=tensor_utility::translation_C(A,B,C,-r);//'-' means moving the center 
	K=tensor_utility::translation_K(A,G,K,-r);
	H=tensor_utility::translation_H(A,B,G,H,-r);
	B=tensor_utility::translation_B(A,B,-r);
	G=tensor_utility::translation_G(A,G,-r);
}

void ResistanceTensor::rotate(const Tensor3x3& QQ){
	A=convert(QQ,A);B=convert(QQ,B);C=convert(QQ,C);
	G=convert(QQ,G);H=convert(QQ,H);K=convert(QQ,K);
}

void ResistanceTensor::rotate(const Tensor3x3& QQ,const Vector3d& center){
	translate(-center);
	A=convert(QQ,A);B=convert(QQ,B);C=convert(QQ,C);
	G=convert(QQ,G);H=convert(QQ,H);K=convert(QQ,K);
	translate(QQ*center);
}

Vector3d ResistanceTensor::getCenterOfResistance()const{
	return tensor_utility::calcCenterOfResistance(A,B);
}

ResistanceTensor& ResistanceTensor::operator=(const ResistanceTensor& r){
	A=r.A;B=r.B;C=r.C;G=r.G;H=r.H;K=r.K;
	return *this;
}
ResistanceTensor& ResistanceTensor::operator*=(double s){
	A*=s;B*=s;C*=s;G*=s;H*=s;K*=s;
	return *this;
}
ResistanceTensor& ResistanceTensor::operator/=(double s){
	A/=s;B/=s;C/=s;G/=s;H/=s;K/=s;
	return *this;
}
ResistanceTensor ResistanceTensor::operator - () const{
	return ResistanceTensor(-A,-B,-C,-G,-H,-K);
}

Vector3d ResistanceTensor::selectCenterOfResistance(){
	Vector3d r=-getCenterOfResistance();
	translate(r);
	return r;
}
Vector3d ResistanceTensor::selectCenterOfMobility(){
	Vector3d r=-getMobilityTensor().getCenterOfMobility();
	translate(r);
	return r;
}
Tensor3x3 ResistanceTensor::select_A_diagonalizedSystem(){
	return diagonalize_tensor(A);
}
Tensor3x3 ResistanceTensor::select_B_diagonalizedSystem(){
	return diagonalize_tensor(B);
}
Tensor3x3 ResistanceTensor::select_C_diagonalizedSystem(){
	return diagonalize_tensor(C);
}
Tensor3x3 ResistanceTensor::select_a_diagonalizedSystem(){
	Tensor3x3 a,b,c;
	tensor_utility::convertABC2abc(A,B,C,a,b,c);
	Vector3d r=tensor_utility::calcCenterOfMobility(b,c);
	Tensor3x3 a_=tensor_utility::translation_a(a,b,c,r);
	return diagonalize_tensor(a_,r);
}
Tensor3x3 ResistanceTensor::select_b_diagonalizedSystem(){
	Tensor3x3 a,b,c;
	tensor_utility::convertABC2abc(A,B,C,a,b,c);
	Vector3d r=tensor_utility::calcCenterOfMobility(b,c);
	Tensor3x3 b_=tensor_utility::translation_b(b,c,r);
	return diagonalize_tensor(b_,r);
}
Tensor3x3 ResistanceTensor::select_c_diagonalizedSystem(){
	Tensor3x3 a,b,c;
	tensor_utility::convertABC2abc(A,B,C,a,b,c);
	Vector3d r=tensor_utility::calcCenterOfMobility(b,c);
	return diagonalize_tensor(c,r);
}

Tensor3x3 ResistanceTensor::diagonalize_tensor(const Tensor3x3& AA){
	double v[3];double R[9];
	for(int i=0;i<3;i++){
		for(int j=i;j<3;j++){
			R[3*i+j]=AA.get(i,j);
		}
	}
	lap::get_eigenvalues_and_vectors_symmetric(3,R,v);
	//eigen vectors
	Vector3d e0(R[3*0+0],R[3*0+1],R[3*0+2]),e1(R[3*1+0],R[3*1+1],R[3*1+2]);
	Vector3d e2=(e0^e1);//keep the right-hand-system
	Tensor3x3 Q=Tensor3x3(e0,e1,e2);rotate(Q);
	return Q;
}

Tensor3x3 ResistanceTensor::diagonalize_tensor(const Tensor3x3& AA,const Vector3d& center){
	double v[3];double R[9];
	for(int i=0;i<3;i++){
		for(int j=i;j<3;j++){
			R[3*i+j]=AA.get(i,j);
		}
	}
	lap::get_eigenvalues_and_vectors_symmetric(3,R,v);
	//eigen vectors
	Vector3d e0(R[3*0+0],R[3*0+1],R[3*0+2]),e1(R[3*1+0],R[3*1+1],R[3*1+2]);
	Vector3d e2=(e0^e1);//keep the right-hand-system
	Tensor3x3 Q=Tensor3x3(e0,e1,e2);rotate(Q,center);
	return Q;
}

//friend function
ResistanceTensor operator+(const ResistanceTensor& r1,const ResistanceTensor& r2){
	return ResistanceTensor(
		r1.A+r2.A,r1.B+r2.B,r1.C+r2.C,r1.G+r2.G,r1.H+r2.H,r1.K+r2.K);
}
ResistanceTensor operator-(const ResistanceTensor& r1,const ResistanceTensor& r2){
	return ResistanceTensor(
		r1.A-r2.A,r1.B-r2.B,r1.C-r2.C,r1.G-r2.G,r1.H-r2.H,r1.K-r2.K);
}
ResistanceTensor operator*(double s, const ResistanceTensor &r){
	return ResistanceTensor(
		s*r.A,s*r.B,s*r.C,s*r.G,s*r.H,s*r.K);
}
ResistanceTensor operator/(const ResistanceTensor& r,double s){
	return ResistanceTensor(
		r.A/s,r.B/s,r.C/s,r.G/s,r.H/s,r.K/s);
}
/////MobilityTensor/////
MobilityTensor::MobilityTensor(const MobilityTensor& m){
	a=m.a;b=m.b;c=m.c;g=m.g;h=m.h;k=m.k;
}
MobilityTensor::MobilityTensor(const ResistanceTensor& r){
	MobilityTensor m=r.getMobilityTensor();
	a=m.a;b=m.b;c=m.c;g=m.g;h=m.h;k=m.k;
}
MobilityTensor::MobilityTensor(
	const Tensor3x3& a_,const Tensor3x3& b_,const Tensor3x3& c_,
	const Tensor3x3x3& g_,const Tensor3x3x3& h_,const Tensor3x3x3x3& k_){
	a=a_;b=b_;c=c_;g=g_;h=h_;k=k_;
}
void MobilityTensor::set(const Tensor3x3& a_,const Tensor3x3& b_,const Tensor3x3& c_,
	                       const Tensor3x3x3& g_,const Tensor3x3x3& h_,const Tensor3x3x3x3& k_){
	a=a_;b=b_;c=c_;g=g_;h=h_;k=k_;
}
ResistanceTensor MobilityTensor::getResistanceTensor()const{
	Tensor3x3 A,B,C;Tensor3x3x3 G,H;Tensor3x3x3x3 K;
	tensor_utility::convertResistance2Mobility(a,b,c,g,h,k,A,B,C,G,H,K);
	return ResistanceTensor(A,B,C,G,H,K);
}

MobilityTensor MobilityTensor::getTranslateTensor(const Vector3d& r)const{
	Tensor3x3 a_=tensor_utility::translation_a(a,b,c,-r);//'-' means moving the center 
	Tensor3x3 b_=tensor_utility::translation_b(b,c,-r);
	Tensor3x3x3 g_=tensor_utility::translation_g(g,h,-r);
	return MobilityTensor(a_,b_,c,g_,h,k);
}

void MobilityTensor::translate(const Vector3d& r){
	a=tensor_utility::translation_a(a,b,c,-r);//'-' means moving the center 
	b=tensor_utility::translation_b(b,c,-r);
	g=tensor_utility::translation_g(g,h,-r);
}

void MobilityTensor::rotate(const Tensor3x3& QQ){
	a=convert(QQ,a);b=convert(QQ,b);c=convert(QQ,c);
	g=convert(QQ,g);h=convert(QQ,h);k=convert(QQ,k);
}
void MobilityTensor::rotate(const Tensor3x3& QQ,const Vector3d& center){
	translate(-center);
	a=convert(QQ,a);b=convert(QQ,b);c=convert(QQ,c);
	g=convert(QQ,g);h=convert(QQ,h);k=convert(QQ,k);
	translate(QQ*center);
}

Vector3d MobilityTensor::getCenterOfMobility()const{
	return tensor_utility::calcCenterOfMobility(b,c);
}

MobilityTensor& MobilityTensor::operator=(const MobilityTensor& m){
	a=m.a;b=m.b;c=m.c;g=m.g;h=m.h;k=m.k;
	return *this;
}
MobilityTensor& MobilityTensor::operator *= (double s){
	a*=s;b*=s;c*=s;g*=s;h*=s;k*=s;
	return *this;
}
MobilityTensor& MobilityTensor::operator /= (double s){
	a/=s;b/=s;c/=s;g/=s;h/=s;k/=s;
	return *this;
}
MobilityTensor MobilityTensor::operator-()const{
	return MobilityTensor(-a,-b,-c,-g,-h,-k);
}

Vector3d MobilityTensor::selectCenterOfResistance(){
	Vector3d r=-getResistanceTensor().getCenterOfResistance();
	translate(r);
	return r;
}
Vector3d MobilityTensor::selectCenterOfMobility(){
	Vector3d r=-getCenterOfMobility();
	translate(r);
	return r;
}
Tensor3x3 MobilityTensor::select_A_diagonalizedSystem(){
	Tensor3x3 A,B,C;
	tensor_utility::convertABC2abc(a,b,c,A,B,C);
	Vector3d r=tensor_utility::calcCenterOfResistance(A,B);
	return diagonalize_tensor(A,r);
}
Tensor3x3 MobilityTensor::select_B_diagonalizedSystem(){
	Tensor3x3 A,B,C;
	tensor_utility::convertABC2abc(a,b,c,A,B,C);
	Vector3d r=tensor_utility::calcCenterOfResistance(A,B);
	Tensor3x3 B_=tensor_utility::translation_B(A,B,r);
	return diagonalize_tensor(B_,r);
}
Tensor3x3 MobilityTensor::select_C_diagonalizedSystem(){
	Tensor3x3 A,B,C;
	tensor_utility::convertABC2abc(a,b,c,A,B,C);
	Vector3d r=tensor_utility::calcCenterOfResistance(A,B);
	Tensor3x3 C_=tensor_utility::translation_C(A,B,C,r);
	return diagonalize_tensor(C_,r);
}
Tensor3x3 MobilityTensor::select_a_diagonalizedSystem(){
	return diagonalize_tensor(a);
}
Tensor3x3 MobilityTensor::select_b_diagonalizedSystem(){
	return diagonalize_tensor(b);
}
Tensor3x3 MobilityTensor::select_c_diagonalizedSystem(){
	return diagonalize_tensor(c);
}

Tensor3x3 MobilityTensor::diagonalize_tensor(const Tensor3x3& AA){
	double v[3];double R[9];
	for(int i=0;i<3;i++){
		for(int j=i;j<3;j++){
			R[3*i+j]=AA.get(i,j);
		}
	}
	lap::get_eigenvalues_and_vectors_symmetric(3,R,v);
	//eigen vectors
	Vector3d e0(R[3*0+0],R[3*0+1],R[3*0+2]),e1(R[3*1+0],R[3*1+1],R[3*1+2]);
	Vector3d e2=(e0^e1);//keep the right-hand-system
	Tensor3x3 Q=Tensor3x3(e0,e1,e2);rotate(Q);
	return Q;
}

Tensor3x3 MobilityTensor::diagonalize_tensor(const Tensor3x3& AA,const Vector3d& center){
	double v[3];double R[9];
	for(int i=0;i<3;i++){
		for(int j=i;j<3;j++){
			R[3*i+j]=AA.get(i,j);
		}
	}
	lap::get_eigenvalues_and_vectors_symmetric(3,R,v);
	//eigen vectors
	Vector3d e0(R[3*0+0],R[3*0+1],R[3*0+2]),e1(R[3*1+0],R[3*1+1],R[3*1+2]);
	Vector3d e2=(e0^e1);//keep the right-hand-system
	Tensor3x3 Q=Tensor3x3(e0,e1,e2);rotate(Q,center);
	return Q;
}

//friend function
MobilityTensor operator+(const MobilityTensor& m1,const MobilityTensor& m2){
	return MobilityTensor(
		m1.a+m2.a,m1.b+m2.b,m1.c+m2.c,m1.g+m2.g,m1.h+m2.h,m1.k+m2.k);
}
MobilityTensor operator-(const MobilityTensor& m1,const MobilityTensor& m2){
	return MobilityTensor(
		m1.a-m2.a,m1.b-m2.b,m1.c-m2.c,m1.g-m2.g,m1.h-m2.h,m1.k-m2.k);
}
MobilityTensor operator*(double s, const MobilityTensor &m){
	return MobilityTensor(
		s*m.a,s*m.b,s*m.c,s*m.g,s*m.h,s*m.k);
}
MobilityTensor operator/(const MobilityTensor& m,double s){
	return MobilityTensor(
		m.a/s,m.b/s,m.c/s,m.g/s,m.h/s,m.k/s);
}
