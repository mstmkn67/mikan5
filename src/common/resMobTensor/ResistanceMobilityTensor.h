#ifndef _RESISTANCE_MOBILITY_TENSOR_H_
#define _RESISTANCE_MOBILITY_TENSOR_H_

#include "TensorUtility.h"
//#include "../solver/SymmetricEigenSystem.h"
#include "../solver/LapackFunctions.h"
class MobilityTensor;

class ResistanceTensor{
public:
	ResistanceTensor(){};
	ResistanceTensor(const ResistanceTensor& r);
	ResistanceTensor(const MobilityTensor& m);
	ResistanceTensor(const Tensor3x3& A,const Tensor3x3& B,const Tensor3x3& C,
	                 const Tensor3x3x3& G,const Tensor3x3x3& H,const Tensor3x3x3x3& K);
	virtual ~ResistanceTensor(){};
	virtual void clear(){A.clear();B.clear();C.clear();G.clear();H.clear();K.clear();};
	virtual void set(const Tensor3x3& A,const Tensor3x3& B,const Tensor3x3& C,
	                 const Tensor3x3x3& G,const Tensor3x3x3& H,const Tensor3x3x3x3& K);
	virtual MobilityTensor getMobilityTensor()const;
	virtual ResistanceTensor getTranslateTensor(const Vector3d& r)const;//translate of object
	virtual void translate(const Vector3d& r);//translate of object
	virtual void rotate(const Tensor3x3& R);
	virtual void rotate(const Tensor3x3& R,const Vector3d& center);
	virtual Vector3d getCenterOfResistance()const;
	virtual ResistanceTensor& operator=(const ResistanceTensor& r);
	friend ResistanceTensor operator+(const ResistanceTensor& r1,const ResistanceTensor& r2);
	friend ResistanceTensor operator-(const ResistanceTensor& r1,const ResistanceTensor& r2);
	friend ResistanceTensor operator*(double s, const ResistanceTensor &r);
	friend ResistanceTensor operator/(const ResistanceTensor& r,double s);
	virtual ResistanceTensor& operator *= (double s);
	virtual ResistanceTensor& operator /= (double s);
	virtual ResistanceTensor operator - () const;
	
	virtual Vector3d selectCenterOfResistance();//return value is traslate vector
	virtual Vector3d selectCenterOfMobility();//return value is traslate vector
	virtual Tensor3x3 select_A_diagonalizedSystem();//return value is rotate tensor
	virtual Tensor3x3 select_B_diagonalizedSystem();//return value is rotate tensor
	virtual Tensor3x3 select_C_diagonalizedSystem();//return value is rotate tensor
	virtual Tensor3x3 select_a_diagonalizedSystem();//return value is rotate tensor
	virtual Tensor3x3 select_b_diagonalizedSystem();//return value is rotate tensor
	virtual Tensor3x3 select_c_diagonalizedSystem();//return value is rotate tensor
	
	virtual Tensor3x3 get_A()const{return A;};
	virtual Tensor3x3 get_B()const{return B;};
	virtual Tensor3x3 get_C()const{return C;};
	virtual Tensor3x3x3 get_G()const{return G;};
	virtual Tensor3x3x3 get_H()const{return H;};
	virtual Tensor3x3x3x3 get_K()const{return K;};
	Tensor3x3 A,B,C;Tensor3x3x3 G,H;Tensor3x3x3x3 K;
protected:
	virtual Tensor3x3 diagonalize_tensor(const Tensor3x3& R);
	virtual Tensor3x3 diagonalize_tensor(const Tensor3x3& R,const Vector3d& center);
private:
};

class MobilityTensor{
public:
	MobilityTensor(){};
	MobilityTensor(const MobilityTensor& m);
	MobilityTensor(const ResistanceTensor& r);
	MobilityTensor(const Tensor3x3& a,const Tensor3x3& b,const Tensor3x3& c,
	               const Tensor3x3x3& g,const Tensor3x3x3& h,const Tensor3x3x3x3& k);
	virtual ~MobilityTensor(){};
	virtual void clear(){a.clear();b.clear();c.clear();g.clear();h.clear();k.clear();};
	virtual void set(const Tensor3x3& a,const Tensor3x3& b,const Tensor3x3& c,
	                 const Tensor3x3x3& g,const Tensor3x3x3& h,const Tensor3x3x3x3& k);
	virtual ResistanceTensor getResistanceTensor()const;
	virtual MobilityTensor getTranslateTensor(const Vector3d& r)const;//translate of object
	virtual void translate(const Vector3d& r);//translate of object
	virtual void rotate(const Tensor3x3& R);
	virtual void rotate(const Tensor3x3& R,const Vector3d& center);
	virtual Vector3d getCenterOfMobility()const;
	virtual MobilityTensor& operator=(const MobilityTensor& m);
	friend MobilityTensor operator+(const MobilityTensor& m1,const MobilityTensor& m2);
	friend MobilityTensor operator-(const MobilityTensor& m1,const MobilityTensor& m2);
	friend MobilityTensor operator*(double s, const MobilityTensor &m);
	friend MobilityTensor operator/(const MobilityTensor& m,double s);
	virtual MobilityTensor& operator *= (double s);
	virtual MobilityTensor& operator /= (double s);
	virtual MobilityTensor operator - () const;
	
	virtual Vector3d selectCenterOfResistance();//return value is traslate vector
	virtual Vector3d selectCenterOfMobility();//return value is traslate vector
	virtual Tensor3x3 select_A_diagonalizedSystem();//return value is rotate tensor
	virtual Tensor3x3 select_B_diagonalizedSystem();//return value is rotate tensor
	virtual Tensor3x3 select_C_diagonalizedSystem();//return value is rotate tensor
	virtual Tensor3x3 select_a_diagonalizedSystem();//return value is rotate tensor
	virtual Tensor3x3 select_b_diagonalizedSystem();//return value is rotate tensor
	virtual Tensor3x3 select_c_diagonalizedSystem();//return value is rotate tensor
	
	virtual Tensor3x3 get_a()const{return a;};
	virtual Tensor3x3 get_b()const{return b;};
	virtual Tensor3x3 get_c()const{return c;};
	virtual Tensor3x3x3 get_g()const{return g;};
	virtual Tensor3x3x3 get_h()const{return h;};
	virtual Tensor3x3x3x3 get_k()const{return k;};
	Tensor3x3 a,b,c;Tensor3x3x3 g,h;Tensor3x3x3x3 k;
protected:
	virtual Tensor3x3 diagonalize_tensor(const Tensor3x3& R);
	virtual Tensor3x3 diagonalize_tensor(const Tensor3x3& R,const Vector3d& center);
private:
};

#endif // _RESISTANCE_MOBILITY_TENSOR_H_
