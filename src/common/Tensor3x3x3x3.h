//Tensor3x3x3x3

#ifndef _TENSOR_3X3X3X3_H_
#define _TENSOR_3X3X3X3_H_
#include "Tensor3x3x3.h"
#include <iomanip>

class Tensor3x3x3x3{
public:
	Tensor3x3x3x3(
		double xxxx=0.0,double xxxy=0.0,double xxxz=0.0,
		double xxyx=0.0,double xxyy=0.0,double xxyz=0.0,
		double xxzx=0.0,double xxzy=0.0,double xxzz=0.0,//
		double xyxx=0.0,double xyxy=0.0,double xyxz=0.0,
		double xyyx=0.0,double xyyy=0.0,double xyyz=0.0,
		double xyzx=0.0,double xyzy=0.0,double xyzz=0.0,//
		double xzxx=0.0,double xzxy=0.0,double xzxz=0.0,
		double xzyx=0.0,double xzyy=0.0,double xzyz=0.0,
		double xzzx=0.0,double xzzy=0.0,double xzzz=0.0,//
		
		double yxxx=0.0,double yxxy=0.0,double yxxz=0.0,
		double yxyx=0.0,double yxyy=0.0,double yxyz=0.0,
		double yxzx=0.0,double yxzy=0.0,double yxzz=0.0,//
		double yyxx=0.0,double yyxy=0.0,double yyxz=0.0,
		double yyyx=0.0,double yyyy=0.0,double yyyz=0.0,
		double yyzx=0.0,double yyzy=0.0,double yyzz=0.0,//
		double yzxx=0.0,double yzxy=0.0,double yzxz=0.0,
		double yzyx=0.0,double yzyy=0.0,double yzyz=0.0,
		double yzzx=0.0,double yzzy=0.0,double yzzz=0.0,//
		
		double zxxx=0.0,double zxxy=0.0,double zxxz=0.0,
		double zxyx=0.0,double zxyy=0.0,double zxyz=0.0,
		double zxzx=0.0,double zxzy=0.0,double zxzz=0.0,//
		double zyxx=0.0,double zyxy=0.0,double zyxz=0.0,
		double zyyx=0.0,double zyyy=0.0,double zyyz=0.0,
		double zyzx=0.0,double zyzy=0.0,double zyzz=0.0,//
		double zzxx=0.0,double zzxy=0.0,double zzxz=0.0,
		double zzyx=0.0,double zzyy=0.0,double zzyz=0.0,
		double zzzx=0.0,double zzzy=0.0,double zzzz=0.0		
	)
	:x(Tensor3x3x3(xxxx,xxxy,xxxz,xxyx,xxyy,xxyz,xxzx,xxzy,xxzz,//
		             xyxx,xyxy,xyxz,xyyx,xyyy,xyyz,xyzx,xyzy,xyzz,//
		             xzxx,xzxy,xzxz,xzyx,xzyy,xzyz,xzzx,xzzy,xzzz)),
	 y(Tensor3x3x3(yxxx,yxxy,yxxz,yxyx,yxyy,yxyz,yxzx,yxzy,yxzz,//
		             yyxx,yyxy,yyxz,yyyx,yyyy,yyyz,yyzx,yyzy,yyzz,//
		             yzxx,yzxy,yzxz,yzyx,yzyy,yzyz,yzzx,yzzy,yzzz)),
	 z(Tensor3x3x3(zxxx,zxxy,zxxz,zxyx,zxyy,zxyz,zxzx,zxzy,zxzz,//
	               zyxx,zyxy,zyxz,zyyx,zyyy,zyyz,zyzx,zyzy,zyzz,//
	               zzxx,zzxy,zzxz,zzyx,zzyy,zzyz,zzzx,zzzy,zzzz)){}
	Tensor3x3x3x3(const Tensor3x3x3& mx,
	              const Tensor3x3x3& my,
	              const Tensor3x3x3& mz)
	:x(mx),y(my),z(mz){}
	Tensor3x3x3 x,y,z;
	void clear(){
		x.clear();y.clear();z.clear();
	}
	double& operator()(unsigned i,unsigned j,unsigned k,unsigned l);
	double get(unsigned i,unsigned j,unsigned k,unsigned l)const;
	vector<vector<vector<vector<double> > > > getSTLVector();
	Tensor3x3x3x3& operator = (const Tensor3x3x3x3& m);
	Tensor3x3x3x3& operator += (const Tensor3x3x3x3& m);
	Tensor3x3x3x3& operator -= (const Tensor3x3x3x3& m);
	Tensor3x3x3x3& operator *= (const double& m);
	Tensor3x3x3x3& operator /= (const double& m);
	Tensor3x3x3x3 operator - () const;

	friend Tensor3x3x3x3 operator + (const Tensor3x3x3x3& m1,const Tensor3x3x3x3& m2);
	friend Tensor3x3x3x3 operator - (const Tensor3x3x3x3& m1,const Tensor3x3x3x3& m2);
	friend Tensor3x3x3x3 operator / (const Tensor3x3x3x3& m,const double& s);
	friend Tensor3x3x3x3 operator * (const double& d,const Tensor3x3x3x3& m1);
	friend Tensor3x3x3x3 operator * (const Tensor3x3x3x3& m1,const double& d);
	
	//����ȑ���
	//(A�EB)_{ijkl}=A_{ijkm}B_{ml}�̑����s��
	friend Tensor3x3x3x3 operator * (const Tensor3x3x3x3& m1,const Tensor3x3& m2);
	friend Tensor3x3x3x3 operator * (const Tensor3x3& m1,const Tensor3x3x3x3& m2);
	//(A�EB)_{ijkl}=A_{ijn}B_{nkl}�̑����s��
	friend Tensor3x3x3x3 operator * (const Tensor3x3x3& A,const Tensor3x3x3& B);
	//(A:B)_{ij}=A_{ijkl}B_{lk}�̂�s��
	friend Tensor3x3  dot2(const Tensor3x3x3x3& m1, const Tensor3x3& m2);
	friend Tensor3x3  dot2(const Tensor3x3& m1, const Tensor3x3x3x3& m2);
	//(AB)_{ijkl}=A_{ijk}B_{l}
	friend Tensor3x3x3x3 dyad(const Tensor3x3x3& m,const Vector3d& v);
	//(BA)_{ijkl}=B_{i}A_{jkl}
	friend Tensor3x3x3x3 dyad(const Vector3d& v,const Tensor3x3x3& m);
	//(AB)_{ijkl}=A_{ij}B_{kl}
	friend Tensor3x3x3x3 dyad(const Tensor3x3& m1,const Tensor3x3& m2);
	//(AB)_{ijkl}=A_{ijl}B_{k}�̑����s��
	friend Tensor3x3x3x3 diamond(const Tensor3x3x3& m,const Vector3d& v);
	//(BA)_{ijkl}=B_{j}A_{ikl}�̑����s��
	friend Tensor3x3x3x3 diamond(const Vector3d& v,const Tensor3x3x3& m);
	//(AB)_{ijkl}=A_{ik}*B_{jl}
	friend Tensor3x3x3x3 diamond(const Tensor3x3& m1,const Tensor3x3& m2);
	//(AB)_{ijkl}=A_{il}*B_{jk}
	friend Tensor3x3x3x3 diamond2(const Tensor3x3& m1,const Tensor3x3& m2);
	//���W�n�̕ϊ��̃e���\��Q��g���č��W�n��ϊ�
	//G'_{abc}=Q_{ar}Q_{bp}Q_{cq}G_{rpq}
	friend Tensor3x3x3x3 convert(const Tensor3x3& Q,const Tensor3x3x3x3& G);
	
	friend ostream& operator << (ostream& os, const Tensor3x3x3x3& m)
	{
		os	<< m.x << endl << m.y << endl << m.z;
		return os;
	}

	friend istream& operator >> (istream& is, Tensor3x3x3x3& m)
	{
		is	>> m.x >> m.y >> m.z;
		return is;
	}
};

inline double& Tensor3x3x3x3::operator()(unsigned i,unsigned j,unsigned k, unsigned l)
{
	if(i==0){
		return x(j,k,l);
	}else if(i==1){
		return y(j,k,l);
	}else{
		return z(j,k,l);
	}
}

inline double Tensor3x3x3x3::get(unsigned i,unsigned j,unsigned k, unsigned l)const
{
	if(i==0){
		return x.get(j,k,l);
	}else if(i==1){
		return y.get(j,k,l);
	}else{
		return z.get(j,k,l);
	}
}

inline vector<vector<vector<vector<double> > > > Tensor3x3x3x3::getSTLVector()
{
	vector<vector<vector<vector<double> > > > a;
	a.resize(3);
	a[0]=x.getSTLVector();
	a[1]=y.getSTLVector();
	a[2]=z.getSTLVector();
	return a;
}

inline Tensor3x3x3x3& Tensor3x3x3x3::operator = (const Tensor3x3x3x3& m)
{
	x=m.x;y=m.y;z=m.z;
	return *this;
}

inline Tensor3x3x3x3& Tensor3x3x3x3::operator += (const Tensor3x3x3x3& m)
{
	x+=m.x;y+=m.y;z+=m.z;
	return *this;
}

inline Tensor3x3x3x3& Tensor3x3x3x3::operator -= (const Tensor3x3x3x3& m)
{
	x-=m.x;y-=m.y;z-=m.z;
	return *this;
}

inline Tensor3x3x3x3& Tensor3x3x3x3::operator *= (const double& c)
{
	x*=c;y*=c;z*=c;
	return *this;
}

inline Tensor3x3x3x3& Tensor3x3x3x3::operator /= (const double& c)
{
	x/=c;y/=c;z/=c;
	return *this;
}

inline Tensor3x3x3x3 Tensor3x3x3x3::operator - () const
{
	return Tensor3x3x3x3(-x, -y, -z);
}

inline Tensor3x3x3x3 operator + (const Tensor3x3x3x3& m1, const Tensor3x3x3x3& m2)
{
	return Tensor3x3x3x3(m1.x+m2.x,m1.y+m2.y,m1.z+m2.z);
}

inline Tensor3x3x3x3 operator - (const Tensor3x3x3x3& m1, const Tensor3x3x3x3& m2)
{
	return Tensor3x3x3x3(m1.x-m2.x,m1.y-m2.y,m1.z-m2.z);
}

inline Tensor3x3x3x3 operator / (const Tensor3x3x3x3& m, const double& s)
{
	return Tensor3x3x3x3(m.x/s,m.y/s,m.z/s);
}

inline Tensor3x3x3x3 operator * (const double& s, const Tensor3x3x3x3& m)
{
	return Tensor3x3x3x3(s*m.x,s*m.y,s*m.z);
}

inline Tensor3x3x3x3 operator * (const Tensor3x3x3x3& m, const double& s)
{
	return Tensor3x3x3x3(m.x*s,m.y*s,m.z*s);
}

inline Tensor3x3x3x3 operator * (const Tensor3x3x3x3& m1,const Tensor3x3& m2)
{
	return Tensor3x3x3x3(m1.x*m2,m1.y*m2,m1.z*m2);
}

inline Tensor3x3x3x3 operator * (const Tensor3x3& m1,const Tensor3x3x3x3& m2)
{
	return Tensor3x3x3x3(
		m1.x.x*m2.x+m1.x.y*m2.y+m1.x.z*m2.z,
		m1.y.x*m2.x+m1.y.y*m2.y+m1.y.z*m2.z,
		m1.z.x*m2.x+m1.z.y*m2.y+m1.z.z*m2.z
	);
}

inline Tensor3x3x3x3 operator * (const Tensor3x3x3& A,const Tensor3x3x3& B)
{
	return Tensor3x3x3x3(A.x*B,A.y*B,A.z*B);
}

inline Tensor3x3 dot2(const Tensor3x3x3x3& m1, const Tensor3x3& m2){
	return Tensor3x3(dot2(m1.x,m2),dot2(m1.y,m2),dot2(m1.z,m2));
}

inline Tensor3x3 dot2(const Tensor3x3& m1, const Tensor3x3x3x3& m2){
	Tensor3x3x3x3 a=m1*m2;
	return Tensor3x3(
		a.x.x.x.x+a.y.y.x.x+a.z.z.x.x, a.x.x.x.y+a.y.y.x.y+a.z.z.x.y, a.x.x.x.z+a.y.y.x.z+a.z.z.x.z,
		a.x.x.y.x+a.y.y.y.x+a.z.z.y.x, a.x.x.y.y+a.y.y.y.y+a.z.z.y.y, a.x.x.y.z+a.y.y.y.z+a.z.z.y.z,
		a.x.x.z.x+a.y.y.z.x+a.z.z.z.x, a.x.x.z.y+a.y.y.z.y+a.z.z.z.y, a.x.x.z.z+a.y.y.z.z+a.z.z.z.z
	);
}

inline Tensor3x3x3x3 dyad(const Tensor3x3x3& m,const Vector3d& v)
{
	return Tensor3x3x3x3(
		m.x.x.x*v.x, m.x.x.x*v.y, m.x.x.x*v.z,
		m.x.x.y*v.x, m.x.x.y*v.y, m.x.x.y*v.z,
		m.x.x.z*v.x, m.x.x.z*v.y, m.x.x.z*v.z,//
		m.x.y.x*v.x, m.x.y.x*v.y, m.x.y.x*v.z,
		m.x.y.y*v.x, m.x.y.y*v.y, m.x.y.y*v.z,
		m.x.y.z*v.x, m.x.y.z*v.y, m.x.y.z*v.z,//
		m.x.z.x*v.x, m.x.z.x*v.y, m.x.z.x*v.z,
		m.x.z.y*v.x, m.x.z.y*v.y, m.x.z.y*v.z,
		m.x.z.z*v.x, m.x.z.z*v.y, m.x.z.z*v.z,//
	
		m.y.x.x*v.x, m.y.x.x*v.y, m.y.x.x*v.z,
		m.y.x.y*v.x, m.y.x.y*v.y, m.y.x.y*v.z,
		m.y.x.z*v.x, m.y.x.z*v.y, m.y.x.z*v.z,//
		m.y.y.x*v.x, m.y.y.x*v.y, m.y.y.x*v.z,
		m.y.y.y*v.x, m.y.y.y*v.y, m.y.y.y*v.z,
		m.y.y.z*v.x, m.y.y.z*v.y, m.y.y.z*v.z,//
		m.y.z.x*v.x, m.y.z.x*v.y, m.y.z.x*v.z,
		m.y.z.y*v.x, m.y.z.y*v.y, m.y.z.y*v.z,
		m.y.z.z*v.x, m.y.z.z*v.y, m.y.z.z*v.z,//
		
		m.z.x.x*v.x, m.z.x.x*v.y, m.z.x.x*v.z,
		m.z.x.y*v.x, m.z.x.y*v.y, m.z.x.y*v.z,
		m.z.x.z*v.x, m.z.x.z*v.y, m.z.x.z*v.z,//
		m.z.y.x*v.x, m.z.y.x*v.y, m.z.y.x*v.z,
		m.z.y.y*v.x, m.z.y.y*v.y, m.z.y.y*v.z,
		m.z.y.z*v.x, m.z.y.z*v.y, m.z.y.z*v.z,//
		m.z.z.x*v.x, m.z.z.x*v.y, m.z.z.x*v.z,
		m.z.z.y*v.x, m.z.z.y*v.y, m.z.z.y*v.z,
		m.z.z.z*v.x, m.z.z.z*v.y, m.z.z.z*v.z
	);
}

inline Tensor3x3x3x3 dyad(const Vector3d& v,const Tensor3x3x3& m)
{
	return Tensor3x3x3x3(
		v.x*m.x.x.x, v.x*m.x.x.y, v.x*m.x.x.z,
		v.x*m.x.y.x, v.x*m.x.y.y, v.x*m.x.y.z,
		v.x*m.x.z.x, v.x*m.x.z.y, v.x*m.x.z.z,//
		v.x*m.y.x.x, v.x*m.y.x.y, v.x*m.y.x.z,
		v.x*m.y.y.x, v.x*m.y.y.y, v.x*m.y.y.z,
		v.x*m.y.z.x, v.x*m.y.z.y, v.x*m.y.z.z,//
		v.x*m.z.x.x, v.x*m.z.x.y, v.x*m.z.x.z,
		v.x*m.z.y.x, v.x*m.z.y.y, v.x*m.z.y.z,
		v.x*m.z.z.x, v.x*m.z.z.y, v.x*m.z.z.z,//
		
		v.y*m.x.x.x, v.y*m.x.x.y, v.y*m.x.x.z,
		v.y*m.x.y.x, v.y*m.x.y.y, v.y*m.x.y.z,
		v.y*m.x.z.x, v.y*m.x.z.y, v.y*m.x.z.z,//
		v.y*m.y.x.x, v.y*m.y.x.y, v.y*m.y.x.z,
		v.y*m.y.y.x, v.y*m.y.y.y, v.y*m.y.y.z,
		v.y*m.y.z.x, v.y*m.y.z.y, v.y*m.y.z.z,//
		v.y*m.z.x.x, v.y*m.z.x.y, v.y*m.z.x.z,
		v.y*m.z.y.x, v.y*m.z.y.y, v.y*m.z.y.z,
		v.y*m.z.z.x, v.y*m.z.z.y, v.y*m.z.z.z,//
		
		v.z*m.x.x.x, v.z*m.x.x.y, v.z*m.x.x.z,
		v.z*m.x.y.x, v.z*m.x.y.y, v.z*m.x.y.z,
		v.z*m.x.z.x, v.z*m.x.z.y, v.z*m.x.z.z,//
		v.z*m.y.x.x, v.z*m.y.x.y, v.z*m.y.x.z,
		v.z*m.y.y.x, v.z*m.y.y.y, v.z*m.y.y.z,
		v.z*m.y.z.x, v.z*m.y.z.y, v.z*m.y.z.z,//
		v.z*m.z.x.x, v.z*m.z.x.y, v.z*m.z.x.z,
		v.z*m.z.y.x, v.z*m.z.y.y, v.z*m.z.y.z,
		v.z*m.z.z.x, v.z*m.z.z.y, v.z*m.z.z.z
	);
}

inline Tensor3x3x3x3 dyad(const Tensor3x3& m1,const Tensor3x3& m2)
{
	return Tensor3x3x3x3(
		m1.x.x*m2.x.x,m1.x.x*m2.x.y,m1.x.x*m2.x.z,
		m1.x.x*m2.y.x,m1.x.x*m2.y.y,m1.x.x*m2.y.z,
		m1.x.x*m2.z.x,m1.x.x*m2.z.y,m1.x.x*m2.z.z,//
		m1.x.y*m2.x.x,m1.x.y*m2.x.y,m1.x.y*m2.x.z,
		m1.x.y*m2.y.x,m1.x.y*m2.y.y,m1.x.y*m2.y.z,
		m1.x.y*m2.z.x,m1.x.y*m2.z.y,m1.x.y*m2.z.z,//
		m1.x.z*m2.x.x,m1.x.z*m2.x.y,m1.x.z*m2.x.z,
		m1.x.z*m2.y.x,m1.x.z*m2.y.y,m1.x.z*m2.y.z,
		m1.x.z*m2.z.x,m1.x.z*m2.z.y,m1.x.z*m2.z.z,//
		
		m1.y.x*m2.x.x,m1.y.x*m2.x.y,m1.y.x*m2.x.z,
		m1.y.x*m2.y.x,m1.y.x*m2.y.y,m1.y.x*m2.y.z,
		m1.y.x*m2.z.x,m1.y.x*m2.z.y,m1.y.x*m2.z.z,//
		m1.y.y*m2.x.x,m1.y.y*m2.x.y,m1.y.y*m2.x.z,
		m1.y.y*m2.y.x,m1.y.y*m2.y.y,m1.y.y*m2.y.z,
		m1.y.y*m2.z.x,m1.y.y*m2.z.y,m1.y.y*m2.z.z,//
		m1.y.z*m2.x.x,m1.y.z*m2.x.y,m1.y.z*m2.x.z,
		m1.y.z*m2.y.x,m1.y.z*m2.y.y,m1.y.z*m2.y.z,
		m1.y.z*m2.z.x,m1.y.z*m2.z.y,m1.y.z*m2.z.z,//
		
		m1.z.x*m2.x.x,m1.z.x*m2.x.y,m1.z.x*m2.x.z,
		m1.z.x*m2.y.x,m1.z.x*m2.y.y,m1.z.x*m2.y.z,
		m1.z.x*m2.z.x,m1.z.x*m2.z.y,m1.z.x*m2.z.z,//
		m1.z.y*m2.x.x,m1.z.y*m2.x.y,m1.z.y*m2.x.z,
		m1.z.y*m2.y.x,m1.z.y*m2.y.y,m1.z.y*m2.y.z,
		m1.z.y*m2.z.x,m1.z.y*m2.z.y,m1.z.y*m2.z.z,//
		m1.z.z*m2.x.x,m1.z.z*m2.x.y,m1.z.z*m2.x.z,
		m1.z.z*m2.y.x,m1.z.z*m2.y.y,m1.z.z*m2.y.z,
		m1.z.z*m2.z.x,m1.z.z*m2.z.y,m1.z.z*m2.z.z//
	);
}

inline Tensor3x3x3x3 diamond(const Tensor3x3x3& m,const Vector3d& v)
{
	return Tensor3x3x3x3(
		m.x.x.x*v.x, m.x.x.y*v.x, m.x.x.z*v.x,
		m.x.x.x*v.y, m.x.x.y*v.y, m.x.x.z*v.y,
		m.x.x.x*v.z, m.x.x.y*v.z, m.x.x.z*v.z,//
		m.x.y.x*v.x, m.x.y.y*v.x, m.x.y.z*v.x,
		m.x.y.x*v.y, m.x.y.y*v.y, m.x.y.z*v.y,
		m.x.y.x*v.z, m.x.y.y*v.z, m.x.y.z*v.z,//
		m.x.z.x*v.x, m.x.z.y*v.x, m.x.z.z*v.x,
		m.x.z.x*v.y, m.x.z.y*v.y, m.x.z.z*v.y,
		m.x.z.x*v.z, m.x.z.y*v.z, m.x.z.z*v.z,//
		
		m.y.x.x*v.x, m.y.x.y*v.x, m.y.x.z*v.x,
		m.y.x.x*v.y, m.y.x.y*v.y, m.y.x.z*v.y,
		m.y.x.x*v.z, m.y.x.y*v.z, m.y.x.z*v.z,//
		m.y.y.x*v.x, m.y.y.y*v.x, m.y.y.z*v.x,
		m.y.y.x*v.y, m.y.y.y*v.y, m.y.y.z*v.y,
		m.y.y.x*v.z, m.y.y.y*v.z, m.y.y.z*v.z,//
		m.y.z.x*v.x, m.y.z.y*v.x, m.y.z.z*v.x,
		m.y.z.x*v.y, m.y.z.y*v.y, m.y.z.z*v.y,
		m.y.z.x*v.z, m.y.z.y*v.z, m.y.z.z*v.z,//
		
		m.z.x.x*v.x, m.z.x.y*v.x, m.z.x.z*v.x,
		m.z.x.x*v.y, m.z.x.y*v.y, m.z.x.z*v.y,
		m.z.x.x*v.z, m.z.x.y*v.z, m.z.x.z*v.z,//
		m.z.y.x*v.x, m.z.y.y*v.x, m.z.y.z*v.x,
		m.z.y.x*v.y, m.z.y.y*v.y, m.z.y.z*v.y,
		m.z.y.x*v.z, m.z.y.y*v.z, m.z.y.z*v.z,//
		m.z.z.x*v.x, m.z.z.y*v.x, m.z.z.z*v.x,
		m.z.z.x*v.y, m.z.z.y*v.y, m.z.z.z*v.y,
		m.z.z.x*v.z, m.z.z.y*v.z, m.z.z.z*v.z
	);
}

inline Tensor3x3x3x3 diamond(const Vector3d& v,const Tensor3x3x3& m)
{
	return Tensor3x3x3x3(
		m.x.x.x*v.x, m.x.x.y*v.x, m.x.x.z*v.x,
		m.x.y.x*v.x, m.x.y.y*v.x, m.x.y.z*v.x,
		m.x.z.x*v.x, m.x.z.y*v.x, m.x.z.z*v.x,//
		m.x.x.x*v.y, m.x.x.y*v.y, m.x.x.z*v.y,
		m.x.y.x*v.y, m.x.y.y*v.y, m.x.y.z*v.y,
		m.x.z.x*v.y, m.x.z.y*v.y, m.x.z.z*v.y,//
		m.x.x.x*v.z, m.x.x.y*v.z, m.x.x.z*v.z,
		m.x.y.x*v.z, m.x.y.y*v.z, m.x.y.z*v.z,
		m.x.z.x*v.z, m.x.z.y*v.z, m.x.z.z*v.z,//
		
		m.y.x.x*v.x, m.y.x.y*v.x, m.y.x.z*v.x,
		m.y.y.x*v.x, m.y.y.y*v.x, m.y.y.z*v.x,
		m.y.z.x*v.x, m.y.z.y*v.x, m.y.z.z*v.x,//
		m.y.x.x*v.y, m.y.x.y*v.y, m.y.x.z*v.y,
		m.y.y.x*v.y, m.y.y.y*v.y, m.y.y.z*v.y,
		m.y.z.x*v.y, m.y.z.y*v.y, m.y.z.z*v.y,//
		m.y.x.x*v.z, m.y.x.y*v.z, m.y.x.z*v.z,
		m.y.y.x*v.z, m.y.y.y*v.z, m.y.y.z*v.z,
		m.y.z.x*v.z, m.y.z.y*v.z, m.y.z.z*v.z,//
		
		m.z.x.x*v.x, m.z.x.y*v.x, m.z.x.z*v.x,
		m.z.y.x*v.x, m.z.y.y*v.x, m.z.y.z*v.x,
		m.z.z.x*v.x, m.z.z.y*v.x, m.z.z.z*v.x,//
		m.z.x.x*v.y, m.z.x.y*v.y, m.z.x.z*v.y,
		m.z.y.x*v.y, m.z.y.y*v.y, m.z.y.z*v.y,
		m.z.z.x*v.y, m.z.z.y*v.y, m.z.z.z*v.y,//
		m.z.x.x*v.z, m.z.x.y*v.z, m.z.x.z*v.z,
		m.z.y.x*v.z, m.z.y.y*v.z, m.z.y.z*v.z,
		m.z.z.x*v.z, m.z.z.y*v.z, m.z.z.z*v.z//
	);
}

inline Tensor3x3x3x3 diamond(const Tensor3x3& m1,const Tensor3x3& m2){
	return Tensor3x3x3x3(
		m1.x.x*m2.x.x,m1.x.x*m2.x.y,m1.x.x*m2.x.z,
		m1.x.y*m2.x.x,m1.x.y*m2.x.y,m1.x.y*m2.x.z,
		m1.x.z*m2.x.x,m1.x.z*m2.x.y,m1.x.z*m2.x.z,//
		m1.x.x*m2.y.x,m1.x.x*m2.y.y,m1.x.x*m2.y.z,
		m1.x.y*m2.y.x,m1.x.y*m2.y.y,m1.x.y*m2.y.z,
		m1.x.z*m2.y.x,m1.x.z*m2.y.y,m1.x.z*m2.y.z,//
		m1.x.x*m2.z.x,m1.x.x*m2.z.y,m1.x.x*m2.z.z,
		m1.x.y*m2.z.x,m1.x.y*m2.z.y,m1.x.y*m2.z.z,
		m1.x.z*m2.z.x,m1.x.z*m2.z.y,m1.x.z*m2.z.z,//
		
		m1.y.x*m2.x.x,m1.y.x*m2.x.y,m1.y.x*m2.x.z,
		m1.y.y*m2.x.x,m1.y.y*m2.x.y,m1.y.y*m2.x.z,
		m1.y.z*m2.x.x,m1.y.z*m2.x.y,m1.y.z*m2.x.z,//
		m1.y.x*m2.y.x,m1.y.x*m2.y.y,m1.y.x*m2.y.z,
		m1.y.y*m2.y.x,m1.y.y*m2.y.y,m1.y.y*m2.y.z,
		m1.y.z*m2.y.x,m1.y.z*m2.y.y,m1.y.z*m2.y.z,//
		m1.y.x*m2.z.x,m1.y.x*m2.z.y,m1.y.x*m2.z.z,
		m1.y.y*m2.z.x,m1.y.y*m2.z.y,m1.y.y*m2.z.z,
		m1.y.z*m2.z.x,m1.y.z*m2.z.y,m1.y.z*m2.z.z,//
		
		m1.z.x*m2.x.x,m1.z.x*m2.x.y,m1.z.x*m2.x.z,
		m1.z.y*m2.x.x,m1.z.y*m2.x.y,m1.z.y*m2.x.z,
		m1.z.z*m2.x.x,m1.z.z*m2.x.y,m1.z.z*m2.x.z,//
		m1.z.x*m2.y.x,m1.z.x*m2.y.y,m1.z.x*m2.y.z,
		m1.z.y*m2.y.x,m1.z.y*m2.y.y,m1.z.y*m2.y.z,
		m1.z.z*m2.y.x,m1.z.z*m2.y.y,m1.z.z*m2.y.z,//
		m1.z.x*m2.z.x,m1.z.x*m2.z.y,m1.z.x*m2.z.z,
		m1.z.y*m2.z.x,m1.z.y*m2.z.y,m1.z.y*m2.z.z,
		m1.z.z*m2.z.x,m1.z.z*m2.z.y,m1.z.z*m2.z.z//
	);
}
	//(AB)_{ijkl}=A_{il}*B_{jk}
inline Tensor3x3x3x3 diamond2(const Tensor3x3& m1,const Tensor3x3& m2){
	return Tensor3x3x3x3(
		m1.x.x*m2.x.x,m1.x.x*m2.y.x,m1.x.x*m2.z.x,
		m1.x.y*m2.x.x,m1.x.y*m2.y.x,m1.x.y*m2.z.x,
		m1.x.z*m2.x.x,m1.x.z*m2.y.x,m1.x.z*m2.z.x,//
		m1.x.x*m2.x.y,m1.x.x*m2.y.y,m1.x.x*m2.z.y,
		m1.x.y*m2.x.y,m1.x.y*m2.y.y,m1.x.y*m2.z.y,
		m1.x.z*m2.x.y,m1.x.z*m2.y.y,m1.x.z*m2.z.y,//
		m1.x.x*m2.z.x,m1.x.x*m2.y.z,m1.x.x*m2.z.z,
		m1.x.y*m2.x.z,m1.x.y*m2.y.z,m1.x.y*m2.z.z,
		m1.x.z*m2.x.z,m1.x.z*m2.y.z,m1.x.z*m2.z.z,//
		
		m1.y.x*m2.x.x,m1.y.x*m2.y.z,m1.y.x*m2.z.x,
		m1.y.y*m2.x.x,m1.y.y*m2.y.x,m1.y.y*m2.z.x,
		m1.y.z*m2.x.x,m1.y.z*m2.y.x,m1.y.z*m2.z.x,//
		m1.y.x*m2.x.y,m1.y.x*m2.y.y,m1.y.x*m2.z.y,
		m1.y.y*m2.x.y,m1.y.y*m2.y.y,m1.y.y*m2.z.y,
		m1.y.z*m2.x.y,m1.y.z*m2.y.y,m1.y.z*m2.z.y,//
		m1.y.x*m2.x.z,m1.y.x*m2.y.z,m1.y.x*m2.z.z,
		m1.y.y*m2.x.z,m1.y.y*m2.y.z,m1.y.y*m2.z.z,
		m1.y.z*m2.x.z,m1.y.z*m2.y.z,m1.y.z*m2.z.z,//
		
		m1.z.x*m2.x.x,m1.z.x*m2.y.x,m1.z.x*m2.z.x,
		m1.z.y*m2.x.x,m1.z.y*m2.y.x,m1.z.y*m2.z.x,
		m1.z.z*m2.x.x,m1.z.z*m2.y.x,m1.z.z*m2.z.x,//
		m1.z.x*m2.x.y,m1.z.x*m2.y.y,m1.z.x*m2.z.y,
		m1.z.y*m2.x.y,m1.z.y*m2.y.y,m1.z.y*m2.z.y,
		m1.z.z*m2.x.y,m1.z.z*m2.y.y,m1.z.z*m2.z.y,//
		m1.z.x*m2.x.z,m1.z.x*m2.y.z,m1.z.x*m2.z.z,
		m1.z.y*m2.x.z,m1.z.y*m2.y.z,m1.z.y*m2.z.z,
		m1.z.z*m2.x.z,m1.z.z*m2.y.z,m1.z.z*m2.z.z//
	);
}

inline Tensor3x3x3x3 convert(const Tensor3x3& Q,const Tensor3x3x3x3& G)
{
	Tensor3x3x3 tx=convert(Q,G.x);
	Tensor3x3x3 ty=convert(Q,G.y);
	Tensor3x3x3 tz=convert(Q,G.z);
	return Q*Tensor3x3x3x3(tx,ty,tz);
}

#endif //_TENSOR_3X3X3X3_H_

