//StaticUdfFlow
//Udfで書かれた時間に依存しない流れを読み込んで速度場とする

#ifndef _UDF_FLOW_H_
#define _UDF_FLOW_H_

#ifdef _WIN32
#pragma warning (disable : 4786)	// VC++ でデバッガーに関する不要の警告を押さえる
#endif


#include "VelocityField.h"
#include "udfmanager.h"
#include "../aprxField/AprxVector2dSquareLatticeField2d.h"
#include "../aprxField/AprxVector3dSquareLatticeField3d.h"
#include <string>
using namespace std;

class UdfFlow2d:public VelocityField{
public:
	UdfFlow2d(const string& udfFileName,double length_unit,double time_unit,double eta_unit);
	virtual ~UdfFlow2d();
	virtual void initial(){readData();};//初期化
	virtual void update(){};//実行
	virtual Vector3d getVelocity(const Vector3d& coord);//場の速度
	virtual Vector3d getAngularVelocity(const Vector3d& coord);//場の角速度
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);//場の変形速度テンソル
protected:
	virtual Vector2d reset_coord(const Vector2d& r);
	virtual void readData();
	UDFManager* udf;
private:
	Vector2d min,size;
	Array2d<Vector2d> field;
	AprxVector2dSquareLatticeField2d* avsl;
	double length_ratio;
	double velocity_ratio;
	//double angular_velocity_ratio;
	//unit of other udf
	//double length_unit;
	//double time_unit;
	//double eta_unit;
};

class UdfFlow3d:public VelocityField{
public:
	UdfFlow3d(const string& udfFileName,double length_unit,double time_unit,double eta_unit);
	virtual ~UdfFlow3d();
	virtual void initial(){readData();};//初期化
	virtual void update(){};//実行
	virtual Vector3d getVelocity(const Vector3d& coord);//場の速度
	virtual Vector3d getAngularVelocity(const Vector3d& coord);//場の角速度
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);//場の変形速度テンソル
protected:
	virtual Vector3d reset_coord(const Vector3d& r);
	virtual void readData();
	UDFManager* udf;
private:
	Vector3d min,size;
	Array3d<Vector3d> field;
	AprxVector3dSquareLatticeField3d* avsl;
	double length_ratio;
	double velocity_ratio;
	//double angular_velocity_ratio;
	//unit of other udf
	//double length_unit;
	//double time_unit;
	//double eta_unit;
};

class DynamicUdfFlow2d:public VelocityField{
public:
	DynamicUdfFlow2d(const string& udfFileName,double length_unit,double time_unit,double eta_unit,double dt);
	virtual ~DynamicUdfFlow2d();
	virtual void initial();//初期化
	virtual void update();//実行
	virtual Vector3d getVelocity(const Vector3d& coord);//場の速度
	virtual Vector3d getAngularVelocity(const Vector3d& coord);//場の角速度
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);//場の変形速度テンソル
protected:
	virtual void readData(int record,Array2d<Vector2d>* field);
	virtual Vector2d reset_coord(const Vector2d& r);
	UDFManager* udf;
private:
	Vector2d min,size;
	//bool order;
	int region;
	Array2d<Vector2d> field1,field2;
	AprxVector2dSquareLatticeField2d *avsl1,*avsl2;
	int timeStepsPerRecord;//レコードあたりの時間ステップ
	double length_ratio;
	double velocity_ratio;
	double dt;
	//double angular_velocity_ratio;
	//unit of other udf
	//double length_unit;
	//double time_unit;
	//double eta_unit;
};

class DynamicUdfFlow3d:public VelocityField{
public:
	DynamicUdfFlow3d(const string& udfFileName,double length_unit,double time_unit,double eta_unit,double dt);
	virtual ~DynamicUdfFlow3d();
	virtual void initial();//初期化
	virtual void update();//実行
	virtual Vector3d getVelocity(const Vector3d& coord);//場の速度
	virtual Vector3d getAngularVelocity(const Vector3d& coord);//場の角速度
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);//場の変形速度テンソル
protected:
	virtual void readData(int record,Array3d<Vector3d>* field);
	virtual Vector3d reset_coord(const Vector3d& r);
	UDFManager* udf;
private:
	Vector3d min,size;
	//bool order;
	int region;
	Array3d<Vector3d> field1,field2;
	AprxVector3dSquareLatticeField3d *avsl1,*avsl2;
	int timeStepsPerRecord;//レコードあたりの時間ステップ
	double power;
	double length_ratio;
	double velocity_ratio;
	double dt;
	//double angular_velocity_ratio;
	//unit of other udf
	//double length_unit;
	//double time_unit;
	//double eta_unit;
};

//以下は時間方向に線形補間されていないバージョン
//class DynamicUdfFlow2d:public UdfFlow2d{
//public:
//	DynamicUdfFlow2d(const string& udfFileName);
//	virtual ~DynamicUdfFlow2d(){};
//	virtual void initial();//初期化
//	virtual void update();//実行
//private:
//	int timePerRecord;//レコードあたりのtime step
//};
//
//class DynamicUdfFlow3d:public UdfFlow3d{
//public:
//	DynamicUdfFlow3d(const string& udfFileName);
//	virtual ~DynamicUdfFlow3d(){};
//	virtual void initial();//初期化
//	virtual void update();//実行
//private:
//	int timePerRecord;//レコードあたりのtime step
//};


#endif // _UDF_FLOW_H_
