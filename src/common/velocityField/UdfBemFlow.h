//UdfBemFlow
//Stokes�ō쐬�����������ǂݍ���

#ifndef _UDF_BEM_FLOW_H_
#define _UDF_BEM_FLOW_H_

#ifdef _WIN32
#pragma warning (disable : 4786)	// VC++ �Ńf�o�b�K�[�Ɋւ���s�v�̌x������������
#endif


#include "VelocityField.h"
#include "udfmanager.h"
#include "../mesh/StokesMesh2d.h"
#include "../mesh/StokesMesh.h"
#include "../integration/LegendreGaussFormula1d.h"
#include "../integration/LegendreGaussFormulaTriangle.h"
#include "../GreenFunc/StokesGreenFunc2d.h"
#include "../GreenFunc/StokesGreenFunc.h"
#include <string>
using namespace std;

class UdfBemFlow2d:public VelocityField{
public:
	UdfBemFlow2d(const string& udfFileName,double length_unit,double time_unit,double eta_unit);
	virtual ~UdfBemFlow2d();
	virtual void initial(){readData();};//������
	virtual void update(){};//���s
	virtual Vector3d getVelocity(const Vector3d& coord);//��̑��x
	virtual Vector3d getAngularVelocity(const Vector3d& coord);//��̊p���x
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);//��̕ό`���x�e���\��
protected:
	virtual Vector2d reset_coord(const Vector2d& r);
	virtual void readData();
	UDFManager* udf;
private:
	StokesMesh2d* mesh;
	ParameterOfLegendreGaussFormula1d *integral;
	double length_ratio;
	double velocity_ratio;
	Vector2d min;
	Vector2d size;
	bool boundary_check;
};

class UdfBemFlow3d:public VelocityField{
public:
	UdfBemFlow3d(const string& udfFileName,double length_unit,double time_unit,double eta_unit);
	virtual ~UdfBemFlow3d();
	virtual void initial(){readData();};//������
	virtual void update(){};//���s
	virtual Vector3d getVelocity(const Vector3d& coord);//��̑��x
	virtual Vector3d getAngularVelocity(const Vector3d& coord);//��̊p���x
	virtual Tensor3x3 getRateOfStrainTensor(const Vector3d& coord);//��̕ό`���x�e���\��
protected:
	virtual Vector3d reset_coord(const Vector3d& r);
	virtual void readData();
	UDFManager* udf;
private:
	StokesMesh* mesh;
	ParameterOfLegendreGaussFormulaTiangle* integral;
	double length_ratio;
	double velocity_ratio;
	Vector3d min;
	Vector3d size;
	bool boundary_check;
};


#endif // _UDF_BEM_FLOW_H_
