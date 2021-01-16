#ifndef _LEGENDRE_GAUSS_FORMULA_1D_H_
#define _LEGENDRE_GAUSS_FORMULA_1D_H_

//base class of 1 dimensional Legendre-Gauss integration formula 
template<class returnType,class argumentType>
class LegendreGaussFormula1d{
public:
	LegendreGaussFormula1d(int point);
	virtual ~LegendreGaussFormula1d();
	virtual returnType integrand(const argumentType& r)const=0;
	virtual returnType integral(const argumentType& rmin,const argumentType& rmax,double length)const;
	double* xi;
	double* w;
	int point;
protected:
};

struct ParameterOfLegendreGaussFormula1d:public LegendreGaussFormula1d<double,double>{
public:
	ParameterOfLegendreGaussFormula1d(int point):LegendreGaussFormula1d<double,double>(point){};
	virtual double integrand(const double& r)const{return 0;};
private:
};

template<class returnType,class argumentType>
LegendreGaussFormula1d<returnType,argumentType>::LegendreGaussFormula1d(int p):point(p){
	xi=new double[point];w=new double[point];
	if(point==1){
		xi[0]=0.0;w[0]=2.0;
	}else if(point==2){
		xi[0]=1.0/sqrt(3.);xi[1]=-xi[0];
		w[0]=w[1]=1.0;
	}else if(point==3){
		xi[0]=0.0;xi[1]=sqrt(0.6);xi[2]=-xi[1];
		w[0]=8.0/9.0;w[1]=w[2]=5.0/9.0;
	}else if(point==4){
		xi[0]=0.861136;xi[1]=-xi[0];xi[2]=0.339981;xi[3]=-xi[2];
		w[0]=w[1]=0.347855;w[2]=w[3]=0.652145;
	}else if(point==5){
		xi[0]=0;xi[1]=0.906180;xi[2]=-xi[1];xi[3]=0.538469;xi[4]=-xi[3];
		w[0]=0.568889;w[1]=w[2]=0.236927;w[3]=w[4]=0.478629;
	}else if(point==6){
		xi[0]=0.932470;xi[1]=-xi[0];xi[2]=0.661209;xi[3]=-xi[2];xi[4]=0.238619;xi[5]=-xi[4];
		w[0]=w[1]=0.171324;w[2]=w[3]=0.360702;w[4]=w[5]=0.467914;
	}else{
		cout << point << " points algorithm is ";
		cout << "not implemented in class LegendreGaussFormula1d" << endl;
		exit(1);
	}
}
template<class returnType,class argumentType>
LegendreGaussFormula1d<returnType,argumentType>::~LegendreGaussFormula1d(){
	delete[] xi;delete[] w;
}

template<class returnType,class argumentType>
returnType LegendreGaussFormula1d<returnType,argumentType>::
integral(const argumentType& rmin,const argumentType& rmax,double length)const{
	returnType s(0);
	for(int i=0;i<point;i++){
		argumentType r=0.5*((xi[i]+1.0)*rmax+(1.0-xi[i])*rmin);
		s+=integrand(r)*w[i];
	}
	return 0.5*length*s;
}

#endif // _LEGENDRE_GAUSS_FORMULA_1D_H_
