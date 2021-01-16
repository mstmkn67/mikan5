#include "denseMatrixBicgstab.h"
#include <iostream>
using namespace std;

void set_matrix(int N,double* A,double* B);

int main(){
	int N=3;
	double* A=new double[N*N];
	double* B=new double[N];
	double* x=new double[N];
	
	double *pivots=new double[N];
	double *ILU_A=new double[N*N];
	double *work[7];double *_work=new double[7*N];
	for(int i=0;i<7;i++){
		work[i]=_work+i*N;
	}
	for(int i=0;i<N;i++){
		pivots[i]=0.0;
		for(int j=0;j<N;j++){
			ILU_A[i*N+j]=0.0;
		}
	}
	
	set_matrix(N,A,B);
	denseMatrixILUdcmp(A,N,pivots,ILU_A);
	int n=dense_matrix_bicgstab(N,x,B,work,A,ILU_A,pivots,1000,1e-5);
	
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			cout << A[i*N+j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	for(int i=0;i<N;i++){
		cout << B[i] << " ";
	}
	cout << endl;
	for(int i=0;i<N;i++){
		cout << x[i] << " ";
	}
	cout << endl;
	
	delete[] pivots,ILU_A,_work;
	delete[] A,B,x;
	return 0;
}

void set_matrix(int N,double* A,double* B){
	A[0*N+0]=2.1;A[0*N+1]=0.1;A[0*N+2]=0.2;
	A[1*N+0]=0.3;A[1*N+1]=1.1;A[1*N+2]=0.1;
	A[2*N+0]=0.4;A[2*N+1]=0.1;A[2*N+2]=2.2;
	B[0]=1.0;B[1]=2.0;B[2]=2.0;
}

