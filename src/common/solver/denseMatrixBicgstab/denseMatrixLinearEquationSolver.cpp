#include "denseMatrixLinearEquationSolver.h"
//----------------------------------------------------------------------------
int dense_matrix_solve_linear_equations(int N,double* A,double* B,double* x,double eps,int ITRMAX)
{
	double *pivots=new double[N];
	double *ILU_A=new double[N*N];
	double *work[7];double *_work=new double[7*N];
	for(int i=0;i<7;i++){
		work[i]=_work+i*N;
	}
	for(int i=0;i<N;i++){
		pivots[i]=0.0;x[i]=0.0;
		for(int j=0;j<N;j++){
			ILU_A[i*N+j]=0.0;
		}
	}
	denseMatrixILUdcmp(A,N,pivots,ILU_A);
	int n=dense_matrix_bicgstab(N,x,B,work,A,ILU_A,pivots,ITRMAX,eps);
	delete[] pivots;delete[] ILU_A;delete[] _work;
	return n;
}

DenseMatrixLinearEquationsSolver::DenseMatrixLinearEquationsSolver(int N_,double* A_,double* B_,
double* x_,double eps_,int ITRMAX_):
N(N_),A(A_),B(B_),x(x_),eps(eps_),ITRMAX(ITRMAX_){
	pivots=new double[N];
	ILU_A=new double[N*N];
	_work=new double[7*N];
	for(int i=0;i<7;i++){	work[i]=_work+i*N;}
	for(int i=0;i<N;i++){
		pivots[i]=0.0;x[i]=0.0;
		for(int j=0;j<N;j++){
			ILU_A[i*N+j]=0.0;
		}
	}
}

DenseMatrixLinearEquationsSolver::~DenseMatrixLinearEquationsSolver(){
	delete[] pivots;delete[] ILU_A;delete[] _work;
}

int DenseMatrixLinearEquationsSolver::solve(){
	for(int i=0;i<N;i++){
		pivots[i]=0.0;x[i]=0.0;
		for(int j=0;j<N;j++){
			ILU_A[i*N+j]=0.0;
		}
	}
	denseMatrixILUdcmp(A,N,pivots,ILU_A);
	return dense_matrix_bicgstab(N,x,B,work,A,ILU_A,pivots,ITRMAX,eps);
}

void DenseMatrixLinearEquationsSolver::change_properties(int N_,double* A_,double* B_,double* x_,double eps_,int ITRMAX_){
	N=N_;A=A_;B=B_;x=x_;eps=eps_;ITRMAX=ITRMAX_;
	delete pivots;delete ILU_A;delete _work;
	pivots=new double[N];ILU_A=new double[N*N];_work=new double[7*N];
	for(int i=0;i<7;i++){	work[i]=_work+i*N;}
	for(int i=0;i<N;i++){
		pivots[i]=0.0;
		for(int j=0;j<N;j++){
			ILU_A[i*N+j]=0.0;
		}
	}
}
