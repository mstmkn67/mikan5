#ifndef _DENSE_MATRIX_LINEAR_EQUATION_SOLVER_H_
#define _DENSE_MATRIX_LINEAR_EQUATION_SOLVER_H_

#include "denseMatrixBicgstab.h"
int dense_matrix_solve_linear_equations(int N,double* A,double* B,double* x,double eps=1.0e-5,int ITRMAX = 1000);

class DenseMatrixLinearEquationsSolver{
public:
	DenseMatrixLinearEquationsSolver(int N,double* A,double* B,double* x,double eps=1.0e-5,int ITRMAX=1000);
	virtual ~DenseMatrixLinearEquationsSolver();
	virtual int solve();//BÇÃílÇ™ïœÇÌÇÈÇÃÇ≈íçà”
	virtual void change_properties(int N,double* A,double* B,double* x,double eps=1.0e-5,int ITRMAX=1000);
private:
	int N;
	double *A,*B,*x;
	double *pivots;
	double *ILU_A;
	double *work[7];
	double *_work;
	double eps;
	int ITRMAX;
};

#endif // _DENSE_MATRIX_LINEAR_EQUATION_SOLVER2_H_
