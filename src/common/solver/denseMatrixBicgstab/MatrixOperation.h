#ifndef __MatrixOperation__
#define __MatrixOperation__
#include <math.h>
double InnerProduct( double *x1, double *x2, int N );
double GetNorm( double *x, int N );
int ArrayCopy( double* destArray, double* srcArray, const int& N );
int MatVec(double alpha, double *x, double beta, double *y, 
	   double *val, int *col_index, int *row_ptr, int N);
#endif
