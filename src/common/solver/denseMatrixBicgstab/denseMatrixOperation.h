#ifndef __denseMatrixOperation__
#define __denseMatrixOperation__
#include <math.h>
int denseMatVec(double alpha, double *x, double beta, double *y, 
	   double *val, /*int *col_index, int *row_ptr,*/ int N);
#endif
