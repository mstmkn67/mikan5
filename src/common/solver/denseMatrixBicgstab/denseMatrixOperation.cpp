#include<stdio.h>
#include "denseMatrixOperation.h"
//
int denseMatVec(
	   double alpha, 
	   double *x, 
	   double beta, 
	   double *y, 
	   double *A_val, 
	   //int *col_index, 
	   //int *row_ptr, 
	   int N 
	   )
{
  int i,j;
  double sum;
  for(i = 0; i < N; i++){
    y[i] = beta * y[i];
    sum = 0.0;
    //for(j = i*N; j <= (i+1)*N-1; j++){
    for(j = 0; j <N; j++){
      //sum += A_val[j] * x[col_index[j]];
      sum += A_val[i*N+j] * x[j];
    }
    y[i] += alpha * sum;
  }
}
//--------------------------------------------------------------------




