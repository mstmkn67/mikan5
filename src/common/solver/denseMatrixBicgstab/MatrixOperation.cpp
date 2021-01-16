#include<stdio.h>

//--------------------------------------------------------------------
// ベクトルx１とx２の内積を計算する内部関数 //
// inner function for computing a inner product of vector x1 and x2 //
#include "MatrixOperation.h"

double InnerProduct( double *x1, double *x2, int N )
{
  int i; double sum = 0.0;
  for(i = 0; i < N; i++){ sum += x1[i]*x2[i]; }
  return sum;
}
//--------------------------------------------------------------------
// inner function for computing a norm of vector //
double GetNorm(double *x, int N)
{
  int i; double norm = 0.0;
  for( i = 0; i < N; i++ ){ norm += x[i] * x[i]; }
  return sqrt(norm);
}
//--------------------------------------------------------------------
int ArrayCopy( double* destArray, double* srcArray, const int& N )
{
  int i; for( i=0; i<N; i++ ){ 
    destArray[i]=srcArray[i]; 
    printf("srcArray[%d]=%f\n",i, srcArray[i]);
    //printf("destArray[%d]=%f\n",i, destArray[i]);
  }
  return(0);
}
//--------------------------------------------------------------------
//
//--------------------------------------------------------------------
//
// y = beta * y + alpha * val * x 
//
int MatVec(
	   double alpha, 
	   double *x, 
	   double beta, 
	   double *y, 
	   double *A_val, 
	   int *col_index, 
	   int *row_ptr, 
	   int N 
	   )
{
  int i,j;
  double sum;
  for(i = 0; i < N; i++){
    y[i] = beta * y[i];
    sum = 0.0;
    for(j = row_ptr[i]; j <= row_ptr[i+1]-1; j++){
      sum += A_val[j] * x[col_index[j]];
    }
    y[i] += alpha * sum;
  }
}
//--------------------------------------------------------------------




