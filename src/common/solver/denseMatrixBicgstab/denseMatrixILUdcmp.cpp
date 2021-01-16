//---------------------------------------------------
//     << ILUdcmp >> Incomplete LU decomposition 
//---------------------------------------------------
#include "denseMatrixILUdcmp.h"
int denseMatrixILUdcmp(
	     double     *val_A, 
	     //int        *col_index, 
	     //int        *row_ptr, 
	     //int        *diag_ptr,
	     const int& N,
	     double     *pivots, 
	     double     *val_ILU
	     )
{
  int i,j,k,found;
  double element;
  //for( i = 0; i < N; i++ ){ pivots[i] = val_A[diag_ptr[i]]; }
  for( i = 0; i < N; i++ ){ pivots[i] = val_A[i*(N+1)]; }
	for( i = 0; i < N; i++ ){
		pivots[i] = 1.0/pivots[i];
		//for( j = diag_ptr[i]+1; j < row_ptr[i+1]; j++ ){
		for( j = i+1; j < N; j++ ){
			found = FALSE;
			//for( k = row_ptr[col_index[j]]; k < diag_ptr[col_index[j]]; k++){
			for( k = 0; k < j; k++){
				//if( col_index[k] == i ){
				if( k == i ){
					found = TURE;
					//element = val_A[k];
					element = val_A[j*N+k];
					break;
				}
			}
			if(found == TURE){
				//val_ILU[ diag_ptr[col_index[j]] ] = 
				//val_A[ diag_ptr[col_index[j]] ] - element * pivots[i] * val_A[j] ;
				val_ILU[j*(N+1)] = val_A[j*(N+1)] - element * pivots[i] * val_A[i*N+j] ;
			}
		}
	}
}
//---------------------------------------------------

