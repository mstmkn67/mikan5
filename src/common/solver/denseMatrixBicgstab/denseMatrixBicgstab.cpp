#include "denseMatrixBicgstab.h"

/*****************************************************************************/
//
//      work[7][N] : work space 
//
//      work[0] : r
//      work[1] : rtld
//      work[2] : p
//      work[3] : phat
//      work[4] : v
//      work[5] : shat
//      work[6] : t
//      work[0] : s
//
/*****************************************************************************/
int dense_matrix_bicgstab(
						 int N, 
						 double *x, 
						 double *b, 
						 double** work,
						 double *val_A, 
						 //int *col_index, 
						 //int *row_ptr,
						 double *val_ILU, 
						 double *pivots, 
						 //int *diag_ptr,
						 int ITRMAX, 
						 double eps
						 )
{
  int n,i,j;
  double omega;
  double alpha,beta;
  double rho1,rho2;
  double bnorm,snorm,xnorm,rnorm;
  double temp1,temp2,temp3;
  //double (*work)[N] = (double (*)[N])(_work);
  /**** alias workspace columns ****/
  int r,rtld,p,phat,v,s,shat,t;
  r    = 0;
  rtld = 1;
  p    = 2;
  phat = 3;
  v    = 4;
  shat = 5;
  t    = 6;
  s    = 0;

  //<< computing norm of vector b >>
  bnorm = GetNorm(b,N);

	if( bnorm == 0.0e0 ){       
		for(i = 0; i < N; i++){ x[i] = 0.0;	}
    //cout << "iter=1 " << endl;
    //cout << "bnorm= " << bnorm << endl;
		return(true);
	}

  //<< initialize >>
  denseMatVec( -1.0, x, 1.0, b, val_A, /*col_index, row_ptr,*/ N ); // b = b - A x 
  memcpy( work[r], b, N*sizeof(double) );                  // r = b
  memcpy( work[rtld], work[r], N*sizeof(double) );         // r~= r

  //-----<< iteration >>-----//

  for(n = 1; n <= ITRMAX; n++){

    rho1 = InnerProduct(work[r], work[rtld], N);      // rho1 = r . r~

    if(n == 1){ 
      memcpy( work[p], work[r], sizeof(double)*N ); // p = r
    }else{
      beta = ( rho1 / rho2 )*( alpha / omega ); 
      for( i=0; i < N; i++ ){
				work[p][i] = work[r][i] + beta*(work[p][i]-omega * work[v][i]);
				// p = r + beta * ( p - omega * v )
      }
    }
    
    //--- _solve phat ---//   
    // ILU_A p^ = p
    denseMatrixLUDecomposedMatrixSolver(
														 val_ILU, 
														 //col_index, 
														 //row_ptr, 
														 pivots, 
														 //diag_ptr, 
														 work[p], 
														 work[phat], 
														 N 
														 );
    
    /**** computing v ****/
    // v = A p^ 
    denseMatVec(ONE, work[phat], ZERO, work[v], val_A, /*col_index, row_ptr,*/ N);
    
    // alpha = rho1 / r~ . v 
    alpha = rho1 / InnerProduct(work[rtld], work[v], N);
    
		// s = r - alpha * v
    for(i =0; i < N; i++){
      work[s][i] = work[r][i] - alpha * work[v][i];
    }
    
    // check norm of s 
    // snorm = || s ||
    snorm = GetNorm( work[s], N );
    
    if(snorm <= eps){
      // x = x + alpha * p^
      for(i = 0; i < N; i++){
				x[i] += alpha * work[phat][i];
      }
      //printf("iter=%d (snorm/bnorm)^2=%f\t", n,(snorm/bnorm)*(snorm/bnorm));
			//cout << "iter=" << n << " snorm/bnorm=" << snorm/bnorm << "\t" ;
			//cout << "bnorm=" << bnorm << " snorm=" << snorm << endl;
      return(true);
    }
    
    // Solve Linear Equation : (ILU) s^ = s 
    denseMatrixLUDecomposedMatrixSolver(
														 val_ILU, /*col_index, row_ptr,*/ pivots, /*diag_ptr,*/ 
														 work[s], work[shat], N );
    
    // t = A s^
    denseMatVec(ONE, work[shat], ZERO, work[t], val_A, /*col_index, row_ptr,*/ N);
    
    // omega = (t.s)/(t.t)
    omega=InnerProduct(work[t],work[s],N)/InnerProduct(work[t],work[t],N);
    
    for(i = 0; i < N; i++){

      // x = x + alpha * p~ + omega * s^ 
      x[i]       = x[i]       + alpha * work[phat][i] + omega * work[shat][i];

      // r = s - omega * t
      work[r][i] = work[s][i] - omega * work[t][i];

    }
    //-------------------------------------------//
    // check convergence //
    xnorm = GetNorm(x,N);
    rnorm = GetNorm(work[r],N);
    //cout << "iter="  << n << " rnorm/bnorm=" << rnorm/bnorm << "\t";
    //cout << "bnorm=" << bnorm << " rnorm="   << rnorm << endl;
    if(rnorm <= eps * (xnorm + bnorm)){ return(true); }
    rho2 = rho1;
  }
  /************************* iteration  end *************************/
  cerr << "bicgsta failed" << endl;
	cerr << "iteration number exceeds Maximum Iteration:" << ITRMAX << endl;
  exit(1);
  return(false);
}
//============================================================================
//----------------------------------------------------------------------------
//
//  When a matrix A is a LU Decomposed Matrix 
//  
//                A y = x
//  
//  will be solved using a following function.
//
//----------------------------------------------------------------------------
int denseMatrixLUDecomposedMatrixSolver(
														 double      *val_LU,    // Input
														 //int         *col_index, // Input
														 //int         *row_ptr,   // Input
														 double      *pivots,    // Input
														 //int         *diag_ptr,  // Input
														 double      *b,         // Input
														 double      *x,         // Output
														 const int&  N           // Input 
														 )
{
  int i,j;  double sum;
  for(i = 0; i < N; i++){
    sum = 0.0;
    //for(j = row_ptr[i]; j < diag_ptr[i]; j++){
    for(j = 0; j < i; j++){
      //sum += val_LU[j] * x[col_index[j]];
      sum += val_LU[i*N+j] * x[j];
    }
    x[i] = pivots[i] * (b[i] - sum);
  }
  for(i = N-1; i >= 0; i--){
    sum = 0.0;
    //for(j = diag_ptr[i]+1; j < row_ptr[i+1]; j++){
    for(j = i+1; j < N; j++){
      //sum  += val_LU[j] * x[col_index[j]];
      sum  += val_LU[i*N+j] * x[j];
      x[i] -= - pivots[i] * sum;
    }
  }
}
//----------------------------------------------------------------------------

