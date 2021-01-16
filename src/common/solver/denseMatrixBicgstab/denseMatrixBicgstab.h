//=======================================================================
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;
//=======================================================================
#include "MatrixOperation.h"
#include "denseMatrixOperation.h"
#include "denseMatrixILUdcmp.h"
//#include "LinearEquationSolver.h"
//#include "Ordering.h"
//=======================================================================
#ifndef __Dense_matrix_Biconjugate_Gradient_with_stabilizing_methd__
#define __Dense_matrix_Biconjugate_Gradient_with_stabilizing_methd__

#define ONE  1
#define ZERO 0

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
						 );

int denseMatrixLUDecomposedMatrixSolver(
														 double      *val_LU,    // Input
														 //int         *col_index, // Input
														 //int         *row_ptr,   // Input
														 double      *pivots,    // Input
														 //int         *diag_ptr,  // Input
														 double      *b,         // Input
														 double      *x,         // Output
														 const int&  N           // Input 
														 );
#endif
//=======================================================================
