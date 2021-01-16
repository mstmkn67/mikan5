#ifndef __dense_Matrix_Imcoplete_LU_decomposition__
#define __dense_Matrix_Imcoplete_LU_decomposition__

//#include "Array.h"
#include <string.h>

#define FALSE 0
#define TURE  1
//---------------------------------------------------
//     << ILUdcmp >> Incomplete LU decomposition 
//---------------------------------------------------
// << Data Structure >>
//  Here, we use the following Data Structure, 
//  CRS (Compressed Row Strage, CRS):
//
//    col:   0 1 2 3 
//         ----------- 
//        A=|1 0 4 0| : 0
//          |3 1 0 0| : 1
//          |0 7 5 9| : 2
//          |4 0 0 6| : 3
//                     row
//
//   The matrix A is described by 
//  the following three arrays : 
//
//      value_A[]=[1,4,3,1,7,5,9,4,6]
//                |   |   |     |   |
//                0   2   4     7   9 = sizeofA + 1
//
//        col_idx[]=[0,2,0,1,1,2,3,0,3]
//        row_ptr[]=[0,2,4,7,9]
//        dia_ptr[]=[0,3,5,8]
//
//  Here, #_Of_Component(row_ptr)=#Of_Row_Of(A) + 1,
//        value_A [ dia_prt[i] ] = A_ii
//
//---------------------------------------------------
// User have to input :
// "val_A", "col_index", "N" = DimOfA (Info:Matrix A)
// then, will get the following values :
// "diag_ptr", "pivots", "val_ILU"
// Using the function << ILUdcmp >>
//---------------------------------------------------
int denseMatrixILUdcmp(
	     double       *val_A,      // Input 
	     //int          *col_index,  // Input 
	     //int          *row_ptr,    // Input
	     //int          *diag_ptr,   // Input
	     const int&   N,           // Input
	     double       *pivots,     // Output
	     double       *val_ILU     // Output
	     );
#endif
