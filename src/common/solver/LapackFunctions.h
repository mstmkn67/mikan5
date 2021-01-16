#ifndef _LAPACK_FUNCTIONS_H_
#define _LAPACK_FUNCTIONS_H_

#include <iostream>
using namespace std;

////////////original subroutines in LAPACK////////////////////////
extern"C"{
	//int ilaenv_(int* ISPEC,char* NAME,char* OPTS,int* N1,int* N2,int* N3,int* N4);
	
	//void dgbsv_(int* N,int* KL,int* KU,int* NRHS,double* AB,int* LDAB,int* IPIV,double* B,int* INFO);
	void dgesv_(int* N,int* NRHS,double* A,int* LDA,int* IPIV,double* B,int* LDB,int* INFO);

	void dgetrf_(int* M,int* N,double* A,int* LDA,int* IPIV,int* INFO);
	void dgetri_(int* N,double* A,int* LDA,int* IPIV,double* WORK,int* LWORK,int* INFO);

	void dpptrf_(char* UPLO,int* N,double* AP,int* INFO);
	void dpptri_(char* UPLO,int* N,double* AP,int* INFO);

	void dsyev_(char* JOBS,char* UPLO,int* N,double* A,int* LDA,double* W,double* WORK,int* LWORK,int* INFO);
}

namespace lap{
//////////////////////////////////////////////////////////////////
///////////interface for c++ ////////////////////////////////
//void solve_linear_equations_band(int* n,int* kl,int* ku,);
int solve_linear_equations_general(int n,double* A,double* B);
int solve_linear_equations_general(int n,double* A,double* B,double* X);//solutions are not written in B

int get_inverse_matrix_symmetric(int n,double* A);//Cholesky decomposition
int get_inverse_matrix_general(int n,double* A);//LU decomposition

int get_LU_decomposition(int n,double* A);
int get_Cholesky_decomposition(int n,double* A);

int get_eigenvalues_symmetric(int n,double* A,double* W);
int get_eigenvalues_and_vectors_symmetric(int n,double* A,double* W);


//general matrix index
// i,j (=0...n-1) index of n x n matrix 
inline int index(int i,int j,int n){
	return n*i+j;
}
//transpose version to use lapack
inline int index_(int i,int j,int n){
	return n*j+i;
}
/*
inline int left(int index){
	return index-1;
}
inline int right(int index){
	return index+1;
}
inline int up(int index,int n){
	return index-n;
}
inline int down(int index,int n){
	return index+n;
}
*/
//a,b (=0,1,2) index of  ei,ej sub matrix of 3n x 3n matrix

inline int mindex(int ei,int ej,int a,int b,int n){
	return 3*n*(3*ei+a)+3*ej+b;
}
//transpose verstion to use lapack
inline int mindex_(int ei,int ej,int a,int b,int n){
	return 3*n*(3*ej+b)+3*ei+a;
}
/*
inline int mleft(int index){
	return index-1;
}
inline int mright(int index){
	return index+1;
}
inline int mup(int index,int n){
	return index-3*n;
}
inline int mdown(int index,int n){
	return index+3*n;
}
*/


//symmetric compressed matrix index
// lower triangle matrix in C or C++ (upper one in Fortran) when you use symmetric case
// i,j (=0...n-1) index of n x n matrix 

inline int sindex(int i,int j){ 
	return i>=j ? (i*(i+1)/2+j):(j*(j+1)/2+i);
}

/*
inline int sleft(int index){
	return index-1;
}
inline int sright(int index){
	return index+1;
}
inline int sup(int index,int i,int n){
	return index-n+i;
}
inline int sdown(int index,int i,int n){
	return index+n-i-1;
}
*/
//a,b (=0,1,2) index of  ei,ej sub matrix of 3n x 3n matrix
inline int msindex(int ei,int ej,int a,int b){
	int i=3*ei+a,j=3*ej+b;
	return i>=j ? (i*(i+1)/2+j):(j*(j+1)/2+i);
}
/*
inline int msleft(int index){
	return index-1;
}
inline int msright(int index){
	return index+1;
}
inline int msup(int index,int i,int a,int n){
	return index-3*n+3*i+a;
}
inline int msdown(int index,int i,int a,int n){
	return index+3*n-3*i-a-1;
}
*/
}
/////////////////////////////////////////

#endif // _LAPACK_FUNCTIONS_H_
