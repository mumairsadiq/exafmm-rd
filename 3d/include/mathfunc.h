#pragma once
#include "type.h"
#include <math.h>
#include <cblas.h>
#include <vector>

extern "C" {
    void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u,
               int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);
}

namespace rtfmm
{

enum class MathType
{
    naive,
    blas,
    cublas
};

/**
 * @brief matrix X vec.
 * 
 * @param m rows
 * @param n cols
 * @param A m by n matrix
 * @param b n by 1 vector
 * @param c m by 1 vector
 * @param type underlaying math function type
 */
void mat_vec(int m, int n, matrix A, matrix b, matrix& c, MathType type); 

/**
 * @brief matrix X matrix. 
 * 
 * @param m rows
 * @param n cols of A / rows of B
 * @param k cols of B
 * @param A m by n matrix
 * @param B n by k matrix
 * @param C m by k matrix
 * @param type underlaying math function type
 */
void mat_mat(int m, int n, int k, matrix A, matrix B, matrix& C, MathType type); 

void svd(int m, int n, matrix A, matrix& U, matrix& S, matrix& VT);

matrix transpose(int m, int n, matrix& A);

void print_matrix(int m, int n, matrix& A);

}