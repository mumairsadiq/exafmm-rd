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
 * @brief c = A x b
 * 
 * @param A m by n matrix
 * @param b n by 1 vector
 * @param type underlaying math function type
 * @return m by 1 vector
 */
Matrix mat_vec(Matrix A, Matrix b, MathType type = MathType::blas); 

/**
 * @brief C = A x B
 * 
 * @param A m by n matrix
 * @param B n by k matrix
 * @param C m by k matrix
 * @param type underlaying math function type
 */
void mat_mat(Matrix A, Matrix B, Matrix& C, MathType type = MathType::blas); 

/**
 * @brief solve A = U x S x VT
 * 
 * @param A input matrix A (m x n)
 * @param U output matrix U (m x m)
 * @param S output vector S (k x 1, k = min(m,n)) 
 * @param VT output matrix VT (n x n)
 */
void svd(Matrix A, Matrix& U, Matrix& S, Matrix& VT);

/**
 * @brief transpose matrix A
 * @param A matrix A
 */
Matrix transpose(Matrix A);

/**
 * @brief solve Ax = b using svd
 * 
 * @param A matrix A (m x n)
 * @param b vector b (m x 1)
 * @return x (n x 1)
 */
Matrix linear_equation_system_svd(Matrix A, Matrix b);

void print_matrix(Matrix& mat);

}