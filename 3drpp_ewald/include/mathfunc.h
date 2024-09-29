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

void mat_scale(Matrix& A, real scale);

/**
 * @brief c = k x A x b
 * 
 * @param A m by n matrix
 * @param b n by 1 vector
 * @param k real number
 * @return m by 1 vector
 */
Matrix mat_vec_mul(const Matrix& A, const Matrix& b, real k = 1.0); 

/**
 * @brief c = k x A x b
 * 
 * @param A m by n matrix
 * @param b n by 1 vector
 * @param c m by 1 vector
 * @param k real number
 */
void mat_vec_mul(const Matrix& A, const Matrix& b, Matrix& c, real k = 1.0); 

/**
 * @brief C = A + B
 * 
 * @param A m by n matrix
 * @param B m by n matrix
 * @return m by n matrix
 */
Matrix mat_mat_add(Matrix& A, Matrix& B); 

/**
 * @brief A += B
 * 
 * @param A m by n matrix
 * @param B m by n matrix
 */
void mat_mat_increment(Matrix& A, Matrix& B); 

/**
 * @brief C = A x B
 * 
 * @param A m by n matrix
 * @param B n by k matrix
 * @return m by k matrix
 */
Matrix mat_mat_mul(const Matrix& A, const Matrix& B); 

Matrix mat_mat_mul_naive(const Matrix& A, const Matrix& B);

/**
 * @brief solve A = U x S x VT
 * 
 * @param A input matrix A (m x n)
 * @param U output matrix U (m x m)
 * @param S output vector S (m x n) 
 * @param VT output matrix VT (n x n)
 */
void svd(Matrix A, Matrix& U, Matrix& S, Matrix& VT);

/**
 * @brief transpose matrix A
 * @param A matrix A
 */
Matrix transpose(Matrix A);

void transpose_inplace(Matrix& A);

/**
 * @brief solve Ax = b using svd
 * 
 * @param A matrix A (m x n)
 * @param b vector b (m x 1)
 * @return x (n x 1)
 */
Matrix linear_equation_system_svd(Matrix A, Matrix b);

/**
 * @brief get pseudo inversed matrix of S
 * @param S diagonal matrix S(m x n)
 * @param 
*/
Matrix pseudo_inverse(Matrix S);

void print_matrix(Matrix& mat);

void print_matriv(Matriv& mav);

real matrix_L2(Matrix& A, Matrix& B);

Matrix identity(int m, int n);

vec2r min_max(Matrix& A);

}