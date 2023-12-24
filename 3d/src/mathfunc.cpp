#include "mathfunc.h"

rtfmm::Matrix rtfmm::mat_vec(Matrix A, Matrix b, MathType type)
{
    int m = A.m, n = A.n;
    Matrix C(m, 1);
    if(type == MathType::naive)
    {
        for(int j = 0; j < m; j++)
        {
            real res = 0.0;
            for(int i = 0; i < n; i++)
            {
                res += A[j * n + i] * b[i];
            }
            C[j] = res;
        }
    }
    else if(type == MathType::blas)
    {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A.d.data(), n, b.d.data(), 1, 0.0, C.d.data(), 1);
    }
    else
    {

    }
    return C;
}

void rtfmm::mat_mat(Matrix A, Matrix B, Matrix& C, MathType type)
{
    int m = A.m, n = A.n, k = B.n;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, 1.0, A.d.data(), n, B.d.data(), k, 0.0, C.d.data(), k);
}

void rtfmm::svd(Matrix A, Matrix& U, Matrix& S, Matrix& VT)
{
    int m = A.m, n = A.n;
    int k = std::min(m,n);
    
    Matrix u(m,m);
    Matrix vt(n,n);
    Matrix s(k,1);

    char jobu = 'S';
    char jobvt = 'S';
    double wkopt;
    int info;
    int lwork = -1;
    dgesvd_(&jobu, &jobvt, &n, &m, A.d.data(), &n, s.d.data(), vt.d.data(), &n, u.d.data(), &k, &wkopt, &lwork,&info);
    lwork = (int)wkopt;
    Matrix wbuff(lwork, 1);
    dgesvd_(&jobu, &jobvt, &n, &m, A.d.data(), &n, s.d.data(), vt.d.data(), &n, u.d.data(), &k, wbuff.d.data(), &lwork,&info);

    U = u;
    S = s;
    VT = vt;
}

rtfmm::Matrix rtfmm::transpose(Matrix A)
{
    int m = A.m;
    int n = A.n;
    Matrix res(n,m);
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            res[i * m + j] = A[j * n + i];
        }
    }
    return res;
}

rtfmm::Matrix rtfmm::linear_equation_system_svd(Matrix A, Matrix b)
{
    int m = A.m, n = A.n;
    int k = std::min(m, n);
    Matrix U(m,m), S(k,1), VT(n,n);
    svd(A, U, S, VT);
    Matrix UT = transpose(U);
    Matrix V = transpose(VT);
    Matrix Sinv(m,n);
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            if(i == j && j < k) Sinv[j * n + i] = 1.0 / S[j];
            else Sinv[j * n + i] = 0;
        }
    }
    Sinv = transpose(Sinv);
    /* x = V * S^-1 * U^T * b */
    Matrix UTp = mat_vec(UT, b, MathType::blas);
    Matrix SinvUTp = mat_vec(Sinv, UTp, MathType::blas);
    Matrix x = mat_vec(V, SinvUTp, MathType::blas);
    return x;
}

void rtfmm::print_matrix(Matrix& A)
{
    printf("{\n");
    for(int j = 0; j < A.m; j++)
    {
        printf("{");
        for(int i = 0; i < A.n; i++)
        {
            printf("%.4f,", A[j * A.n + i]);
        }
        printf("}\n");
    }
    printf("}\n");
}