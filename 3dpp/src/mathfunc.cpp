#include "mathfunc.h"

rtfmm::Matrix rtfmm::mat_vec_mul(Matrix A, Matrix b)
{
    int m = A.m, n = A.n;
    Matrix C(m, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A.d.data(), n, b.d.data(), 1, 0.0, C.d.data(), 1);
    return C;
}

rtfmm::Matrix rtfmm::mat_mat_add(Matrix A, Matrix B)
{
    assert_exit(A.m == B.m && A.n == B.n, "mat_mat_add inconsitent size");
    int m = A.m, n = A.n;
    Matrix C(m,n);
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            C[j * n + i] = A[j * n + i] + B[j * n + i];
        }
    }
    return C;
}

rtfmm::Matrix rtfmm::mat_mat_mul(Matrix A, Matrix B)
{
    int m = A.m, n = A.n, k = B.n;
    Matrix C(m,k);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, 1.0, A.d.data(), n, B.d.data(), k, 0.0, C.d.data(), k);
    return C;
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
    VT = vt;
    S = Matrix(m,n);

    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            if(i == j && j < k) S[j * n + i] = s[j];
            else S[j * n + i] = 0;
        }
    }
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
    Matrix U, S, VT;
    svd(A, U, S, VT);
    Matrix UT = transpose(U);
    Matrix V = transpose(VT);
    Matrix Sinv = pseudo_inverse(S);
    Matrix UTb = mat_vec_mul(UT, b);
    Matrix SinvUTb = mat_vec_mul(Sinv, UTb);
    Matrix x = mat_vec_mul(V, SinvUTb);
    return x;
}

rtfmm::Matrix rtfmm::pseudo_inverse(Matrix S)
{
    int m = S.m, n = S.n;
    int k = std::min(m, n);
    real s_max = 0;
    const real EPS = 4 * 1e-16;
    for(int i = 0; i < k; i++)
    {
        s_max = std::max(s_max, S[i * n + i]);
    }
    for(int i = 0; i < k; i++)
    {
        S[i * n + i] = S[i * n + i] > s_max * EPS ? 1.0 / S[i * n + i] : 0.0;
    }
    return transpose(S);
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

void rtfmm::print_matriv(Matriv& A)
{
    printf("{\n");
    for(int j = 0; j < A.m; j++)
    {
        printf("{");
        for(int i = 0; i < A.n; i++)
        {
            printf("(%.4f,%.4f,%.4f),", A[j * A.n + i][0], A[j * A.n + i][1], A[j * A.n + i][2]);
        }
        printf("}\n");
    }
    printf("}\n");
}