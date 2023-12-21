#include "mathfunc.h"

void rtfmm::mat_vec(int m, int n, matrix A, matrix b, matrix& c, MathType type)
{
    matrix C(m * 1);
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
        //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, 1, n, 1.0, A.data(), n, b.data(), 1, 0.0, c.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A.data(), n, b.data(), 1, 0.0, C.data(), 1);
    }
    else
    {

    }
    c = C;
}

void rtfmm::mat_mat(int m, int n, int k, matrix A, matrix B, matrix& C, MathType type)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, 1.0, A.data(), n, B.data(), k, 0.0, C.data(), k);
}

void rtfmm::svd(int m, int n, matrix A, matrix& U, matrix& S, matrix& VT)
{
    int k = std::min(m,n);
    
    matrix u(m * m);
    matrix vt(n * n);
    matrix s(k);

    char jobu = 'S';
    char jobvt = 'S';
    double wkopt;
    int info;
    int lwork = -1;
    dgesvd_(&jobu, &jobvt, &n, &m, A.data(), &n, s.data(), vt.data(), &n, u.data(), &k, &wkopt, &lwork,&info);
    lwork = (int)wkopt;
    matrix wbuff(lwork);
    dgesvd_(&jobu, &jobvt, &n, &m, A.data(), &n, s.data(), vt.data(), &n, u.data(), &k, wbuff.data(), &lwork,&info);

    U = u;
    S = s;
    VT = vt;
}

rtfmm::matrix rtfmm::transpose(int m, int n, matrix& A)
{
    matrix res(n * m);
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            res[i * m + j] = A[j * n + i];
        }
    }
    return res;
}

void rtfmm::print_matrix(int m, int n, matrix& A)
{
    printf("{\n");
    for(int j = 0; j < m; j++)
    {
        printf("{");
        for(int i = 0; i < n; i++)
        {
            printf("%.4f,", A[j * n + i]);
        }
        printf("}\n");
    }
    printf("}\n");
}