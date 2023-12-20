#include "mathfunc.h"

void rtfmm::mat_vec(int m, int n, std::vector<real>& A, std::vector<real>& b, std::vector<real>& c, MathType type)
{
    if(type == MathType::naive)
    {
        for(int j = 0; j < m; j++)
        {
            real res = 0.0;
            for(int i = 0; i < n; i++)
            {
                res += A[j * n + i] * b[i];
            }
            c[j] = res;
        }
    }
    else if(type == MathType::blas)
    {
        //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, 1, n, 1.0, A.data(), n, b.data(), 1, 0.0, c.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A.data(), n, b.data(), 1, 0.0, c.data(), 1);
    }
    else
    {

    }
}