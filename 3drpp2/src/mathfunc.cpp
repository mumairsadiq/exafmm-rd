#include "mathfunc.h"

void rtfmm::mat_scale(Matrix& A, real scale)
{
    for(int i = 0; i < A.m * A.n; i++)
    {
        A.d[i] *= scale;
    }
}

rtfmm::Matrix rtfmm::mat_vec_mul(const Matrix& A, const Matrix& b, real k)
{
    int m = A.m, n = A.n;
    Matrix C(m, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, k, A.d.data(), n, b.d.data(), 1, 0.0, C.d.data(), 1);
    return C;
}

void rtfmm::mat_vec_mul(const Matrix& A, const Matrix& b, Matrix& c, real k)
{
    int m = A.m, n = A.n;
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, k, A.d.data(), n, b.d.data(), 1, 0.0, c.d.data(), 1);
}

rtfmm::Matrix rtfmm::mat_mat_add(Matrix& A, Matrix& B)
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

void rtfmm::mat_mat_increment(Matrix& A, Matrix& B)
{
    assert_exit(A.m == B.m && A.n == B.n, "mat_mat_add inconsitent size");
    int m = A.m, n = A.n;
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            A[j * n + i] += B[j * n + i];
        }
    }
}

rtfmm::Matrix rtfmm::mat_mat_mul(const Matrix& A, const Matrix& B)
{
    assert_exit(A.n == B.m, "mat_mat_mul size error");
    int m = A.m, n = A.n, k = B.n;
    Matrix C(m,k);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, 1.0, A.d.data(), n, B.d.data(), k, 0.0, C.d.data(), k);
    return C;
}

rtfmm::Matrix rtfmm::mat_mat_mul_naive(const Matrix& A, const Matrix& B)
{
    int m = A.m, n = A.n, k = B.n;
    Matrix C(m,k);
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < k; i++)
        {
            real res = 0.0;
            for(int z = 0; z < n; z++)
            {
                res += A.d[j * n + z] * B.d[z * k + i];
            }
            C[j * k + i] = res;
        }
    }
    return C;
}

void rtfmm::svd(Matrix A, Matrix& U, Matrix& S, Matrix& VT)
{
    int m = A.m, n = A.n;
    int k = std::min(m,n);

    U = Matrix(m,m);
    VT = Matrix(n,n);
    Matrix s(k,1);

    char jobu = 'S';
    char jobvt = 'S';
    double wkopt;
    int info;
    int lwork = -1;
    dgesvd_(&jobu, &jobvt, &n, &m, A.d.data(), &n, s.d.data(), VT.d.data(), &n, U.d.data(), &k, &wkopt, &lwork,&info);
    lwork = (int)wkopt;
    Matrix wbuff(lwork, 1);
    dgesvd_(&jobu, &jobvt, &n, &m, A.d.data(), &n, s.d.data(), VT.d.data(), &n, U.d.data(), &k, wbuff.d.data(), &lwork,&info);

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

void rtfmm::transpose_inplace(Matrix& A)
{
    assert_exit(A.m == A.n, "inplace transpose matrix not a square");
    int m = A.m;
    int n = A.n;
    for(int j = 0; j < m; j++)
    {
        for(int i = j + 1; i < n; i++)
        {
            std::swap(A[j * n + i],A[i * n + j]);
        }
    }
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
        s_max = std::max(s_max, std::abs(S[i * n + i]));
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

rtfmm::real rtfmm::matrix_L2(Matrix& A, Matrix& B)
{
    assert_exit(A.m == B.m && A.n == B.n, "matrix_L2 size inconsisent");
    real res = 0;
    for(int j = 0; j < A.m; j++)
    {
        for(int i = 0; i < A.n; i++)
        {
            res += std::pow(A.d[j * A.n + i] - B.d[j * A.n + i],2);
        }
    }
    return std::sqrt(res / A.n / A.m);
}

rtfmm::Matrix rtfmm::identity(int m, int n)
{
    Matrix res(m,n);
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            if(i == j)
            {
                res[j * n + i] = 1;
            }
            else
            {
                res[j * n + i] = 0;
            }
        }
    }
    return res;
}

rtfmm::vec2r rtfmm::min_max(Matrix& A)
{
    int m = A.m;
    int n = A.n;
    real minv = A.d[0];
    real maxv = A.d[0];
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            real v = A[j * n + i];
            minv = std::min(minv, v);
            maxv = std::max(maxv, v);
        }
    }
    return {minv, maxv};
}