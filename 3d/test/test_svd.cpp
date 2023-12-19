#include "type.h"
#include "body.h"
#include "surface.h"
#include <cblas.h>
#include <math.h>
#include "stdlib.h"

#define M 5
#define N 6

void print_matrix(std::string desc, int m, int n, double* a) 
{
    int i, j;
    printf( "\n%s : \n", desc.c_str());
    for( j = 0; j < m; j++ ) 
    {
        for( i = 0; i < n; i++ ) 
        {
            printf("%6.4f ", a[j*n+i]);
        }
        printf( "\n" );
    }
}

extern "C" {
    void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u,
               int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);
}

double err(double* c1, double* c2, int m, int n)
{  
    double err = 0;
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            err += std::pow(c1[j * n + i] - c2[j * n + i], 2);
        }
    }
    return std::sqrt(err / M / N);
}

void check_result(double* u, double* s, double* vt, double* a)
{
    double* S = new double[M * N];
    for(int j = 0; j < M; j++)
    {
        for(int i = 0; i < N; i++)
        {
            if(i == j)
            {
                S[j * N + i] = s[j];
            }
            else
            {
                S[j * N + i] = 0;
            }
        }
    }
    print_matrix("S", M, N, S);

    double* SVT = new double[M * N];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, N, 1.0, S, N, vt, N, 0.0, SVT, N);
    print_matrix("SVT", M, N, SVT);

    double* USVT = new double[M * N];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, M, 1.0, u, M, SVT, N, 0.0, USVT, N);
    print_matrix("USVT", M, N, USVT);

    std::cout<<"\nerror = "<<err(a,USVT,M,N)<<std::endl;

    delete[] S;
    delete[] SVT;
    delete[] USVT;
}

int main(int argc, char* argv[])
{
    std::cout<<"rtfmm_3d_test_svd"<<std::endl;

    int m = M, n = N, info, lwork;
    double wkopt;
    double* work;
    /* Local arrays */
    double a[M*N] = {
        8.79, 6.11, -9.15, 9.57, -3.49, 9.84,
        9.93, 6.91, -7.93, 1.64, 4.02, 0.15,
        9.83, 5.04, 4.86, 8.83, 9.80, -8.99,
        5.45, -0.27, 4.85, 0.74, 10.00, -6.02,
        3.16, 7.98, 3.01, 5.80, 4.27, -5.31
    };
    double* A = new double[M * N];
    for(int j = 0; j < M; j++)
    {
        for(int i = 0; i < N; i++)
        {
            A[j * N + i] = a[j * N + i];
        }
    }
    int k = std::min(m,n);
    char jobu = 'S';
    char jobvt = 'S';
    double u[M*M], vt[N*N];
    double* s = new double[k];
    lwork = -1;
    dgesvd_(&jobu, &jobvt, &n, &m, a, &n, s, vt, &n, u, &k, &wkopt, &lwork,&info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );

    /* Compute SVD */
    dgesvd_(&jobu, &jobvt, &n, &m, a, &n, s, vt, &n, u, &k, work, &lwork,&info );

    /* Check for convergence */
    if(info > 0) 
    {
        printf( "The algorithm computing SVD failed to converge.\n" );
        exit( 1 );
    }

    print_matrix( "S", 1, k, s);
    print_matrix( "U", m, m, u);
    print_matrix( "V^t", n, n, vt);
    print_matrix( "a", m, n, a); // initial a was destory
    check_result(u,s,vt,A);

    free((void*)work);

    return 0;
}

