#include <cblas.h>
#include <iostream>
#include <math.h>
#include <omp.h>

#define N 1024

void matmul(double* a, double* b, double* c)
{
    #pragma omp parallel for collapse(2)
    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < N; i++)
        {
            double r = 0;
            for(int k = 0; k < N; k++)
            {
                r += a[j * N + k] * b[k * N + i];
            }
            c[j * N + i] = r;
        }
    }
}

double err(double* c1, double* c2)
{  
    double err = 0;
    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < N; i++)
        {
            err += std::pow(c1[j * N + i] - c2[j * N + i], 2);
        }
    }
    return std::sqrt(err / N / N);
}

int main()
{
	std::cout<<"rtfmm_3d_test_surface"<<std::endl;
    std::cout<<omp_get_max_threads()<<std::endl;

	double* A = new double[N * N];
    double* B = new double[N * N];
    double* C = new double[N * N];
    double* C_ = new double[N * N];

	for (int i = 0; i < N * N; i++) 
    {
		A[i] = (double)rand() / RAND_MAX;
		B[i] = (double)rand() / RAND_MAX;
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, A, N, B, N, 0.0, C, N);

    matmul(A,B,C_);

    std::cout<<"error = "<<err(C, C_)<<std::endl;

	delete[] A;
    delete[] B;
    delete[] C;

	return 0;
}