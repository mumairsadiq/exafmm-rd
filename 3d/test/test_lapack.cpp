#include <stdio.h>

extern "C" {
  void dgesv_ ( const int& N, const int& NRHS,
		double** A, const int& LDA, int* IPIV,
		double** B, const int& LDB, int& INFO );
};

template <int N> int dgesv( double A[N][N], double x[N], const double b[N] )
{
  int i, j, info;
  static int ipiv[N];
  static double U[N][N];

  for( j=0; j<N; j++ ){
    for( i=0; i<N; i++ ){
      U[i][j] = A[j][i];
    }
    x[j] = b[j];
  }

  dgesv_( N, 1, (double**)U, N, ipiv, (double**)x, N, info );

  return info;
}

int main(void)
{
  const int N = 4;
  int i, j, info;

  double A[N][N];
  double x[N], b[N];

  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){
      if( i==j ){
	A[i][j] = 2.0;
      }else{
	A[i][j] = 1.0;
      }
    }
    b[i] = 5.0;
  }

  info = dgesv( A, x, b );

  printf("# info=%d.\n", info );

  printf("# Solution.\n");
  for( i=0; i<N; i++ ){
    printf("%+f\n", x[i] );
  }

  printf("# Error.\n");
  for( i=0; i<N; i++ ){
    double sum=0.0;
    for( j=0; j<N; j++ ){
      sum += A[i][j]*x[j];
    }
    sum -= b[i];

    printf("%+f\n", sum );
  }

  return 0;
}