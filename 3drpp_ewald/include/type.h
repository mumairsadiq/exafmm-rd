#pragma once
#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <chrono>
#include <assert.h>
#include "emmintrin.h"
#include "immintrin.h"

#define TIME_BEGIN(a) auto time_begin_##a = std::chrono::high_resolution_clock::now()

#define TIME_END(a)   auto time_end_##a = std::chrono::high_resolution_clock::now();\
					  auto elapse_##a = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_##a - time_begin_##a);\
                      RTLOG("[%s time measured : %.5f seconds.]\n", #a, elapse_##a.count() * 1e-9)

#define tbegin(a) TIME_BEGIN(a)
#define tend(a) TIME_END(a)

#define RTLOG(...) fprintf (stderr, __VA_ARGS__)

namespace rtfmm
{

template <typename ... Args>
std::string format(const std::string& fmt, Args ... args )
{
    size_t len = std::snprintf( nullptr, 0, fmt.c_str(), args ... );
    std::vector<char> buf(len + 1);
    std::snprintf(&buf[0], len + 1, fmt.c_str(), args ... );
    return std::string(&buf[0], &buf[0] + len);
}

inline void title(std::string s)
{
    std::cout<<"<"<<s<<">\n\n";
}

inline void title(char* s)
{
    printf("%s\n", s);
}

using Indices = std::vector<int>;

using real = double;

struct complexr
{
    real r;
    real i;
};

/**
 * @brief alert error message and exit if conditon is false
*/
inline void assert_exit(bool condition, std::string err_msg)
{
    if(!condition)
    {
        std::cout<<"error : "<<err_msg<<std::endl;
        exit(1);
    }
}

template<int N, typename T>
class vec
{
private:
    T data[N];
public:
    __host__  __device__
    vec(){}
    __host__  __device__
    vec(T v)
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] = v;
        }
    }
    __host__  __device__
    vec(T v0, T v1)
    {  
        static_assert(N == 2);
        data[0] = v0;
        data[1] = v1;
    }
    __host__  __device__
    vec(T v0, T v1, T v2)
    {  
        static_assert(N == 3);
        data[0] = v0;
        data[1] = v1;
        data[2] = v2;
    }
    __host__  __device__
    vec(const vec &v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] = v[i];
        }
    }
    __host__  __device__
    T &operator[](int i) 
    {
        return data[i];
    }
    __host__  __device__
    const T &operator[](int i) const
    {
        return data[i];
    }
    __host__  __device__
    vec operator+(const T v) const 
    {
        return vec(*this) += v;
    }
    __host__  __device__
    vec operator-(const T v) const 
    {
        return vec(*this) -= v;
    }
    __host__  __device__
    vec operator-() const
    {
        vec res;
        for(int i = 0; i < N; i++)
        {
            res.data[i] = -this->data[i];
        }
        return res;
    }
    __host__  __device__
    vec operator*(const T v) const 
    {
        return vec(*this) *= v;
    }
    __host__  __device__
    vec operator/(const T v) const 
    {
        return vec(*this) /= v;
    }
    __host__  __device__
    vec operator+(const vec v) const 
    {
        return vec(*this) += v;
    }
    __host__  __device__
    vec operator-(const vec v) const 
    {
        return vec(*this) -= v;
    }
    __host__  __device__
    vec operator*(const vec v) const 
    {
        return vec(*this) *= v;
    }
    __host__  __device__
    vec operator/(const vec v) const 
    {
        return vec(*this) /= v;
    }
    __host__  __device__
    const vec &operator+=(const T v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] += v;
        }
        return *this;
    }
    __host__  __device__
    const vec &operator-=(const T v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] -= v;
        }
        return *this;
    }
    __host__  __device__
    const vec &operator*=(const T v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] *= v;
        }
        return *this;
    }
    __host__  __device__
    const vec &operator/=(const T v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] /= v;
        }
        return *this;
    }
    __host__  __device__
    const vec &operator+=(const vec v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] += v.data[i];
        }
        return *this;
    }
    __host__  __device__
    const vec &operator-=(const vec v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] -= v.data[i];
        }
        return *this;
    }
    __host__  __device__
    const vec &operator*=(const vec v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] *= v.data[i];
        }
        return *this;
    }
    __host__  __device__
    const vec &operator/=(const vec v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] /= v.data[i];
        }
        return *this;
    }
    __host__  __device__
    bool operator==(const vec v) 
    {
        for(int i = 0; i < N; i++)
        {
            if(this->data[i] != v.data[i])
                return false;
        }
        return true;
    }
    __host__  __device__
    bool operator!=(const vec v) 
    {
        return !(*this == v);
    }
    __host__  __device__
    friend vec operator*(const T v, const vec& _vec)  // value * this
    {
        vec res;
        for(int i = 0; i < N; i++)
        {
            res[i] = v * _vec[i];
        }
        return res;
    }

    friend std::ostream &operator<<(std::ostream & os, const vec & v) 
    {
        os << "[";
        for(int i = 0; i < N; i++)
        {
            os << v.data[i];
            if(i < N - 1) os << ",";
        }
        os << "]";
        return os;
    }
    __host__  __device__
    T norm()
    {
        T s = 0;
        for(int i = 0; i < N; i++)
        {
            s += data[i] * data[i];
        }
        return s;
    }
    vec abs()
    {
        vec res;
        for(int i = 0; i < N; i++)
        {
            res.data[i] = fabs(this->data[i]);
        }
        return res;
    }
    __host__  __device__
    T r()
    {
        return sqrt(norm());
    }
    T sum()
    {
        T res = 0;
        for(int i = 0; i < N; i++)
        {
            res += this->data[i];
        }
        return res;
    }
    real mul()
    {
        real res = 1;
        for(int i = 0; i < N; i++)
        {
            res *= this->data[i];
        }
        return res;
    }
    vec<N, int> round()
    {
        vec<N,int> res;
        for(int i = 0; i < N; i++)
        {
            res[i] = std::round(this->data[i]);
        }
        return res;
    }
};

template<>
class vec<4, double>
{
private:
    union 
    {
        __m256d d;
        double data[4];
    };
public:
    vec(){}
    vec(double v)
    {
        d = _mm256_set1_pd(v);
    }
    vec(double v0, double v1, double v2, double v3)
    {  
        //although 'r' means reversed, 
        //the memory order is also reversed, 
        //so it is more intuitive than _mm256_set_pd
        d = _mm256_setr_pd(v0, v1, v2, v3);  
    }
    vec(const __m256d v) 
    {
        d = v;
    }
    vec(const vec &v) 
    {
        d = v.d;
    }
    double &operator[](int i) 
    {
        return data[i];
    }
    const double &operator[](int i) const
    {
        return data[i];
    }
    vec operator+(const double& v) const 
    {
        return vec(*this) += v;
    }
    vec operator-(const double& v) const 
    {
        return vec(*this) -= v;
    }
    vec operator-() const
    {
        return vec(_mm256_sub_pd(_mm256_set1_pd(0), d));
    }
    vec operator*(const double& v) const 
    {
        return vec(*this) *= v;
    }
    vec operator/(const double& v) const 
    {
        return vec(*this) /= v;
    }
    vec operator&(const vec& v) const 
    {
        return vec(*this) &= v;
    }
    vec operator+(const vec& v) const 
    {
        return vec(*this) += v;
    }
    vec operator-(const vec& v) const 
    {
        return vec(*this) -= v;
    }
    vec operator*(const vec& v) const 
    {
        return vec(*this) *= v;
    }
    vec operator/(const vec& v) const 
    {
        return vec(*this) /= v;
    }
    const vec &operator+=(const double& v) 
    {
        d = _mm256_add_pd(d, _mm256_set1_pd(v));
        return *this;
    }
    const vec &operator-=(const double& v) 
    {
        d = _mm256_sub_pd(d, _mm256_set1_pd(v));
        return *this;
    }
    const vec &operator*=(const double& v) 
    {
        d = _mm256_mul_pd(d, _mm256_set1_pd(v));
        return *this;
    }
    const vec &operator/=(const double& v) 
    {
        d = _mm256_div_pd(d, _mm256_set1_pd(v));
        return *this;
    }
    const vec &operator+=(const vec& v) 
    {
        d = _mm256_add_pd(d, v.d);
        return *this;
    }
    const vec &operator-=(const vec& v) 
    {
        d = _mm256_sub_pd(d, v.d);
        return *this;
    }
    const vec &operator*=(const vec& v) 
    {
        d = _mm256_mul_pd(d, v.d);
        return *this;
    }
    const vec &operator/=(const vec& v) 
    {
        d = _mm256_div_pd(d, v.d);
        return *this;
    }
    const vec &operator&=(const vec& v) 
    {
        d = _mm256_and_pd(d, v.d);
        return *this;
    }
    vec operator==(const vec& v) 
    {
        return vec(_mm256_cmp_pd(d,v.d,_CMP_EQ_OQ));
    }
    vec operator>(const vec& v) 
    {
        return vec(_mm256_cmp_pd(d,v.d,_CMP_GT_OQ));
    }
    vec operator>=(const vec& v) 
    {
        return vec(_mm256_cmp_pd(d,v.d,_CMP_GE_OQ));
    }
    vec operator<(const vec& v) 
    {
        return vec(_mm256_cmp_pd(d,v.d,_CMP_LT_OQ));
    }
    vec operator<=(const vec& v) 
    {
        return vec(_mm256_cmp_pd(d,v.d,_CMP_LE_OQ));
    }
    vec operator!=(const vec& v) 
    {
        return vec(_mm256_cmp_pd(d,v.d,_CMP_NEQ_OQ));
    }
    friend vec operator*(const double v, const vec& _vec)  // value * this
    {
        return vec(_mm256_mul_pd(_mm256_set1_pd(v), _vec.d));
    }

    friend std::ostream &operator<<(std::ostream & os, const vec & v) 
    {
        os << "__m256d[";
        os << v.data[0] << ",";
        os << v.data[1] << ",";
        os << v.data[2] << ",";
        os << v.data[3] << "]";
        return os;
    }
    double norm()
    {
        union
        {
            __m256d res;
            double data[4];
        };
        __m256d d2 = _mm256_mul_pd(d,d);
        res = _mm256_permute2f128_pd(d2,d2,1);;
        res = _mm256_add_pd(res,d2);
        res = _mm256_hadd_pd(res,res);
        return data[0];
    }
    double r()
    {
        return std::sqrt(norm());
    }
    double sum()
    {
        union
        {
            __m256d res;
            double data[4];
        };
        res = _mm256_permute2f128_pd(d,d,1);;
        res = _mm256_add_pd(res,d);
        res = _mm256_hadd_pd(res,res);
        return data[0];
    }
    vec sqrt()
    {
        return vec(_mm256_sqrt_pd(d));
    }
};



inline void print_simd(__m128d v)
{
    union {
            __m128d temp;
            double data[2];
        };
    temp = v;
    printf("(%.4f, %.4f)\n", data[0], data[1]);
}

inline void print_simd(__m256d v)
{
    union {
            __m256d temp;
            double data[4];
        };
    temp = v;
    printf("(%.4f, %.4f, %.4f, %.4f)\n", data[0], data[1], data[2], data[3]);
}

using vec2r = vec<2, real>;
using vec3r = vec<3, real>;
using vec2i = vec<2, int>;
using vec3i = vec<3, int>;
using vec4d = vec<4, double>;

template<typename T>
struct MatrixBase
{
    MatrixBase() {}
    MatrixBase(int rows, int cols)
    {
        m = rows;
        n = cols;
        d.resize(m * n);
    }
    int m;
    int n;
    std::vector<T> d;

    T& operator[](int idx)
    {
        return d[idx];
    }
};

using Matrix = MatrixBase<real>;
using Matriv = MatrixBase<vec3r>;

struct OffsetAndNumber
{
    int offset;
    int number;

    OffsetAndNumber(){}
    OffsetAndNumber(int offset_, int number_) : offset(offset_), number(number_){}

    friend std::ostream &operator<<(std::ostream & os, const OffsetAndNumber & on) 
    {
        os << "[" << on.offset << "," << on.number << "]";
        return os;
    }
};

using Range = OffsetAndNumber;


}