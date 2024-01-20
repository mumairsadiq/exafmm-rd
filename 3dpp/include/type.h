#pragma once
#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <chrono>
#include <assert.h>

#define TIME_BEGIN(a) auto time_begin_##a = std::chrono::high_resolution_clock::now()

#define TIME_END(a)   auto time_end_##a = std::chrono::high_resolution_clock::now();\
					  auto elapse_##a = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_##a - time_begin_##a);\
                      printf("[%s time measured : %.5f seconds.]\n", #a, elapse_##a.count() * 1e-9)

#define tbegin(a) TIME_BEGIN(a)
#define tend(a) TIME_END(a)

namespace rtfmm
{

inline void title(std::string s)
{
    std::cout<<"<"<<s<<">\n\n";
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
};

using vec2r = vec<2, real>;
using vec3r = vec<3, real>;
using vec3i = vec<3, int>;

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