#pragma once
#include "type.h"
#include <math.h>
#include <cblas.h>
#include <vector>

namespace rtfmm
{

enum class MathType
{
    naive,
    blas,
    cublas
};
/**
 * @brief matrix X vec.
 * 
 * @param m rows
 * @param n cols
 * @param A m by n matrix
 * @param b n by 1 vector
 * @param c m by 1 vector
 * @param type underlaying math function type
 */
void mat_vec(int m, int n, std::vector<real>& A, std::vector<real>& b, std::vector<real>& c, MathType type = MathType::naive); 
}