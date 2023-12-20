#pragma once
#include "type.h"
#include "body.h"
#include "cell.h"

namespace rtfmm
{
class LaplaceKernel
{
public:
    enum class KernelType
    {
        naive,
        matrix,
        fmm
    };
    LaplaceKernel();
    void p2p(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar, KernelType type);

private:
    void get_p2p_matrix(
        std::vector<vec3r>& x_src, 
        std::vector<vec3r>& x_tar, 
        std::vector<real>& matrix_p2p
    );
};
};