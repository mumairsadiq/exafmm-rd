#pragma once
#include "type.h"
#include "body.h"
#include "tree.h"

namespace rtfmm
{
class LaplaceKernel
{
public:
    LaplaceKernel();

    void p2p(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar);

    void p2p_matrix(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar);
    
    void p2m(int P, Bodies3& bs_src, Cell3& cell_src);

    void m2m(int P, Cell3& cell_src);

    void m2l(int P, Cell3& cell_src, Cell3& cell_tar);

    void l2l(int P, Cell3& cell_tar);

    void l2p(int P, Bodies3& bs_tar, Cell3& cell_tar);

private:
    Matrix get_p2p_matrix(
        std::vector<vec3r>& x_src, 
        std::vector<vec3r>& x_tar
    );

    Matriv get_force_naive(
        std::vector<vec3r> x_src, 
        std::vector<vec3r> x_tar, 
        Matrix q_src
    );
};
};