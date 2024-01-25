#pragma once
#include "type.h"
#include "body.h"
#include "tree.h"
#include <map>
#include "traverser.h"
#include <fftw3.h>
#include "align.h"

namespace rtfmm
{
    void dipole_correction(Bodies3& bs, real cycle);
class LaplaceKernel
{
public:
    LaplaceKernel();

    //static void p2p_p(std::vector<vec3r>& xs_src, std::vector<real>& qs_src, std::vector<vec3r>& xs_tar, std::vector<real>& ps_tar, const vec3r& offset = vec3r(0,0,0));

    static void p2p_pf(std::vector<vec3r>& xs_src, std::vector<real>& qs_src, std::vector<vec3r>& xs_tar, std::vector<real>& ps_tar, std::vector<vec3r>& fs_tar, const vec3r& offset = vec3r(0,0,0));

    void p2p(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar, vec3r offset = vec3r(0,0,0), int use_simd = 1);

    void p2p_1toN_128(Bodies3& bs_src, Bodies3& bs_tar, Cells3& cs, std::vector<std::pair<int, vec3r>>& p2ps, Cell3& cell_tar);

    void p2p_1toN_256(Bodies3& bs_src, Bodies3& bs_tar, Cells3& cs, std::vector<std::pair<int, vec3r>>& p2ps, Cell3& cell_tar);

    void p2m(int P, Bodies3& bs_src, Cell3& cell_src);

    void p2m_precompute(int P, Bodies3& bs_src, Cell3& cell_src);

    void m2m(int P, Cell3& cell_parent, Cells3& cs);

    void m2m_precompute(int P, Cell3& cell_parent, Cells3& cs);

    void m2m_img(int P, Cell3& cell_parent, Cells3& cs, real cycle);

    void m2m_img_precompute(int P, Cell3& cell_parent, Cells3& cs, real cycle);

    void m2l(int P, Cell3& cell_src, Cell3& cell_tar, vec3r offset = vec3r(0,0,0));

    void m2l_fft(int P, Cell3& cell_src, Cell3& cell_tar, vec3r offset = vec3r(0,0,0));

    void m2l_fft_precompute_naive(int P, Cells3& cs, PeriodicInteractionMap& m2l_map, PeriodicInteractionPairs& m2l_pairs);

    void m2l_fft_precompute_advanced(int P, Cells3& cs, PeriodicInteractionMap& m2l_map, PeriodicInteractionPairs& m2l_pairs);

    void m2l_fft_precompute_advanced2(int P, Cells3& cs, PeriodicInteractionMap& m2l_map);

    void m2l_fft_precompute_advanced3(int P, Cells3& cs, PeriodicInteractionMap& m2l_map, PeriodicInteractionPairs& m2l_pairs);

    void m2l_fft_precompute_t(int P, Cells3& cs, PeriodicInteractionMap& m2l_parent_map);

    void l2l(int P, Cell3& cell_parent, Cells3& cs);

    void l2l_precompute(int P, Cell3& cell_parent, Cells3& cs);

    void l2l_img_precompute(int P, Cell3& cell_parent, Cells3& cs);

    void l2p(int P, Bodies3& bs_tar, Cell3& cell_tar);

    void l2p_precompute(int P, Bodies3& bs_tar, Cell3& cell_tar);

    void m2p(int P, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar, vec3r offset = vec3r(0,0,0));

    void p2l(int P, Bodies3& bs_src, Cell3& cell_src, Cell3& cell_tar, vec3r offset = vec3r(0,0,0));

    void precompute(int P, real r0, int images);

    void precompute_m2l(int P, real r0, Cells3 cs, PeriodicInteractionMap m2l_map, int images);

    Matrix get_p2p_matrix(
        std::vector<vec3r>& x_src, 
        std::vector<vec3r>& x_tar
    );

    Matriv get_force_naive(
        std::vector<vec3r>& x_src, 
        std::vector<vec3r>& x_tar, 
        Matrix& q_src
    );

    std::vector<rtfmm::real> get_G_matrix(std::vector<rtfmm::vec3r>& grid, int N);
    std::vector<rtfmm::real> get_Q_matrix(Matrix& surface_q, int N, std::map<int,rtfmm::vec3i>& surf_conv_map);
    void get_Q_matrix(real* Q, Matrix& surface_q, int N, std::map<int,rtfmm::vec3i>& surf_conv_map);

    /**
     * @brief convert vec3i to hashed value
     * @return hash value, for example, (-1,-2,3) -> 1236, where 6=110 means negative for 1 and 2
    */
    static int hash(rtfmm::vec3i v)
    {
        int signs = ((v[0] < 0 ? 1 : 0) << 2) + ((v[1] < 0 ? 1 : 0) << 1) + (v[2] < 0 ? 1 : 0);
        int coord = std::abs(v[0]) * 1000 + std::abs(v[1]) * 100 + std::abs(v[2]) * 10;
        return coord + signs;
    }

    /**
     * @brief Get a hash(relx) -> index map.
     * For each coordinate x in [-range,+range]^3 >= minv, give it a unique index for hash(x).
     * @note The main usage of the returned map is to link an interaction pair to a precomputed interaction matrix.
     * For example, when coordinates of two cells in a M2L pair is given, say xsrc and xtar, the M2L matrix is map[hash(xsrc-xtar)].
    */
    static std::map<int,int> get_relx_idx_map(int range, int minv)
    {
        std::map<int, int> res;
        int cnt = 0;
        for(int k = -range; k <= range; k++)
        {
            for(int j = -range; j <= range; j++)
            {
                for(int i = -range; i <= range; i++)
                {
                    if(std::abs(i) >= minv || std::abs(j) >= minv || std::abs(k) >= minv)
                    {
                        res[hash(vec3i(i,j,k))] = cnt;
                        cnt++;
                    }
                }
            }
        }
        return res;
    }

    /**
     * @brief Get coorinates in [-range,+range]^3 >= minv.
     * @note when same parameters are passed, coordinates returned by this function and indices returned by get_relx_idx_map are one-to-one mathced.
     * For example, [1,-1,-1] is the 2th coordinates generated by get_relx(1,1), 
     * at the same time, index corresponding to hash([1,-1,-1]) is 2 for map from get_relx_idx_map(1,1)
    */
    template<typename T>
    static std::vector<vec<3,T>> get_relx(int range, int minv)
    {
        std::vector<vec<3,T>> res;
        for(int k = -range; k <= range; k++)
        {
            for(int j = -range; j <= range; j++)
            {
                for(int i = -range; i <= range; i++)
                {
                    if(std::abs(i) >= minv || std::abs(j) >= minv || std::abs(k) >= minv)
                    {
                        res.push_back(vec<3,T>(i,j,k));
                    }
                }
            }
        }
        return res;
    }

private:
    Matrix UT_p2m_precompute;
    Matrix V_p2m_precompute;
    Matrix Sinv_p2m_precompute;
    Matrix VSinv_p2m_precompute;

    Matrix UT_l2p_precompute;
    Matrix V_l2p_precompute;
    Matrix Sinv_l2p_precompute;
    Matrix VSinv_l2p_precompute;
    Matrix VSinvUT_l2p_precompute;

    Matrix matrix_m2m[8];
    Matrix matrix_m2m_img[27];
    Matrix matrix_l2l[8];
    Matrix matrix_l2l_img;

    std::vector<int> m2l_tars;
    std::vector<int> m2l_srcs;
    std::map<int, std::vector<complexr>> m2l_Gks;
    std::map<int, std::pair<std::vector<complexr>, int>> m2l_Gk_idx;

    void matmult_8x8x1(real*& M_, real*& IN0, real*& OUT0);
    
    void matmult_8x8x1_naive(real*& M_, real*& IN0, real*& OUT0);

    void matmult_8x8x2(real*& M_, real*& IN0, real*& IN1, real*& OUT0, real*& OUT1);

    void matmult_8x8x2_avx(double*& M_, double*& IN0, double*& IN1, double*& OUT0, double*& OUT1);

    /**
     * @brief matrix vector storing child-child interaction Gk of 2 neighbour cells
     * @note ### number of matrix = 26, size of each matrix = N_freq * 8 * 8 * complex
    */
    std::vector<AlignedVec> ccGks;
};
};