#pragma once
#include "type.h"
#include "body.h"
#include "tree.h"
#include <map>
#include "traverser.h"
#include <fftw3.h>
#include "align.h"
#include "surface.h"

namespace rtfmm
{
class LaplaceKernel
{
public:
    LaplaceKernel(){}

    /**
     * @brief simple direct method use AVX enable PBC
     * @note image cells are created of indices = [lower, dm] and [-dm, -lower]
    */
    void direct(Bodies3& bs_src, Bodies3& bs_tar, int images, real cycle, int lower = 0);

    /**
     * @brief direct method use AVX without PBC but offset
    */
    void direct_w(Bodies3& bs_src, Bodies3& bs_tar, std::vector<real>& ws_src, std::vector<real>& ws_tar, vec3r offset = vec3r(0,0,0));

    void direct(Cell3& cell_src, Cell3& cell_tar, int images, real cycle);

    void direct_reg(Cells3& cell_srcs, Cell3& cell_tar, std::vector<vec3r>& offsets);

    /**
     * @brief cell_src -> cell_tar P2P
     * @param offset offset of cell_src
     * @param use_simd true to use AVX
    */
    void p2p(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar, vec3r offset = vec3r(0,0,0), int use_simd = 1);

    /**
     * @brief (cells in p2ps) -> cell_tar P2P use SSE
    */
    void p2p_1toN_128(Bodies3& bs_src, Bodies3& bs_tar, Cells3& cs, std::vector<std::pair<int, vec3r>>& p2ps, Cell3& cell_tar);

    /**
     * @brief (cells in p2ps) -> cell_tar P2P use AVX
    */
    void p2p_1toN_256(Cells3& cs, std::vector<std::pair<int, vec3r>>& p2ps, Cell3& cell_tar);

    void p2m(int P, Bodies3& bs_src, Cell3& cell_src);

    void p2m_precompute(int P, Cell3& cell_src);

    void m2m(int P, Cell3& cell_parent, Cells3& cs);

    std::vector<real> m2m_ewald_root0(int P, Cell3& cell_parent, Cells3& cs, real cycle);

    void m2m_precompute(int P, Cell3& cell_parent, Cells3& cs);

    void m2m_img(int P, Cell3& cell_parent, Cells3& cs, real cycle);

    void m2m_img_precompute(int P, Cell3& cell_parent, Cells3& cs, real cycle);

    /**
     * @brief naive M2L performing O(N^2) potential computation
    */
    void m2l(int P, Cell3& cell_src, Cell3& cell_tar, vec3r offset = vec3r(0,0,0));

    /**
     * @brief M2L using naive FFT for each pair(cell_src -> cell_tar)
    */
    void m2l_fft(int P, Cell3& cell_src, Cell3& cell_tar, vec3r offset = vec3r(0,0,0));

    /**
     * @brief M2L using naive FFT for all pairs
    */
    void m2l_fft_precompute_naive(int P, Cells3& cs, PeriodicInteractionPairs& m2l_pairs);

    /**
     * @brief M2L using naive FFT for all pairs, reuse fftw_plan
    */
    void m2l_fft_precompute_advanced(int P, Cells3& cs, PeriodicInteractionPairs& m2l_pairs);

    /**
     * @brief M2L using FFT in the order of m2l_map.
     * Precompute all possible Gks, 
     * compute Qks for all src cells, 
     * then hadamard routine picks up Gk and Qk according to m2l, add Gk*Qk to corresponding Pk, 
     * finally convert Pk to all tar cells.
     * Since Gk for each real M2L is queried from map by m2l relative coordinate, 
     * we don't need to distinguish non-periodic/periodic cells in hadamard product stage(but the scale is different at the final store stage).
     * @note 2nd fast
    */
    void m2l_fft_precompute_advanced2(int P, Cells3& cs, PeriodicInteractionMap& m2l_map);

    /**
     * @brief Hadamard transposed version of m2l_fft_precompute_advanced2, frequency-wise parallelization.
     * @note 3rd fast
    */
    void m2l_fft_precompute_advanced3(int P, Cells3& cs, PeriodicInteractionMap& m2l_map, PeriodicInteractionPairs& m2l_pairs);

    /**
     * @brief M2L with FFT like exafmm-t.
     * Instead of computing M2L directly, it traverses neighbouring parent cell pairs, then compute M2L of children inside them,
     * namely, 8x8 child-to-child M2L interaction for non-periodic cells, 1x27 for periodic cells.
     * @note 1st fast
    */
    void m2l_fft_precompute_t(int P, Cells3& cs, PeriodicM2LMap& m2l_parent_map);

    void l2l(int P, Cell3& cell_parent, Cells3& cs);

    void l2l_precompute(int P, Cell3& cell_parent, Cells3& cs);

    void l2l_img_precompute(int P, Cell3& cell_parent, Cells3& cs);

    void l2p(int P, Bodies3& bs_tar, Cell3& cell_tar);

    void l2p_precompute(int P, Cell3& cell_tar);

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

    /**
     * @brief matrix vector storing 8x8 child-child interaction Gk of 2 neighbour cells
     * @note ### number of matrix = 26, size of each matrix = N_freq * 8 * 8 * complex
    */
    std::vector<AlignedVec> ccGks_8x8;

    /**
     * @brief matrix vector storing 1x7 child-child interaction Gk of 2 neighbour cells
     * @note ### number of matrix = 26, size of each matrix = N_freq * 1 * 28 * complex
     * 27 is padding to 28 for simd(because 28 % 4 = 0)
    */
    std::vector<AlignedVec> ccGks_1x27;

    void hadamard_8x8(
        int fft_size, 
        conv_grid_setting& cgrid,
        std::vector<size_t>& interaction_count_offset_8x8,
        std::vector<size_t>& interaction_offset_f_8x8,
        AlignedVec& Qk_all, 
        AlignedVec& Pk_all
    );

    void hadamard_1x27(
        int fft_size, 
        conv_grid_setting& cgrid,
        std::vector<size_t>& interaction_count_offset_1x27,
        std::vector<size_t>& interaction_offset_f_1x27,
        AlignedVec& Qk_all, 
        AlignedVec& Pk_all
    );
};
};