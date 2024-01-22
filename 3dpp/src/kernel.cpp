#include "kernel.h"
#include "mathfunc.h"
#include "surface.h"
#include "argument.h"
#include <omp.h>
#include <set>
#include "align.h"
#include "emmintrin.h"
#include "immintrin.h"

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


rtfmm::LaplaceKernel::LaplaceKernel()
{
    
}

static void print_simd(__m128d v)
{
    union {
            __m128d temp;
            double data[2];
        };
    temp = v;
    printf("(%.4f, %.4f)\n", data[0], data[1]);
}

void rtfmm::LaplaceKernel::p2p(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar, vec3r offset, int use_simd)
{
    if(!use_simd)
    {
        for(int j = 0; j < cell_tar.brange.number; j++)
        {
            Body3& btar = bs_tar[cell_tar.brange.offset + j];
            real potential = 0;
            vec3r force(0,0,0);
            for(int i = 0; i < cell_src.brange.number; i++)
            {
                Body3& bsrc = bs_src[cell_src.brange.offset + i];
                vec3r dx = btar.x - bsrc.x - offset;
                real r = dx.r();
                real invr = r == 0 ? 0 : 1 / r;
                potential += bsrc.q * invr;
                force += bsrc.q * invr * invr * invr * (-dx);
            }
            btar.p += potential;
            btar.f += force;
        }
    }
    else
    {
        __m128d zero = _mm_setr_pd(0,0);
        __m128d one = _mm_setr_pd(1,1);
        for(int j = 0; j < cell_tar.brange.number; j+=2)
        {
            Body3* btar = &bs_tar[cell_tar.brange.offset + j];
            __m128d xj,yj,zj;
            if(j >= cell_tar.brange.number - 1)
            {
                xj = _mm_setr_pd(btar->x[0],0);
                yj = _mm_setr_pd(btar->x[1],0);
                zj = _mm_setr_pd(btar->x[2],0);
            }
            else
            {
                xj = _mm_setr_pd(btar->x[0],(btar+1)->x[0]);
                yj = _mm_setr_pd(btar->x[1],(btar+1)->x[1]);
                zj = _mm_setr_pd(btar->x[2],(btar+1)->x[2]);
            }
            
            union{
                __m128d pj = _mm_setr_pd(0,0);
                double pdata[2];
            };
            union{
                __m128d fxj = _mm_setr_pd(0,0);
                double fxjdata[2];
            };
            union{
                __m128d fyj = _mm_setr_pd(0,0);
                double fyjdata[2];
            };
            union{
                __m128d fzj = _mm_setr_pd(0,0);
                double fzjdata[2];
            };

            __m128d xi, yi, zi, qi, invr, f, qinvr;

            for(int i = 0; i < cell_src.brange.number; i++)
            {
                Body3& bsrc = bs_src[cell_src.brange.offset + i];
                vec3r bsrc_xi = bsrc.x + offset;
                xi = _mm_setr_pd(bsrc_xi[0],bsrc_xi[0]);
                yi = _mm_setr_pd(bsrc_xi[1],bsrc_xi[1]);
                zi = _mm_setr_pd(bsrc_xi[2],bsrc_xi[2]);
                qi = _mm_setr_pd(bsrc.q,bsrc.q);
                xi = _mm_sub_pd(xj, xi);
                yi = _mm_sub_pd(yj, yi);
                zi = _mm_sub_pd(zj, zi);
                invr = _mm_mul_pd(xi, xi);
                invr = _mm_add_pd(invr, _mm_mul_pd(yi, yi));
                invr = _mm_add_pd(invr, _mm_mul_pd(zi, zi));
                f = _mm_cmpgt_pd(invr, zero);
                invr = _mm_sqrt_pd(invr);
                invr = _mm_div_pd(one, invr);
                invr = _mm_and_pd(invr, f);
                qinvr = _mm_mul_pd(qi, invr);
                pj = _mm_add_pd(pj, qinvr);
                qinvr = _mm_mul_pd(qinvr, invr);
                qinvr = _mm_mul_pd(qinvr, invr);
                fxj = _mm_add_pd(fxj, _mm_mul_pd(qinvr, xi));
                fyj = _mm_add_pd(fyj, _mm_mul_pd(qinvr, yi));
                fzj = _mm_add_pd(fzj, _mm_mul_pd(qinvr, zi));
            }
            btar->p += pdata[0];
            btar->f[0] -= fxjdata[0];
            btar->f[1] -= fyjdata[0];
            btar->f[2] -= fzjdata[0];
            if(j < cell_tar.brange.number - 1)
            {
                (btar + 1)->p += pdata[1];
                (btar + 1)->f[0] -= fxjdata[1];
                (btar + 1)->f[1] -= fyjdata[1];
                (btar + 1)->f[2] -= fzjdata[1];
            }
        }
    }
}

void rtfmm::LaplaceKernel::p2p_1toN(Bodies3& bs_src, Bodies3& bs_tar, Cells3& cs, std::vector<std::pair<int, vec3r>>& p2ps, Cell3& cell_tar)
{
    #if 0
    __m128d zero = _mm_set1_pd(0);
    __m128d one = _mm_set1_pd(1);
    for(int j = 0; j < cell_tar.brange.number; j+=2)
    {
        Body3* btar = &bs_tar[cell_tar.brange.offset + j];
        __m128d xj,yj,zj;
        if(j >= cell_tar.brange.number - 1)
        {
            xj = _mm_set1_pd(btar->x[0]);
            yj = _mm_set1_pd(btar->x[1]);
            zj = _mm_set1_pd(btar->x[2]);
        }
        else
        {
            xj = _mm_setr_pd(btar->x[0],(btar+1)->x[0]);
            yj = _mm_setr_pd(btar->x[1],(btar+1)->x[1]);
            zj = _mm_setr_pd(btar->x[2],(btar+1)->x[2]);
        }
        
        union{
            __m128d pj = _mm_set1_pd(0);
            double pdata[2];
        };
        union{
            __m128d fxj = _mm_set1_pd(0);
            double fxjdata[2];
        };
        union{
            __m128d fyj = _mm_set1_pd(0);
            double fyjdata[2];
        };
        union{
            __m128d fzj = _mm_set1_pd(0);
            double fzjdata[2];
        };
        __m128d xi, yi, zi, qi, invr, f, qinvr;
        for(int si = 0; si < p2ps.size(); si++)
        {
            Cell3& cell_src = cs[p2ps[si].first];
            vec3r& offset = p2ps[si].second;
            for(int i = 0; i < cell_src.brange.number; i++)
            {
                Body3& bsrc = bs_src[cell_src.brange.offset + i];
                vec3r bsrc_xi = bsrc.x + offset;
                xi = _mm_set1_pd(bsrc_xi[0]);
                yi = _mm_set1_pd(bsrc_xi[1]);
                zi = _mm_set1_pd(bsrc_xi[2]);
                qi = _mm_set1_pd(bsrc.q);
                xi = _mm_sub_pd(xj, xi);
                yi = _mm_sub_pd(yj, yi);
                zi = _mm_sub_pd(zj, zi);
                invr = _mm_mul_pd(xi, xi);
                invr = _mm_add_pd(invr, _mm_mul_pd(yi, yi));
                invr = _mm_add_pd(invr, _mm_mul_pd(zi, zi));
                f = _mm_cmpgt_pd(invr, zero);
                invr = _mm_sqrt_pd(invr);
                invr = _mm_div_pd(one, invr);
                invr = _mm_and_pd(invr, f);
                qinvr = _mm_mul_pd(qi, invr);
                pj = _mm_add_pd(pj, qinvr);
                qinvr = _mm_mul_pd(qinvr, invr);
                qinvr = _mm_mul_pd(qinvr, invr);
                fxj = _mm_add_pd(fxj, _mm_mul_pd(qinvr, xi));
                fyj = _mm_add_pd(fyj, _mm_mul_pd(qinvr, yi));
                fzj = _mm_add_pd(fzj, _mm_mul_pd(qinvr, zi));
            }
        }
        btar->p += pdata[0];
        btar->f[0] -= fxjdata[0];
        btar->f[1] -= fyjdata[0];
        btar->f[2] -= fzjdata[0];
        if(j < cell_tar.brange.number - 1)
        {
            (btar + 1)->p += pdata[1];
            (btar + 1)->f[0] -= fxjdata[1];
            (btar + 1)->f[1] -= fyjdata[1];
            (btar + 1)->f[2] -= fzjdata[1];
        }
    }
    #else
    __m256d zero = _mm256_set1_pd(0);
    __m256d one = _mm256_set1_pd(1);
    for(int j = 0; j < cell_tar.brange.number; j+=4)
    {
        Body3* btar = &bs_tar[cell_tar.brange.offset + j];
        __m256d xj,yj,zj;
        int padding = j + 4 - cell_tar.brange.number;
        padding = padding <= 0 ? 0 : padding;
        if(padding <= 0)
        {
            xj = _mm256_setr_pd(btar->x[0],(btar+1)->x[0],(btar+2)->x[0],(btar+3)->x[0]);
            yj = _mm256_setr_pd(btar->x[1],(btar+1)->x[1],(btar+2)->x[1],(btar+3)->x[1]);
            zj = _mm256_setr_pd(btar->x[2],(btar+1)->x[2],(btar+2)->x[2],(btar+3)->x[2]);
        }
        else if(padding == 1)
        {
            xj = _mm256_setr_pd(btar->x[0],(btar+1)->x[0],(btar+2)->x[0],0);
            yj = _mm256_setr_pd(btar->x[1],(btar+1)->x[1],(btar+2)->x[1],0);
            zj = _mm256_setr_pd(btar->x[2],(btar+1)->x[2],(btar+2)->x[2],0);
        }
        else if(padding == 2)
        {
            xj = _mm256_setr_pd(btar->x[0],(btar+1)->x[0],0,0);
            yj = _mm256_setr_pd(btar->x[1],(btar+1)->x[1],0,0);
            zj = _mm256_setr_pd(btar->x[2],(btar+1)->x[2],0,0);
        }
        else if(padding == 3)
        {
            xj = _mm256_setr_pd(btar->x[0],0,0,0);
            yj = _mm256_setr_pd(btar->x[1],0,0,0);
            zj = _mm256_setr_pd(btar->x[2],0,0,0);
        }
        
        union{
            __m256d pj = _mm256_set1_pd(0);
            double pdata[2];
        };
        union{
            __m256d fxj = _mm256_set1_pd(0);
            double fxjdata[2];
        };
        union{
            __m256d fyj = _mm256_set1_pd(0);
            double fyjdata[2];
        };
        union{
            __m256d fzj = _mm256_set1_pd(0);
            double fzjdata[2];
        };
        __m256d xi, yi, zi, qi, invr, f, qinvr;
        for(int si = 0; si < p2ps.size(); si++)
        {
            Cell3& cell_src = cs[p2ps[si].first];
            vec3r& offset = p2ps[si].second;
            for(int i = 0; i < cell_src.brange.number; i++)
            {
                Body3& bsrc = bs_src[cell_src.brange.offset + i];
                vec3r bsrc_xi = bsrc.x + offset;
                xi = _mm256_set1_pd(bsrc_xi[0]);
                yi = _mm256_set1_pd(bsrc_xi[1]);
                zi = _mm256_set1_pd(bsrc_xi[2]);
                qi = _mm256_set1_pd(bsrc.q);
                xi = _mm256_sub_pd(xj, xi);
                yi = _mm256_sub_pd(yj, yi);
                zi = _mm256_sub_pd(zj, zi);
                invr = _mm256_mul_pd(xi, xi);
                invr = _mm256_add_pd(invr, _mm256_mul_pd(yi, yi));
                invr = _mm256_add_pd(invr, _mm256_mul_pd(zi, zi));
                f = _mm256_cmp_pd(invr, zero, _CMP_GT_OQ);
                invr = _mm256_sqrt_pd(invr);
                invr = _mm256_div_pd(one, invr);
                invr = _mm256_and_pd(invr, f);
                qinvr = _mm256_mul_pd(qi, invr);
                pj = _mm256_add_pd(pj, qinvr);
                qinvr = _mm256_mul_pd(qinvr, invr);
                qinvr = _mm256_mul_pd(qinvr, invr);
                fxj = _mm256_add_pd(fxj, _mm256_mul_pd(qinvr, xi));
                fyj = _mm256_add_pd(fyj, _mm256_mul_pd(qinvr, yi));
                fzj = _mm256_add_pd(fzj, _mm256_mul_pd(qinvr, zi));
            }
        }
        for(int k = 0; k < 4 - padding; k++)
        {
            (btar + k)->p    += pdata[k];
            (btar + k)->f[0] -= fxjdata[k];
            (btar + k)->f[1] -= fyjdata[k];
            (btar + k)->f[2] -= fzjdata[k];
        }
    }
    #endif
}

void rtfmm::LaplaceKernel::p2p_1toN(ManyBody& bs_src, ManyBody& bs_tar, Cells3& cs, std::vector<std::pair<int, vec3r>>& p2ps, Cell3& cell_tar)
{
    __m128d zero = _mm_set1_pd(0);
    __m128d one = _mm_set1_pd(1);
    for(int j = 0; j < cell_tar.brange.number; j+=2)
    {
        int tidx = cell_tar.brange.offset + j;
        __m128d xj,yj,zj;
        if(j >= cell_tar.brange.number - 1)
        {
            xj = _mm_set1_pd(bs_tar.xs[tidx]);
            yj = _mm_set1_pd(bs_tar.ys[tidx]);
            zj = _mm_set1_pd(bs_tar.zs[tidx]);
        }
        else
        {
            xj = _mm_setr_pd(bs_tar.xs[tidx],bs_tar.xs[tidx+1]);
            yj = _mm_setr_pd(bs_tar.ys[tidx],bs_tar.ys[tidx+1]);
            zj = _mm_setr_pd(bs_tar.zs[tidx],bs_tar.zs[tidx+1]);
        }
        
        union{
            __m128d pj = _mm_set1_pd(0);
            double pdata[2];
        };
        union{
            __m128d fxj = _mm_set1_pd(0);
            double fxjdata[2];
        };
        union{
            __m128d fyj = _mm_set1_pd(0);
            double fyjdata[2];
        };
        union{
            __m128d fzj = _mm_set1_pd(0);
            double fzjdata[2];
        };
        __m128d xi, yi, zi, qi, invr, f, qinvr;
        for(int si = 0; si < p2ps.size(); si++)
        {
            Cell3& cell_src = cs[p2ps[si].first];
            vec3r& offset = p2ps[si].second;
            for(int i = 0; i < cell_src.brange.number; i++)
            {
                int sidx = cell_src.brange.offset + i;
                xi = _mm_set1_pd(bs_src.xs[sidx]+offset[0]);
                yi = _mm_set1_pd(bs_src.ys[sidx]+offset[1]);
                zi = _mm_set1_pd(bs_src.zs[sidx]+offset[2]);
                qi = _mm_set1_pd(bs_src.qs[sidx]);
                xi = _mm_sub_pd(xj, xi);
                yi = _mm_sub_pd(yj, yi);
                zi = _mm_sub_pd(zj, zi);
                invr = _mm_mul_pd(xi, xi);
                invr = _mm_add_pd(invr, _mm_mul_pd(yi, yi));
                invr = _mm_add_pd(invr, _mm_mul_pd(zi, zi));
                f = _mm_cmpgt_pd(invr, zero);
                invr = _mm_sqrt_pd(invr);
                invr = _mm_div_pd(one, invr);
                invr = _mm_and_pd(invr, f);
                qinvr = _mm_mul_pd(qi, invr);
                pj = _mm_add_pd(pj, qinvr);
                qinvr = _mm_mul_pd(qinvr, invr);
                qinvr = _mm_mul_pd(qinvr, invr);
                fxj = _mm_add_pd(fxj, _mm_mul_pd(qinvr, xi));
                fyj = _mm_add_pd(fyj, _mm_mul_pd(qinvr, yi));
                fzj = _mm_add_pd(fzj, _mm_mul_pd(qinvr, zi));
            }
        }
        bs_tar.ps[tidx] += pdata[0];
        bs_tar.fxs[tidx] -= fxjdata[0];
        bs_tar.fys[tidx] -= fyjdata[0];
        bs_tar.fzs[tidx] -= fzjdata[0];
        if(j < cell_tar.brange.number - 1)
        {
            bs_tar.ps[tidx+1] += pdata[1];
            bs_tar.fxs[tidx+1] -= fxjdata[1];
            bs_tar.fys[tidx+1] -= fyjdata[1];
            bs_tar.fzs[tidx+1] -= fzjdata[1];
        }
    }
}

void rtfmm::LaplaceKernel::p2m(int P, Bodies3& bs_src, Cell3& cell_src)
{
    /* get source to check matrix */
    std::vector<vec3r> x_check = get_surface_points(P, cell_src.r * 2.95, cell_src.x);
    std::vector<vec3r> x_src = get_bodies_x(bs_src, cell_src.brange);
    Matrix s2c = get_p2p_matrix(x_src, x_check);

    /* get check potential */
    Matrix q_src = get_bodies_q(bs_src, cell_src.brange);
    Matrix p_check = mat_vec_mul(s2c, q_src);

    /* get equivalent to check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_src.r * 1.05, cell_src.x);
    Matrix e2c = get_p2p_matrix(x_equiv, x_check);

    /* get equivalent charge */
    cell_src.q_equiv = linear_equation_system_svd(e2c, p_check);
}

void rtfmm::LaplaceKernel::p2m_precompute(int P, Bodies3& bs_src, Cell3& cell_src)
{
    int num_surf = get_surface_point_num(P);
    std::vector<vec3r> x_check = get_surface_points(P, cell_src.r * 2.95, cell_src.x);
    Matrix p_check(num_surf, 1);
    for(int j = 0; j < num_surf; j++)
    {
        vec3r& xj = x_check[j];
        real potential = 0;
        for(int i = 0; i < cell_src.brange.number; i++)
        {
            Body3& bi = bs_src[cell_src.brange.offset + i];
            real r = (xj - bi.x).r();
            real invr = r == 0 ? 0 : 1 / r;
            potential += bi.q * invr;
        }
        p_check[j] = potential;
    }

    /* get equivalent charge */
    real scale = std::pow(0.5, cell_src.depth);
    Matrix buffer(num_surf, 1);
    mat_vec_mul(UT_p2m_precompute, p_check, buffer);
    mat_vec_mul(VSinv_p2m_precompute, buffer, cell_src.q_equiv, scale);
}


void rtfmm::LaplaceKernel::m2m(int P, Cell3& cell_parent, Cells3& cs)
{
    /* parent equiv to parent check matrix */
    std::vector<rtfmm::vec3r> x_equiv_parent = get_surface_points(P, cell_parent.r * 1.05, cell_parent.x);
    std::vector<rtfmm::vec3r> x_check_parent = get_surface_points(P, cell_parent.r * 2.95, cell_parent.x);
    Matrix pe2pc = get_p2p_matrix(x_equiv_parent, x_check_parent);

    for(int i = 0; i < cell_parent.crange.number; i++)
    {
        Cell3& cell_child = cs[cell_parent.crange.offset + i];

        /* child equivalent to parent check */
        std::vector<rtfmm::vec3r> x_equiv_child = get_surface_points(P, cell_child.r * 1.05, cell_child.x);
        Matrix ce2pc = get_p2p_matrix(x_equiv_child, x_check_parent);
        Matrix p_check_parent = mat_vec_mul(ce2pc, cell_child.q_equiv);

        /* parent check potential to parent equivalent charge */
        Matrix q_equiv_parent = linear_equation_system_svd(pe2pc, p_check_parent);
        cell_parent.q_equiv = mat_mat_add(cell_parent.q_equiv, q_equiv_parent);
    }
}

void rtfmm::LaplaceKernel::m2m_precompute(int P, Cell3& cell_parent, Cells3& cs)
{
    for(int i = 0; i < cell_parent.crange.number; i++)
    {
        Cell3& cell_child = cs[cell_parent.crange.offset + i];
        Matrix q_equiv_parent(get_surface_point_num(P), 1);
        mat_vec_mul(matrix_m2m[cell_child.octant], cell_child.q_equiv, q_equiv_parent);
        mat_mat_increment(cell_parent.q_equiv, q_equiv_parent);
    }
}

void rtfmm::LaplaceKernel::m2m_img(int P, Cell3& cell_parent, Cells3& cs, real cycle)
{
    /* parent equiv to parent check matrix */
    std::vector<rtfmm::vec3r> x_equiv_parent = get_surface_points(P, cell_parent.r * 1.05, cell_parent.x);
    std::vector<rtfmm::vec3r> x_check_parent = get_surface_points(P, cell_parent.r * 2.95, cell_parent.x);
    Matrix pe2pc = get_p2p_matrix(x_equiv_parent, x_check_parent);

    for(int pz = -1; pz <= 1; pz++)
    {
        for(int py = -1; py <= 1; py++)
        {
            for(int px = -1; px <= 1; px++)
            {
                Cell3& cell_child = cs[cell_parent.crange.offset];

                /* child equivalent to parent check */
                std::vector<rtfmm::vec3r> x_equiv_child = get_surface_points(P, cell_child.r * 1.05, cell_child.x + vec3r(px,py,pz) * cycle);
                Matrix ce2pc = get_p2p_matrix(x_equiv_child, x_check_parent);
                Matrix p_check_parent = mat_vec_mul(ce2pc, cell_child.q_equiv);

                /* parent check potential to parent equivalent charge */
                Matrix q_equiv_parent = linear_equation_system_svd(pe2pc, p_check_parent);
                cell_parent.q_equiv = mat_mat_add(cell_parent.q_equiv, q_equiv_parent);
            }
        }
    }
}


void rtfmm::LaplaceKernel::m2m_img_precompute(int P, Cell3& cell_parent, Cells3& cs, real cycle)
{
    Cell3& cell_child = cs[cell_parent.crange.offset];
    for(int octant = 0; octant < 27; octant++)
    {
        Matrix q_equiv_parent = mat_vec_mul(matrix_m2m_img[octant], cell_child.q_equiv);
        cell_parent.q_equiv = mat_mat_add(cell_parent.q_equiv, q_equiv_parent);
    }
}


void rtfmm::LaplaceKernel::m2l(int P, Cell3& cell_src, Cell3& cell_tar, vec3r offset)
{
    /* get src equivalent to tar check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_src.r * 1.05, cell_src.x + offset);
    std::vector<rtfmm::vec3r> x_check = get_surface_points(P, cell_tar.r * 1.05, cell_tar.x);
    Matrix e2c = get_p2p_matrix(x_equiv, x_check);

    // source cell equivalent charge to target cell check potential
    Matrix p_check = mat_vec_mul(e2c, cell_src.q_equiv);
    cell_tar.p_check = mat_mat_add(cell_tar.p_check, p_check);
}

void rtfmm::LaplaceKernel::m2l_fft(int P, Cell3& cell_src, Cell3& cell_tar, vec3r offset)
{
    assert_exit(cell_src.depth == cell_tar.depth, "depth not equal");
    real r = cell_src.r * 1.05;
    vec3r relative_pos = cell_src.x + offset - cell_tar.x;
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, r, relative_pos);

    rtfmm::real delta = 2 * r / (P - 1);
    int N = 2 * P - 1;
    rtfmm::real gsize = r * 2;
    std::vector<rtfmm::vec3r> grid = get_conv_grid(N, gsize, delta, relative_pos);
    std::vector<rtfmm::real> G = get_G_matrix(grid, N);
    std::map<int,rtfmm::vec3i> surf_conv_map = get_surface_conv_map(P);
    std::vector<rtfmm::real> Q = get_Q_matrix(cell_src.q_equiv, N, surf_conv_map);

    std::vector<fftw_complex> Gk(N * N * N);
    std::vector<fftw_complex> Qk(N * N * N);
    std::vector<fftw_complex> pk(N * N * N);
    std::vector<rtfmm::real> p_fft_grid(N * N * N);
    fftw_plan plan_G = fftw_plan_dft_r2c_3d(N, N, N, G.data(), Gk.data(), FFTW_ESTIMATE);
    fftw_plan plan_Q = fftw_plan_dft_r2c_3d(N, N, N, Q.data(), Qk.data(), FFTW_ESTIMATE);
    fftw_execute(plan_G);
    fftw_execute(plan_Q);
    for(int i = 0; i < N * N * N; i++)
    {
        rtfmm::real a = Gk[i][0];
        rtfmm::real b = Gk[i][1];
        rtfmm::real c = Qk[i][0];
        rtfmm::real d = Qk[i][1];
        pk[i][0] = a * c - b * d;
        pk[i][1] = a * d + b * c;
    }
    fftw_plan plan_P = fftw_plan_dft_c2r_3d(N, N, N, pk.data(), p_fft_grid.data(), FFTW_ESTIMATE);
    fftw_execute(plan_P);
    for(int i = 0; i < N * N * N; i++)
    {
        p_fft_grid[i] /= N * N * N;
    }

    int surface_number = rtfmm::get_surface_point_num(P);
    Matrix ps_fft(surface_number, 1);
    for(int i = 0; i < surface_number; i++)
    {
        rtfmm::vec3i idx3 = surf_conv_map[i] + rtfmm::vec3i(P-1,P-1,P-1);
        ps_fft.d[i] = p_fft_grid[idx3[2] * N * N + idx3[1] * N + idx3[0]];
    }

    fftw_destroy_plan(plan_P);
    fftw_destroy_plan(plan_G);
    fftw_destroy_plan(plan_Q);

    cell_tar.p_check = mat_mat_add(cell_tar.p_check, ps_fft);
}

void rtfmm::LaplaceKernel::m2l_fft_precompute_naive(int P, Cells3& cs, PeriodicInteractionMap& m2l_map, PeriodicInteractionPairs& m2l_pairs)
{
    printf("m2l_fft_precompute_naive\n");
    int num_fft = m2l_pairs.size();
    int N = 2 * P - 1;
    int N3 = N * N * N;
    int surface_number = rtfmm::get_surface_point_num(P);
    std::map<int,rtfmm::vec3i> surf_conv_map = get_surface_conv_map(P);

    TIME_BEGIN(init_mem);
    std::vector<std::vector<real>> G_all(num_fft, std::vector<real>(N3));
    std::vector<std::vector<real>> Q_all(num_fft, std::vector<real>(N3));
    std::vector<std::vector<real>> p_all(num_fft, std::vector<real>(N3));
    std::vector<std::vector<complexr>> Gk_all(num_fft, std::vector<complexr>(N3));
    std::vector<std::vector<complexr>> Qk_all(num_fft, std::vector<complexr>(N3));
    std::vector<std::vector<complexr>> pk_all(num_fft, std::vector<complexr>(N3));
    std::vector<fftw_plan> plans_G(num_fft);
    std::vector<fftw_plan> plans_Q(num_fft);
    std::vector<fftw_plan> plans_P(num_fft);
    std::vector<int> indices_tar(num_fft);
    TIME_END(init_mem);

    TIME_BEGIN(setup);
    for(int i = 0; i < num_fft; i++)
    {
        PeriodicInteractionPair m2l = m2l_pairs[i];
        Cell3& cell_tar = cs[m2l.first];
        Cell3& cell_src = cs[m2l.second.first];
        vec3r offset_src = m2l.second.second;
        real r = cell_src.r * 1.05;
        vec3r relative_pos = cell_src.x + offset_src - cell_tar.x;
        std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, r, relative_pos);
        rtfmm::real delta = 2 * r / (P - 1);
        rtfmm::real gsize = r * 2;
        std::vector<rtfmm::vec3r> grid = get_conv_grid(N, gsize, delta, relative_pos);
        std::vector<rtfmm::real>& G = G_all[i];
        std::vector<rtfmm::real>& Q = Q_all[i];
        std::vector<rtfmm::real>& p_fft_grid = p_all[i];
        std::vector<complexr>& Gk = Gk_all[i];
        std::vector<complexr>& Qk = Qk_all[i];
        std::vector<complexr>& pk = pk_all[i];
        G = get_G_matrix(grid, N);
        Q = get_Q_matrix(cell_src.q_equiv, N, surf_conv_map);
        plans_G[i] = fftw_plan_dft_r2c_3d(N, N, N, G.data(), reinterpret_cast<fftw_complex*>(Gk.data()), FFTW_ESTIMATE);
        plans_Q[i] = fftw_plan_dft_r2c_3d(N, N, N, Q.data(), reinterpret_cast<fftw_complex*>(Qk.data()), FFTW_ESTIMATE);
        plans_P[i] = fftw_plan_dft_c2r_3d(N, N, N, reinterpret_cast<fftw_complex*>(pk.data()), p_fft_grid.data(), FFTW_ESTIMATE);
        indices_tar[i] = m2l.first;
    }
    TIME_END(setup);
    
    int num_plan = plans_G.size();
    printf("num_plan = %d\n", num_plan);

    TIME_BEGIN(GQ);
    #pragma omp parallel for
    for(int i = 0; i < num_plan; i++)
    {
        fftw_plan& plan_G = plans_G[i];
        fftw_plan& plan_Q = plans_Q[i];
        fftw_execute(plan_G);
        fftw_execute(plan_Q);
        fftw_destroy_plan(plan_G);
        fftw_destroy_plan(plan_Q);
    }
    TIME_END(GQ);

    TIME_BEGIN(P);
    #pragma omp parallel for
    for(int i = 0; i < num_plan; i++)
    {
        std::vector<complexr>& Gk = Gk_all[i];
        std::vector<complexr>& Qk = Qk_all[i];
        std::vector<complexr>& pk = pk_all[i];
        fftw_plan& plan_P = plans_P[i];

        for(int i = 0; i < N * N * N; i++)
        {
            rtfmm::real a = Gk[i].r;
            rtfmm::real b = Gk[i].i;
            rtfmm::real c = Qk[i].r;
            rtfmm::real d = Qk[i].i;
            pk[i].r = a * c - b * d;
            pk[i].i = a * d + b * c;
        }
        fftw_execute(plan_P);
        fftw_destroy_plan(plan_P);
    }
    TIME_END(P);

    for(int i = 0; i < num_plan; i++)
    {
        Cell3& cell_tar = cs[indices_tar[i]];
        std::vector<rtfmm::real>& p_fft_grid = p_all[i];
        for(int i = 0; i < surface_number; i++)
        {
            rtfmm::vec3i idx3 = surf_conv_map[i] + rtfmm::vec3i(P-1,P-1,P-1);
            cell_tar.p_check.d[i] += p_fft_grid[idx3[2] * N * N + idx3[1] * N + idx3[0]] / (N * N * N);
        }
    }
}

void rtfmm::LaplaceKernel::m2l_fft_precompute_advanced(int P, Cells3& cs, PeriodicInteractionMap& m2l_map, PeriodicInteractionPairs& m2l_pairs)
{
    printf("m2l_fft_precompute_advanced\n");
    TIME_BEGIN(init_map);
    int num_fft = m2l_pairs.size();
    int N = 2 * P - 1;
    int N3 = N * N * N;
    int surface_number = rtfmm::get_surface_point_num(P);
    std::map<int,rtfmm::vec3i> surf_conv_map = get_surface_conv_map(P);
    
    printf("num_fft = %d\n", num_fft);

    int num_src = m2l_srcs.size();
    std::map<int,int> m2l_src_map;
    for(int i = 0; i < num_src; i++)
    {
        m2l_src_map[m2l_srcs[i]] = i;
    }
    TIME_END(init_map);

    TIME_BEGIN(init_mem);
    std::vector<std::vector<real>> G_all(num_fft, std::vector<real>(N3));
    std::vector<std::vector<real>> Q_all(num_src, std::vector<real>(N3));
    std::vector<std::vector<real>> p_all(num_fft, std::vector<real>(N3));
    std::vector<std::vector<complexr>> Gk_all(num_fft, std::vector<complexr>(N3));
    std::vector<std::vector<complexr>> Qk_all(num_src, std::vector<complexr>(N3));
    std::vector<std::vector<complexr>> pk_all(num_fft, std::vector<complexr>(N3));
    TIME_END(init_mem);

    TIME_BEGIN(setup);
    #pragma omp parallel for
    for(int i = 0; i < num_src; i++)
    {
        Cell3& cell_src = cs[m2l_srcs[i]];
        std::vector<rtfmm::real>& Q = Q_all[i];
        Q = get_Q_matrix(cell_src.q_equiv, N, surf_conv_map);
    }
    #pragma omp parallel for
    for(int i = 0; i < num_fft; i++)
    {
        PeriodicInteractionPair m2l = m2l_pairs[i];
        Cell3& cell_tar = cs[m2l.first];
        Cell3& cell_src = cs[m2l.second.first];
        vec3r offset_src = m2l.second.second;
        real r = cell_src.r * 1.05;
        vec3r relative_pos = (cell_src.x + offset_src - cell_tar.x);
        rtfmm::real delta = 2 * r / (P - 1);
        rtfmm::real gsize = r * 2;
        std::vector<rtfmm::vec3r> grid = get_conv_grid(N, gsize, delta, relative_pos);
        std::vector<rtfmm::real>& G = G_all[i];
        G = get_G_matrix(grid, N);
    }
    TIME_END(setup);

    TIME_BEGIN(create_plan);
    std::vector<rtfmm::real>& Q = Q_all[0];
    std::vector<rtfmm::real>& G = G_all[0];
    std::vector<rtfmm::real>& p_fft_grid = p_all[0];
    std::vector<complexr>& Qk = Qk_all[0];
    std::vector<complexr>& Gk = Gk_all[0];
    std::vector<complexr>& pk = pk_all[0];
    fftw_plan plan_G = fftw_plan_dft_r2c_3d(N, N, N, G.data(), reinterpret_cast<fftw_complex*>(Gk.data()), FFTW_ESTIMATE);
    fftw_plan plan_Q = fftw_plan_dft_r2c_3d(N, N, N, Q.data(), reinterpret_cast<fftw_complex*>(Qk.data()), FFTW_ESTIMATE);
    fftw_plan plan_P = fftw_plan_dft_c2r_3d(N, N, N, reinterpret_cast<fftw_complex*>(pk.data()), p_fft_grid.data(), FFTW_ESTIMATE);
    TIME_END(create_plan);

    TIME_BEGIN(G);
    #pragma omp parallel for
    for(int i = 0; i < num_fft; i++)
    {
        std::vector<rtfmm::real>& G = G_all[i];
        std::vector<complexr>& Gk = Gk_all[i];
        fftw_execute_dft_r2c(plan_G, G.data(), reinterpret_cast<fftw_complex*>(Gk.data()));
    }
    TIME_END(G);
    TIME_BEGIN(Q);
    #pragma omp parallel for
    for(int i = 0; i < num_src; i++)
    {
        std::vector<rtfmm::real>& Q = Q_all[i];
        std::vector<complexr>& Qk = Qk_all[i];
        fftw_execute_dft_r2c(plan_Q, Q.data(), reinterpret_cast<fftw_complex*>(Qk.data()));
    }
    TIME_END(Q);

    TIME_BEGIN(P);
    #pragma omp parallel for
    for(int i = 0; i < num_fft; i++)
    {
        std::vector<complexr>& Gk = Gk_all[i];
        std::vector<complexr>& Qk = Qk_all[m2l_src_map[m2l_pairs[i].second.first]];
        std::vector<complexr>& pk = pk_all[i];
        std::vector<rtfmm::real>& p_fft_grid = p_all[i];

        for(int i = 0; i < N * N * N; i++)
        {
            rtfmm::real a = Gk[i].r;
            rtfmm::real b = Gk[i].i;
            rtfmm::real c = Qk[i].r;
            rtfmm::real d = Qk[i].i;
            pk[i].r = a * c - b * d;
            pk[i].i = a * d + b * c;
        }
        fftw_execute_dft_c2r(plan_P, reinterpret_cast<fftw_complex*>(pk.data()), p_fft_grid.data());
    }
    TIME_END(P);

    TIME_BEGIN(add);
    for(int i = 0; i < num_fft; i++)
    {
        Cell3& cell_tar = cs[m2l_pairs[i].first];
        std::vector<rtfmm::real>& p_fft_grid = p_all[i];
        for(int i = 0; i < surface_number; i++)
        {
            rtfmm::vec3i idx3 = surf_conv_map[i] + rtfmm::vec3i(P-1,P-1,P-1);
            cell_tar.p_check.d[i] += p_fft_grid[idx3[2] * N * N + idx3[1] * N + idx3[0]] / (N * N * N);
        }
    }
    TIME_END(add);

    TIME_BEGIN(destroy_plan);
    fftw_destroy_plan(plan_G);
    fftw_destroy_plan(plan_Q);
    fftw_destroy_plan(plan_P);
    TIME_END(destroy_plan);
}

void rtfmm::LaplaceKernel::m2l_fft_precompute_advanced2(int P, Cells3& cs, PeriodicInteractionMap& m2l_map, PeriodicInteractionPairs& m2l_pairs)
{
    #define USE_MANY
    if(verbose) printf("m2l_fft_precompute_advanced2\n");
    TIME_BEGIN(init);
    int N = 2 * P - 1;
    int N3 = N * N * N;
    int N_freq = N * N * (N / 2 + 1);
    printf("N = %d, N3 = %d, N_freq = %d\n", N, N3, N_freq);
    int surface_number = rtfmm::get_surface_point_num(P);
    std::map<int,rtfmm::vec3i> surf_conv_map = get_surface_conv_map(P);
    int num_src = m2l_srcs.size();
    int num_tar = m2l_tars.size();
    printf("num_src = %d, num_tar = %d\n", num_src, num_tar);
    std::map<int,int> m2l_src_map;
    for(int i = 0; i < num_src; i++)
    {
        m2l_src_map[m2l_srcs[i]] = i;
    }
    std::vector<real> Q_all;
    std::vector<real> p_all;
    std::vector<complexr> Qk_all;
    std::vector<complexr> pk_all;
    #ifdef USE_MANY
    int batch_size = 8;
    int batch_num = (num_src + batch_size - 1) / batch_size;
    printf("num_src = %d, batch_size = %d, batch_num = %d\n", num_src, batch_size, batch_num);
    int fft_size[3] = {N,N,N};
    Q_all.reserve(batch_num * batch_size * N3);
    Qk_all.reserve(batch_num * batch_size * N_freq);
    #else
    Q_all.reserve(num_src * N3);
    Qk_all.reserve(num_src * N_freq);
    #endif
    p_all.reserve(num_tar * N3);
    pk_all.reserve(num_tar * N_freq);
    TIME_END(init);

    TIME_BEGIN(up);
    #pragma omp parallel for
    for(int i = 0; i < num_src; i++)
    {
        Cell3& cell_src = cs[m2l_srcs[i]];
        real* Q = &Q_all[i * N3];
        std::memset(Q, 0, N3 * sizeof(real));
        get_Q_matrix(Q, cell_src.q_equiv, N, surf_conv_map);
    }
    #ifdef USE_MANY
    fftw_plan plan_Q = fftw_plan_many_dft_r2c(3, fft_size, batch_size, &Q_all[0], nullptr, 1, N3, reinterpret_cast<fftw_complex*>(&Qk_all[0]), nullptr, 1, N_freq, FFTW_ESTIMATE);
    #pragma omp parallel for
    for(int batch_idx = 0; batch_idx < batch_num; batch_idx++)
    {
        int i = batch_idx * batch_size;
        real* Q = &Q_all[i * N3];
        complexr* Qk = &Qk_all[i * N_freq];
        fftw_execute_dft_r2c(plan_Q, Q, reinterpret_cast<fftw_complex*>(&Qk[0]));
    }
    #else
    fftw_plan plan_Q = fftw_plan_dft_r2c_3d(N, N, N, &Q_all[0], reinterpret_cast<fftw_complex*>(&Qk_all[0]), FFTW_ESTIMATE);
    #pragma omp parallel for
    for(int i = 0; i < num_src; i++)
    {
        real* Q = &Q_all[i * N3];
        complexr* Qk = &Qk_all[i * N_freq];
        fftw_execute_dft_r2c(plan_Q, Q, reinterpret_cast<fftw_complex*>(Qk));
    }
    #endif
    fftw_destroy_plan(plan_Q);
    TIME_END(up);

    TIME_BEGIN(hadamard);
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < num_tar; i++)
    {
        auto m2l = m2l_map.begin();
        std::advance(m2l, i);
        Cell3& cell_tar = cs[m2l->first];
        complexr* pk = &pk_all[i * N_freq];
        std::memset(pk, 0, N_freq * sizeof(complexr));
        int m2l_len = m2l->second.size();
        for(int j = 0; j < m2l_len; j++)
        {
            Cell3& cell_src = cs[m2l->second[j].first];
            vec3r rv = (cell_src.x + m2l->second[j].second - cell_tar.x) / (cell_src.r * 2);
            vec3i rvi(std::round(rv[0]), std::round(rv[1]), std::round(rv[2]));
            std::vector<complexr>& Gk = m2l_Gks[hash(rvi)];
            complexr* Qk = &Qk_all[m2l_src_map[m2l->second[j].first] * N_freq];
            for(int f = 0; f < N_freq; f++)
            {
                real a = Gk[f].r;
                real b = Gk[f].i;
                real c = Qk[f].r;
                real d = Qk[f].i;
                pk[f].r += a * c - b * d;
                pk[f].i += a * d + b * c;
            }
        }
    }
    /*std::map<int, PeriodicInteractionPairs> pair_of_matrix;
    for(auto m2l : m2l_pairs)
    {
        Cell3& cell_tar = cs[m2l.first];
        Cell3& cell_src = cs[m2l.second.first];
        vec3r offset_src = m2l.second.second;
        vec3r relative_pos = cell_src.x + offset_src - cell_tar.x;
        vec3r rv = relative_pos / (cell_src.r * 2);
        vec3i rvi(std::round(rv[0]), std::round(rv[1]), std::round(rv[2]));
        pair_of_matrix[hash(rvi)].push_back(m2l);
    }
    std::cout<<"m2l_Gks.size() = "<<m2l_Gks.size()<<std::endl;
    std::cout<<"pair_of_matrix.size() = "<<pair_of_matrix.size()<<std::endl;
    #pragma omp parallel for
    for(int i = 0; i < m2l_Gks.size(); i++)
    {
        auto m2l_Gk = m2l_Gks.begin();
        auto m2l_pairs = pair_of_matrix.begin();
        std::advance(m2l_Gk, i);
        std::advance(m2l_pairs, i);
        std::vector<complexr>& Gk = m2l_Gk->second;
        for(int j = 0; j < m2l_pairs->second.size(); j++)
        {
            auto m2l_pair = m2l_pairs->second.at(j);
            complexr* pk = &pk_all[i * N_freq];
            complexr* Qk = &Qk_all[m2l_src_map[m2l_pair.second.first] * N_freq];
            for(int f = 0; f < N_freq; f++)
            {
                rtfmm::real a = Gk[f].r;
                rtfmm::real b = Gk[f].i;
                rtfmm::real c = Qk[f].r;
                rtfmm::real d = Qk[f].i;
                pk[f].r += a * c - b * d;
                pk[f].i += a * d + b * c;
            }
        }
    }*/
    TIME_END(hadamard);

    TIME_BEGIN(down);
    fftw_plan plan_P = fftw_plan_dft_c2r_3d(N, N, N, reinterpret_cast<fftw_complex*>(&pk_all[0]), &p_all[0], FFTW_ESTIMATE);
    #pragma omp parallel for
    for(int i = 0; i < num_tar; i++)
    {
        complexr* pk = &pk_all[i * N_freq];
        real* p_fft_grid = &p_all[i * N3];
        fftw_execute_dft_c2r(plan_P, reinterpret_cast<fftw_complex*>(pk), p_fft_grid);
        auto m2l = m2l_map.begin();
        std::advance(m2l, i);
        Cell3& cell_tar = cs[m2l->first];
        real scale = std::pow(cell_tar.depth > 0 ? 2 : 3, cell_tar.depth);
        for(int i = 0; i < surface_number; i++)
        {
            rtfmm::vec3i idx3 = surf_conv_map[i] + rtfmm::vec3i(P-1,P-1,P-1);
            cell_tar.p_check.d[i] += p_fft_grid[idx3[2] * N * N + idx3[1] * N + idx3[0]] / N3 * scale;
        }
    }
    fftw_destroy_plan(plan_P);
    TIME_END(down);
}


void rtfmm::LaplaceKernel::l2l(int P, Cell3& cell_parent, Cells3& cs)
{
    /* parent check to parent equiv */
    std::vector<rtfmm::vec3r> x_equiv_parent = get_surface_points(P, cell_parent.r * 2.95, cell_parent.x);
    std::vector<rtfmm::vec3r> x_check_parent = get_surface_points(P, cell_parent.r * 1.05, cell_parent.x);
    Matrix pe2pc = get_p2p_matrix(x_equiv_parent, x_check_parent);
    Matrix q_equiv_parent = linear_equation_system_svd(pe2pc, cell_parent.p_check);

    for(int i = 0; i < cell_parent.crange.number; i++)
    {
        Cell3& cell_child = cs[cell_parent.crange.offset + i];

        /* parent equiv to child check */
        std::vector<rtfmm::vec3r> x_check_child = get_surface_points(P, cell_child.r * 1.05, cell_child.x);
        Matrix pe2cc = get_p2p_matrix(x_equiv_parent, x_check_child);
        Matrix p_check_child = mat_vec_mul(pe2cc, q_equiv_parent);
        cell_child.p_check = mat_mat_add(cell_child.p_check, p_check_child);
    }
}

void rtfmm::LaplaceKernel::l2l_precompute(int P, Cell3& cell_parent, Cells3& cs)
{
    for(int i = 0; i < cell_parent.crange.number; i++)
    {
        Cell3& cell_child = cs[cell_parent.crange.offset + i];
        Matrix p_check_child(get_surface_point_num(P), 1);
        mat_vec_mul(matrix_l2l[cell_child.octant], cell_parent.p_check, p_check_child);
        mat_mat_increment(cell_child.p_check, p_check_child);
    }
}

void rtfmm::LaplaceKernel::l2l_img_precompute(int P, Cell3& cell_parent, Cells3& cs)
{
    Cell3& cell_child = cs[cell_parent.crange.offset];
    Matrix p_check_child = mat_vec_mul(matrix_l2l_img, cell_parent.p_check);
    cell_child.p_check = mat_mat_add(cell_child.p_check, p_check_child);
}

void rtfmm::LaplaceKernel::l2p(int P, Bodies3& bs_tar, Cell3& cell_tar)
{
    /* get equivalent to check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_tar.r * 2.95, cell_tar.x);
    std::vector<rtfmm::vec3r> x_check = get_surface_points(P, cell_tar.r * 1.05, cell_tar.x);
    Matrix e2c = get_p2p_matrix(x_equiv, x_check);

    /* get equivalent charge */
    Matrix q_equiv = linear_equation_system_svd(e2c, cell_tar.p_check);

    /* get equiv to tar matrix */
    std::vector<vec3r> x_tar = get_bodies_x(bs_tar, cell_tar.brange);
    Matrix e2t = get_p2p_matrix(x_equiv, x_tar);

    /* get target potential and force */
    Matrix p_tar = mat_vec_mul(e2t, q_equiv);
    Matriv f_tar = get_force_naive(x_equiv, x_tar, q_equiv);
    add_boides_p(bs_tar, p_tar, cell_tar.brange);
    add_boides_f(bs_tar, f_tar, cell_tar.brange);
}

void rtfmm::LaplaceKernel::l2p_precompute(int P, Bodies3& bs_tar, Cell3& cell_tar)
{
    int num_surf = get_surface_point_num(P);
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_tar.r * 2.95, cell_tar.x);

    real scale = std::pow(0.5, cell_tar.depth);
    Matrix UTb = mat_vec_mul(UT_l2p_precompute, cell_tar.p_check);
    Matrix q_equiv = mat_vec_mul(VSinv_l2p_precompute, UTb, scale);

    for(int j = 0; j < cell_tar.brange.number; j++)
    {
        Body3& bj = bs_tar[cell_tar.brange.offset + j];
        real potential = 0;
        vec3r force(0,0,0);
        for(int i = 0; i < num_surf; i++)
        {
            vec3r& xi = x_equiv[i];
            real qi = q_equiv[i];
            vec3r dx = bj.x - xi;
            real r = dx.r();
            real invr = r == 0 ? 0 : 1 / r;
            potential += qi * invr;
            force += qi * invr * invr * invr * (-dx);
        }
        bj.p += potential;
        bj.f += force;
    }
}

void rtfmm::LaplaceKernel::m2p(int P, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar, vec3r offset)
{
    std::vector<rtfmm::vec3r> x_equiv_src = get_surface_points(P, cell_src.r * 1.05, cell_src.x + offset);
    std::vector<vec3r> x_tar = get_bodies_x(bs_tar, cell_tar.brange);
    Matrix e2t = get_p2p_matrix(x_equiv_src, x_tar);
    Matrix p_tar = mat_vec_mul(e2t, cell_src.q_equiv);
    Matriv f_tar = get_force_naive(x_equiv_src, x_tar, cell_src.q_equiv);
    add_boides_p(bs_tar, p_tar, cell_tar.brange);
    add_boides_f(bs_tar, f_tar, cell_tar.brange);
}

void rtfmm::LaplaceKernel::p2l(int P, Bodies3& bs_src, Cell3& cell_src, Cell3& cell_tar, vec3r offset)
{
    std::vector<vec3r> x_src = get_bodies_x(bs_src, cell_src.brange, offset);
    Matrix q_src = get_bodies_q(bs_src, cell_src.brange);
    std::vector<rtfmm::vec3r> x_check_tar = get_surface_points(P, cell_tar.r * 1.05, cell_tar.x);
    Matrix s2c = get_p2p_matrix(x_src, x_check_tar);
    Matrix p_check = mat_vec_mul(s2c, q_src);
    cell_tar.p_check = mat_mat_add(cell_tar.p_check, p_check);
}


rtfmm::Matrix rtfmm::LaplaceKernel::get_p2p_matrix(
    std::vector<vec3r>& x_src, 
    std::vector<vec3r>& x_tar
)
{
    int num_src = x_src.size();
    int num_tar = x_tar.size();
    Matrix matrix_p2p(num_tar, num_src);
    for(int j = 0; j < num_tar; j++)
    {
        vec3r xtar = x_tar[j];
        for(int i = 0; i < num_src; i++)
        {
            vec3r xsrc = x_src[i];
            vec3r dx = xtar - xsrc;
            real r = dx.r();
            real invr = r == 0 ? 0 : 1 / r;
            matrix_p2p[j * num_src + i] = invr;
        }
    }
    return matrix_p2p;
}

rtfmm::Matriv rtfmm::LaplaceKernel::get_force_naive(
    std::vector<vec3r>& x_src, 
    std::vector<vec3r>& x_tar, 
    Matrix& q_src
)
{
    Matriv res(x_tar.size(), 1);
    for(int j = 0; j < x_tar.size(); j++)
    {   
        vec3r force(0,0,0);
        for(int i = 0; i < x_src.size(); i++)
        {
            vec3r dx = x_tar[j] - x_src[i];
            real r = dx.r();
            real invr = r == 0 ? 0 : 1 / r;
            force += q_src[i] * invr * invr * invr * (-dx);
        }
        res[j] = force;
    }
    return res;
}

void rtfmm::dipole_correction(Bodies3& bs, real cycle)
{
    if(verbose) std::cout<<"dipole correction"<<std::endl;
    int num_body = bs.size();
    real coef = 4 * M_PI / (3 * cycle * cycle * cycle);
    vec3r dipole(0,0,0);
    for (int i = 0; i < num_body; i++) 
    {
        dipole += bs[i].x * bs[i].q;
    }
    real dnorm = dipole.norm();
    for (int i = 0; i < num_body; i++) 
    { 
        bs[i].p -= coef * dnorm / num_body / bs[i].q; 
        bs[i].f -= coef * dipole;
    }
}

std::vector<rtfmm::real> rtfmm::LaplaceKernel::get_G_matrix(std::vector<rtfmm::vec3r>& grid, int N)
{
    std::vector<rtfmm::real> G;
    for(int k = 0; k < N; k++)
    {
        for(int j = 0; j < N; j++)
        {
            for(int i = 0; i < N; i++)
            {
                rtfmm::real r = grid[k * N * N + j * N + i].r();
                rtfmm::real invr = r == 0 ? 0 : 1 / r;
                G.push_back(invr);
            }
        }
    }
    return G;
}

std::vector<rtfmm::real> rtfmm::LaplaceKernel::get_Q_matrix(Matrix& surface_q, int N, std::map<int,rtfmm::vec3i>& surf_conv_map)
{
    std::vector<rtfmm::real> Q(N * N * N);
    int surface_number = surface_q.m * surface_q.n;
    for(int i = 0; i < surface_number; i++)
    {
        rtfmm::vec3i idx = surf_conv_map[i];
        Q[idx[2] * N * N + idx[1] * N + idx[0]] = surface_q[i];
    }
    return Q;
}

void rtfmm::LaplaceKernel::get_Q_matrix(real* Q, Matrix& surface_q, int N, std::map<int,rtfmm::vec3i>& surf_conv_map)
{
    int surface_number = surface_q.m * surface_q.n;
    for(int i = 0; i < surface_number; i++)
    {
        rtfmm::vec3i idx = surf_conv_map[i];
        Q[idx[2] * N * N + idx[1] * N + idx[0]] = surface_q[i];
    }
}

std::map<int,rtfmm::vec3i> rtfmm::LaplaceKernel::get_surface_conv_map(int p)
{
    int num = rtfmm::get_surface_point_num(p);
    std::map<int,rtfmm::vec3i> map;
    int cnt = 0;
    for(int i = 0; i < p; i++)
    {
        for(int j = 0; j < p; j++)
        {
            for(int k = 0; k < p; k++)
            {
                if(i == 0 || i == p - 1 || j == 0 || j == p - 1 || k == 0 || k == p - 1)
                {
                    map[cnt] = rtfmm::vec3i(i,j,k);
                    cnt++;
                }
            }
        }
    }
    return map;
}

void rtfmm::LaplaceKernel::precompute(int P, real r0, int images)
{
    if(verbose) printf("precompute\n");

    /* p2m */
    std::vector<vec3r> x_check_up = get_surface_points(P, r0 * 2.95);
    std::vector<vec3r> x_equiv_up = get_surface_points(P, r0 * 1.05);
    Matrix e2c_up_precompute = get_p2p_matrix(x_equiv_up, x_check_up);
    Matrix U, S, VT;
    svd(e2c_up_precompute, U, S, VT);
    UT_p2m_precompute = transpose(U);
    V_p2m_precompute = transpose(VT);
    Sinv_p2m_precompute = pseudo_inverse(S);
    VSinv_p2m_precompute = mat_mat_mul(V_p2m_precompute, Sinv_p2m_precompute);

    /* m2m */
    std::vector<vec3r>& x_equiv_up_parent = x_equiv_up;
    std::vector<vec3r>& x_check_up_parent = x_check_up;
    Matrix& UT_m2m_precompute = UT_p2m_precompute;
    Matrix& V_m2m_precompute = V_p2m_precompute;
    Matrix& Sinv_m2m_precompute = Sinv_p2m_precompute;
    Matrix& VSinv_m2m_precompute = VSinv_p2m_precompute;
    #pragma omp parallel for
    for(int octant = 0; octant < 8; octant++)
    {
        vec3r offset_child = Tree::get_child_cell_x(vec3r(0,0,0), r0, octant, 0);
        std::vector<vec3r> x_equiv_child_up = get_surface_points(P, r0 / 2 * 1.05, offset_child);
        Matrix ce2pc_up_precompute = get_p2p_matrix(x_equiv_child_up, x_check_up_parent);
        Matrix m1 = mat_mat_mul(UT_m2m_precompute, ce2pc_up_precompute);
        matrix_m2m[octant] = mat_mat_mul(VSinv_m2m_precompute, m1); // (VSinv)(UTG)
    }

    /* m2m image */
    if(images > 0)
    {
        #pragma omp parallel for
        for(int octant = 0; octant < 27; octant++)
        {
            vec3r offset_child = Tree::get_child_cell_x(vec3r(0,0,0), r0, octant, 1);
            std::vector<vec3r> x_equiv_child_up = get_surface_points(P, r0 / 3 * 1.05, offset_child);
            Matrix ce2pc_up_precompute = get_p2p_matrix(x_equiv_child_up, x_check_up_parent);

            Matrix m1 = mat_mat_mul(UT_m2m_precompute, ce2pc_up_precompute);
            matrix_m2m_img[octant] = mat_mat_mul(VSinv_m2m_precompute, m1); // (VSinv)(UTG)
        }
    }

    /* l2l */
    #pragma omp parallel for
    for(int octant = 0; octant < 8; octant++)
    {
        matrix_l2l[octant] = transpose(matrix_m2m[octant]); // (G(VSinv))UT
    }

    /* l2l image */
    if(images > 0)
    {
        matrix_l2l_img = transpose(matrix_m2m_img[13]); // (G(VSinv))UT
    }

    /* l2p */
    UT_l2p_precompute = transpose(V_p2m_precompute);
    V_l2p_precompute = transpose(UT_p2m_precompute);
    Sinv_l2p_precompute = Sinv_p2m_precompute;
    VSinv_l2p_precompute = mat_mat_mul(V_l2p_precompute, Sinv_l2p_precompute);
    VSinvUT_l2p_precompute = mat_mat_mul(VSinv_l2p_precompute, UT_l2p_precompute);
}

void rtfmm::LaplaceKernel::precompute_m2l(int P, real r0, Cells3 cs, PeriodicInteractionMap m2l_map, int images)
{
    std::set<int> m2l_tar_set, m2l_src_set;
    for(auto m2l_list : m2l_map)
    {
        m2l_tar_set.insert(m2l_list.first);
        for(int i = 0; i < m2l_list.second.size(); i++)
        {
            m2l_src_set.insert(m2l_list.second[i].first);
        }
    }
    m2l_tars = std::vector<int>(m2l_tar_set.begin(), m2l_tar_set.end());
    m2l_srcs = std::vector<int>(m2l_src_set.begin(), m2l_src_set.end());

    if(verbose)
    {
        printf("# of cells = %ld\n", cs.size());
        printf("# of m2l_tar = %ld\n", m2l_tars.size());
        printf("# of m2l_src = %ld\n", m2l_srcs.size());
    }

    for(int i = 0; i < m2l_tars.size(); i++)
    {
        assert_exit(m2l_tars[i] == m2l_srcs[i], "m2l inconsisent");
    }

    int N = 2 * P - 1;
    rtfmm::real delta = 2 * r0 / (P - 1) * 1.05;
    rtfmm::real gsize = r0 * 2 * 1.05;
    int N3 = N * N * N;
    int N_freq = N * N * (N / 2 + 1);
    printf("N_freq = %d\n", N_freq);
    std::vector<real> G0(N3);
    std::vector<complexr> Gk0(N_freq);
    fftw_plan plan_G = fftw_plan_dft_r2c_3d(N, N, N, G0.data(), reinterpret_cast<fftw_complex*>(Gk0.data()), FFTW_ESTIMATE);
    int range = images >= 2 ? 4 : 3;
    #pragma omp parallel for collapse(3)
    for(int k = -range; k <= range; k++)
    {
        for(int j = -range; j <= range; j++)
        {
            for(int i = -range; i <= range; i++)
            {
                if(std::abs(i) > 1 || std::abs(j) > 1 || std::abs(k) > 1)
                {
                    vec3i relative_idx(i,j,k);
                    vec3r relative_pos = vec3r(i,j,k) * r0 * 2;
                    std::vector<rtfmm::vec3r> grid = get_conv_grid(N, gsize, delta, relative_pos);
                    std::vector<real> G = get_G_matrix(grid, N); 
                    std::vector<complexr> Gk(N_freq);                   
                    fftw_execute_dft_r2c(plan_G, G.data(), reinterpret_cast<fftw_complex*>(Gk.data()));
                    #pragma omp critical // std::map is not thread-safe
                    m2l_Gks[hash(relative_idx)] = Gk;
                }
            }
        }
    }
    fftw_destroy_plan(plan_G);
}