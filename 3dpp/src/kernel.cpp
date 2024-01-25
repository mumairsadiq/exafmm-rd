#include "kernel.h"
#include "mathfunc.h"
#include "surface.h"
#include "argument.h"
#include <omp.h>
#include <set>
#include "align.h"
#include <algorithm>

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

rtfmm::LaplaceKernel::LaplaceKernel(){}

/*void rtfmm::LaplaceKernel::p2p_p(std::vector<vec3r>& xs_src, std::vector<real>& qs_src, std::vector<vec3r>& xs_tar, std::vector<real>& ps_tar, const vec3r& offset)
{

}*/

void rtfmm::LaplaceKernel::p2p_pf(std::vector<vec3r>& xs_src, std::vector<real>& qs_src, std::vector<vec3r>& xs_tar, std::vector<real>& ps_tar, std::vector<vec3r>& fs_tar, const vec3r& offset)
{
    int num_tar = xs_tar.size();
    int num_src = xs_src.size();
    vec4d zero(0),one(1), xj, yj, zj, dx, dy, dz, qi, invr, f, qinvr;
    for(int j = 0; j < num_tar; j+=4)
    {
        vec3r* xtar = &xs_tar[j];
        int padding = j + 4 - num_tar;
        padding = padding <= 0 ? 0 : padding;
        if(padding <= 0)
        {
            xj = vec4d((*(xtar+0))[0],(*(xtar+1))[0],(*(xtar+2))[0],(*(xtar+3))[0]);
            yj = vec4d((*(xtar+0))[1],(*(xtar+1))[1],(*(xtar+2))[1],(*(xtar+3))[1]);
            zj = vec4d((*(xtar+0))[2],(*(xtar+1))[2],(*(xtar+2))[2],(*(xtar+3))[2]);
        }
        else if(padding == 1)
        {
            xj = vec4d((*(xtar+0))[0],(*(xtar+1))[0],(*(xtar+2))[0],0);
            yj = vec4d((*(xtar+0))[1],(*(xtar+1))[1],(*(xtar+2))[1],0);
            zj = vec4d((*(xtar+0))[2],(*(xtar+1))[2],(*(xtar+2))[2],0);
        }
        else if(padding == 2)
        {
            xj = vec4d((*(xtar+0))[0],(*(xtar+1))[0],0,0);
            yj = vec4d((*(xtar+0))[1],(*(xtar+1))[1],0,0);
            zj = vec4d((*(xtar+0))[2],(*(xtar+1))[2],0,0);
        }
        else if(padding == 3)
        {
            xj = vec4d((*(xtar+0))[0],0,0,0);
            yj = vec4d((*(xtar+0))[1],0,0,0);
            zj = vec4d((*(xtar+0))[2],0,0,0);
        }
        vec4d pj(0),fxj(0),fyj(0),fzj(0);
        for(int i = 0; i < num_src; i++)
        {
            vec3r bsrc_xi = xs_src[i] + offset;
            dx = vec4d(bsrc_xi[0]);
            dy = vec4d(bsrc_xi[1]);
            dz = vec4d(bsrc_xi[2]);
            qi = vec4d(qs_src[i]);
            dx = xj - dx;
            dy = yj - dy;
            dz = zj - dz;
            invr = dx * dx + dy * dy + dz * dz;
            f = invr > zero;
            invr = one / invr.sqrt();
            invr &= f;
            qinvr = qi * invr;
            pj += qinvr;
            qinvr = qinvr * invr * invr;
            fxj += qinvr * dx;
            fyj += qinvr * dy;
            fzj += qinvr * dz;
        }
        for(int k = 0; k < 4 - padding; k++)
        {
            ps_tar[j + k] += pj[k];
            fs_tar[j + k][0] -= fxj[k];
            fs_tar[j + k][1] -= fyj[k];
            fs_tar[j + k][2] -= fzj[k];
        }
    }
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
        vec4d zero(0);
        vec4d one(1);
        for(int j = 0; j < cell_tar.brange.number; j+=4)
        {
            Body3* btar = &bs_tar[cell_tar.brange.offset + j];
            vec4d xj,yj,zj;
            int padding = j + 4 - cell_tar.brange.number;
            padding = padding <= 0 ? 0 : padding;
            if(padding <= 0)
            {
                xj = vec4d(btar->x[0],(btar+1)->x[0],(btar+2)->x[0],(btar+3)->x[0]);
                yj = vec4d(btar->x[1],(btar+1)->x[1],(btar+2)->x[1],(btar+3)->x[1]);
                zj = vec4d(btar->x[2],(btar+1)->x[2],(btar+2)->x[2],(btar+3)->x[2]);
            }
            else if(padding == 1)
            {
                xj = vec4d(btar->x[0],(btar+1)->x[0],(btar+2)->x[0],0);
                yj = vec4d(btar->x[1],(btar+1)->x[1],(btar+2)->x[1],0);
                zj = vec4d(btar->x[2],(btar+1)->x[2],(btar+2)->x[2],0);
            }
            else if(padding == 2)
            {
                xj = vec4d(btar->x[0],(btar+1)->x[0],0,0);
                yj = vec4d(btar->x[1],(btar+1)->x[1],0,0);
                zj = vec4d(btar->x[2],(btar+1)->x[2],0,0);
            }
            else if(padding == 3)
            {
                xj = vec4d(btar->x[0],0,0,0);
                yj = vec4d(btar->x[1],0,0,0);
                zj = vec4d(btar->x[2],0,0,0);
            }
            
            vec4d pj(0);
            vec4d fxj(0);
            vec4d fyj(0);
            vec4d fzj(0);

            vec4d dx, dy, dz, qi, invr, f, qinvr;
            for(int i = 0; i < cell_src.brange.number; i++)
            {
                Body3& bsrc = bs_src[cell_src.brange.offset + i];
                vec3r bsrc_xi = bsrc.x + offset;
                dx = vec4d(bsrc_xi[0]);
                dy = vec4d(bsrc_xi[1]);
                dz = vec4d(bsrc_xi[2]);
                qi = vec4d(bsrc.q);
                dx = xj - dx;
                dy = yj - dy;
                dz = zj - dz;
                invr = dx * dx + dy * dy + dz * dz;
                f = invr > zero;
                invr = one / invr.sqrt();
                invr &= f;
                qinvr = qi * invr;
                pj += qinvr;
                qinvr = qinvr * invr * invr;
                fxj += qinvr * dx;
                fyj += qinvr * dy;
                fzj += qinvr * dz;
            }
            for(int k = 0; k < 4 - padding; k++)
            {
                (btar + k)->p    += pj[k];
                (btar + k)->f[0] -= fxj[k];
                (btar + k)->f[1] -= fyj[k];
                (btar + k)->f[2] -= fzj[k];
            }
        }
    }
}

void rtfmm::LaplaceKernel::p2p_1toN_128(Bodies3& bs_src, Bodies3& bs_tar, Cells3& cs, std::vector<std::pair<int, vec3r>>& p2ps, Cell3& cell_tar)
{
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
}

void rtfmm::LaplaceKernel::p2p_1toN_256(Bodies3& bs_src, Bodies3& bs_tar, Cells3& cs, std::vector<std::pair<int, vec3r>>& p2ps, Cell3& cell_tar)
{
    /*__m256d zero = _mm256_set1_pd(0);
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
            double pdata[4];
        };
        union{
            __m256d fxj = _mm256_set1_pd(0);
            double fxjdata[4];
        };
        union{
            __m256d fyj = _mm256_set1_pd(0);
            double fyjdata[4];
        };
        union{
            __m256d fzj = _mm256_set1_pd(0);
            double fzjdata[4];
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
    }*/

    vec4d zero(0);
    vec4d one(1);
    for(int j = 0; j < cell_tar.brange.number; j+=4)
    {
        Body3* btar = &bs_tar[cell_tar.brange.offset + j];
        vec4d xj,yj,zj;
        int padding = j + 4 - cell_tar.brange.number;
        padding = padding <= 0 ? 0 : padding;
        if(padding <= 0)
        {
            xj = vec4d(btar->x[0],(btar+1)->x[0],(btar+2)->x[0],(btar+3)->x[0]);
            yj = vec4d(btar->x[1],(btar+1)->x[1],(btar+2)->x[1],(btar+3)->x[1]);
            zj = vec4d(btar->x[2],(btar+1)->x[2],(btar+2)->x[2],(btar+3)->x[2]);
        }
        else if(padding == 1)
        {
            xj = vec4d(btar->x[0],(btar+1)->x[0],(btar+2)->x[0],0);
            yj = vec4d(btar->x[1],(btar+1)->x[1],(btar+2)->x[1],0);
            zj = vec4d(btar->x[2],(btar+1)->x[2],(btar+2)->x[2],0);
        }
        else if(padding == 2)
        {
            xj = vec4d(btar->x[0],(btar+1)->x[0],0,0);
            yj = vec4d(btar->x[1],(btar+1)->x[1],0,0);
            zj = vec4d(btar->x[2],(btar+1)->x[2],0,0);
        }
        else if(padding == 3)
        {
            xj = vec4d(btar->x[0],0,0,0);
            yj = vec4d(btar->x[1],0,0,0);
            zj = vec4d(btar->x[2],0,0,0);
        }
        
        vec4d pj(0);
        vec4d fxj(0);
        vec4d fyj(0);
        vec4d fzj(0);

        vec4d dx, dy, dz, qi, invr, f, qinvr;
        for(int si = 0; si < p2ps.size(); si++)
        {
            Cell3& cell_src = cs[p2ps[si].first];
            vec3r& offset = p2ps[si].second;
            for(int i = 0; i < cell_src.brange.number; i++)
            {
                Body3& bsrc = bs_src[cell_src.brange.offset + i];
                vec3r bsrc_xi = bsrc.x + offset;
                dx = vec4d(bsrc_xi[0]);
                dy = vec4d(bsrc_xi[1]);
                dz = vec4d(bsrc_xi[2]);
                qi = vec4d(bsrc.q);
                dx = xj - dx;
                dy = yj - dy;
                dz = zj - dz;
                invr = dx * dx + dy * dy + dz * dz;
                f = invr > zero;
                invr = one / invr.sqrt();
                invr &= f;
                qinvr = qi * invr;
                pj += qinvr;
                qinvr = qinvr * invr * invr;
                fxj += qinvr * dx;
                fyj += qinvr * dy;
                fzj += qinvr * dz;
            }
        }
        for(int k = 0; k < 4 - padding; k++)
        {
            (btar + k)->p    += pj[k];
            (btar + k)->f[0] -= fxj[k];
            (btar + k)->f[1] -= fyj[k];
            (btar + k)->f[2] -= fzj[k];
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

void rtfmm::LaplaceKernel::m2l_fft_precompute_advanced2(int P, Cells3& cs, PeriodicInteractionMap& m2l_map)
{
    if(verbose)
    {
        int interaciton_number = 0;
        for(auto m2l : m2l_map)
        {
            std::cout<<m2l.first<<","<<m2l.second.size()<<std::endl;
            interaciton_number += m2l.second.size();
        }
        std::cout<<"interaciton_number = "<<interaciton_number<<std::endl;
    }
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
    std::vector<real,AlignAllocator> Q_all;
    std::vector<real,AlignAllocator> p_all;
    std::vector<complexr,AlignAllocator> Qk_all;
    std::vector<complexr,AlignAllocator> pk_all;
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

void rtfmm::LaplaceKernel::m2l_fft_precompute_advanced3(int P, Cells3& cs, PeriodicInteractionMap& m2l_map, PeriodicInteractionPairs& m2l_pairs)
{
    if(verbose) printf("m2l_fft_precompute_advanced3\n");
    int num_src = m2l_srcs.size();
    int num_tar = m2l_tars.size();
    int num_gk = m2l_Gk_idx.size();
    int N = 2 * P - 1;
    int N3 = N * N * N;
    int N_freq = N * N * (N / 2 + 1);
    int surface_number = rtfmm::get_surface_point_num(P);
    std::map<int,rtfmm::vec3i> surf_conv_map = get_surface_conv_map(P);

    // make local pairs 
    tbegin(make_local_pairs);
    std::map<int,int> m2l_src_map;
    std::map<int,int> m2l_tar_map;
    for(int i = 0; i < num_src; i++)
    {
        m2l_src_map[m2l_srcs[i]] = i;
    }
    for(int i = 0; i < num_tar; i++)
    {
        m2l_tar_map[m2l_tars[i]] = i;
    }
    //std::vector<vec3i> tar_src_gk;
    std::vector<std::vector<vec2i>> gk_tar_src(num_gk);
    for(auto m2l : m2l_pairs)
    {
        int tar_idx = m2l.first;
        int src_idx = m2l.second.first;
        int tar_local_idx = m2l_tar_map[tar_idx];
        int src_local_idx = m2l_src_map[src_idx];
        Cell3& cell_tar = cs[tar_idx];
        Cell3& cell_src = cs[src_idx];
        vec3r offset = m2l.second.second;
        vec3r rv = (cell_src.x + offset - cell_tar.x) / (cell_src.r * 2);
        vec3i rvi(std::round(rv[0]), std::round(rv[1]), std::round(rv[2]));
        int gk_idx = m2l_Gk_idx[hash(rvi)].second;
        //tar_src_gk.push_back(vec3i(tar_local_idx, src_local_idx, gk_idx));
        gk_tar_src[gk_idx].push_back({tar_local_idx, src_local_idx});
    }
    tend(make_local_pairs);

    // make all gks
    tbegin(make_all_gks);
    std::vector<complexr> all_gks;
    all_gks.reserve(N_freq * num_gk);
    for(auto gk_idx : m2l_Gk_idx)
    {
        std::vector<complexr>& gk = gk_idx.second.first;
        int idx = gk_idx.second.second;
        for(int f = 0; f < N_freq; f++)
        {
            all_gks[f * num_gk + idx] = gk[f];
        }
    }
    tend(make_all_gks);

    // up
    tbegin(up);
    std::vector<real,AlignAllocator> Q_all;
    std::vector<complexr,AlignAllocator> Qk_all;
    Q_all.reserve(num_src * N3); 
    Qk_all.reserve(N_freq * num_src);
    #pragma omp parallel for
    for(int i = 0; i < num_src; i++)
    {
        Cell3& cell_src = cs[m2l_srcs[i]];
        real* Q = &Q_all[i * N3];
        std::memset(Q, 0, N3 * sizeof(real));
        get_Q_matrix(Q, cell_src.q_equiv, N, surf_conv_map);
    }
    fftw_plan plan_Q = fftw_plan_dft_r2c_3d(N, N, N, &Q_all[0], reinterpret_cast<fftw_complex*>(&Qk_all[0]), FFTW_ESTIMATE);
    #pragma omp parallel for
    for(int i = 0; i < num_src; i++)
    {
        real* Q = &Q_all[i * N3];
        std::vector<complexr> buffer;
        buffer.reserve(N_freq);
        fftw_execute_dft_r2c(plan_Q, Q, reinterpret_cast<fftw_complex*>(&buffer[0]));
        for(int f = 0; f < N_freq; f++)
        {
            Qk_all[f * num_src + i] = buffer[f];
        }
    }
    fftw_destroy_plan(plan_Q);
    tend(up);

    std::vector<complexr,AlignAllocator> pk_all;
    pk_all.resize(N_freq * num_tar);

    tbegin(hadamard);
    #pragma omp parallel for
    for(int f = 0; f < N_freq; f++)
    {
        for(int idx = 0; idx < gk_tar_src.size(); idx++)
        {
            complexr* Gkf = &all_gks[f * num_gk + idx];
            real a = Gkf->r;
            real b = Gkf->i;
            std::vector<vec2i> pairs = gk_tar_src[idx];
            for(int pidx = 0; pidx < pairs.size(); pidx+=1)
            {
                complexr* Qkf1 = &Qk_all[f * num_src + pairs[pidx][1]];
                complexr* pkf1 = &pk_all[f * num_tar + pairs[pidx][0]];
                pkf1->r += a * Qkf1->r - b * Qkf1->i;
                pkf1->i += a * Qkf1->i + b * Qkf1->r;
            }
        }
    }
    tend(hadamard);

    tbegin(down);
    std::vector<real,AlignAllocator> p_all;
    p_all.reserve(num_tar * N3);
    fftw_plan plan_P = fftw_plan_dft_c2r_3d(N, N, N, reinterpret_cast<fftw_complex*>(&pk_all[0]), &p_all[0], FFTW_ESTIMATE);
    #pragma omp parallel for
    for(int i = 0; i < num_tar; i++)
    {
        std::vector<complexr> buffer;
        buffer.reserve(N_freq);
        for(int f = 0; f < N_freq; f++)
        {
            buffer[f] = pk_all[f * num_tar + i];
        }
        real* p_fft_grid = &p_all[i * N3];
        fftw_execute_dft_c2r(plan_P, reinterpret_cast<fftw_complex*>(&buffer[0]), p_fft_grid);
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
    tend(down);
}

void rtfmm::LaplaceKernel::matmult_8x8x1(real*& M_, real*& IN0, real*& OUT0)
{
    real out_reg000, out_reg001, out_reg010, out_reg011;
    real  in_reg000,  in_reg001,  in_reg010,  in_reg011;
    real   m_reg000,   m_reg001,   m_reg010,   m_reg011;
    real   m_reg100,   m_reg101,   m_reg110,   m_reg111;
    for(int i1=0;i1<8;i1+=2){
      real* IN0_=IN0;
      out_reg000=OUT0[ 0]; out_reg001=OUT0[ 1];
      out_reg010=OUT0[ 2]; out_reg011=OUT0[ 3];
      for(int i2=0;i2<8;i2+=2){
        m_reg000=M_[ 0]; m_reg001=M_[ 1];
        m_reg010=M_[ 2]; m_reg011=M_[ 3];
        m_reg100=M_[16]; m_reg101=M_[17];
        m_reg110=M_[18]; m_reg111=M_[19];

        in_reg000=IN0_[0]; in_reg001=IN0_[1];
        in_reg010=IN0_[2]; in_reg011=IN0_[3];

        out_reg000 += m_reg000*in_reg000 - m_reg001*in_reg001;
        out_reg001 += m_reg000*in_reg001 + m_reg001*in_reg000;
        out_reg010 += m_reg010*in_reg000 - m_reg011*in_reg001;
        out_reg011 += m_reg010*in_reg001 + m_reg011*in_reg000;

        out_reg000 += m_reg100*in_reg010 - m_reg101*in_reg011;
        out_reg001 += m_reg100*in_reg011 + m_reg101*in_reg010;
        out_reg010 += m_reg110*in_reg010 - m_reg111*in_reg011;
        out_reg011 += m_reg110*in_reg011 + m_reg111*in_reg010;

        M_+=32; // Jump to (column+2).
        IN0_+=4;
      }
      OUT0[ 0]=out_reg000; OUT0[ 1]=out_reg001;
      OUT0[ 2]=out_reg010; OUT0[ 3]=out_reg011;
      M_+=4-64*2; // Jump back to first column (row+2).
      OUT0+=4;
    }
}

void rtfmm::LaplaceKernel::matmult_8x8x1_naive(real*& M_, real*& IN, real*& OUT)
{
    real pkr, pki;
    real qkr, qki;
    real gkr, gki;
    for(int idx_out = 0; idx_out < 8; idx_out++)
    {
        real* IN0_ = IN;
        pkr = OUT[0]; 
        pki = OUT[1];
        for(int idx_in = 0; idx_in < 8; idx_in++)
        {
            gkr = M_[ 0]; 
            gki = M_[ 1];
            qkr = IN0_[0]; 
            qki = IN0_[1];
            pkr += gkr * qkr - gki * qki;
            pki += gkr * qki + gki * qkr; 
            M_ += 2 * 8; // move to next row
            IN0_ += 2; // move to next in complex
        }
        OUT[0] = pkr; 
        OUT[1] = pki;
        M_ += -2 * 8 * 8 + 2; // move to next colum
        OUT += 2; // move to next out complex
    }
}

void rtfmm::LaplaceKernel::matmult_8x8x2(real*& M_, real*& IN0, real*& IN1, real*& OUT0, real*& OUT1)
{
    real out_reg000, out_reg001, out_reg010, out_reg011;
    real out_reg100, out_reg101, out_reg110, out_reg111;
    real  in_reg000,  in_reg001,  in_reg010,  in_reg011;
    real  in_reg100,  in_reg101,  in_reg110,  in_reg111;
    real   m_reg000,   m_reg001,   m_reg010,   m_reg011;
    real   m_reg100,   m_reg101,   m_reg110,   m_reg111;
    for(int i1=0;i1<8;i1+=2){
      real* IN0_=IN0;
      real* IN1_=IN1;
      out_reg000=OUT0[ 0]; out_reg001=OUT0[ 1];
      out_reg010=OUT0[ 2]; out_reg011=OUT0[ 3];
      out_reg100=OUT1[ 0]; out_reg101=OUT1[ 1];
      out_reg110=OUT1[ 2]; out_reg111=OUT1[ 3];
      for(int i2=0;i2<8;i2+=2){
        m_reg000=M_[ 0]; m_reg001=M_[ 1];
        m_reg010=M_[ 2]; m_reg011=M_[ 3];
        m_reg100=M_[16]; m_reg101=M_[17];
        m_reg110=M_[18]; m_reg111=M_[19];

        in_reg000=IN0_[0]; in_reg001=IN0_[1];
        in_reg010=IN0_[2]; in_reg011=IN0_[3];
        in_reg100=IN1_[0]; in_reg101=IN1_[1];
        in_reg110=IN1_[2]; in_reg111=IN1_[3];

        out_reg000 += m_reg000*in_reg000 - m_reg001*in_reg001;
        out_reg001 += m_reg000*in_reg001 + m_reg001*in_reg000;
        out_reg010 += m_reg010*in_reg000 - m_reg011*in_reg001;
        out_reg011 += m_reg010*in_reg001 + m_reg011*in_reg000;

        out_reg000 += m_reg100*in_reg010 - m_reg101*in_reg011;
        out_reg001 += m_reg100*in_reg011 + m_reg101*in_reg010;
        out_reg010 += m_reg110*in_reg010 - m_reg111*in_reg011;
        out_reg011 += m_reg110*in_reg011 + m_reg111*in_reg010;

        out_reg100 += m_reg000*in_reg100 - m_reg001*in_reg101;
        out_reg101 += m_reg000*in_reg101 + m_reg001*in_reg100;
        out_reg110 += m_reg010*in_reg100 - m_reg011*in_reg101;
        out_reg111 += m_reg010*in_reg101 + m_reg011*in_reg100;

        out_reg100 += m_reg100*in_reg110 - m_reg101*in_reg111;
        out_reg101 += m_reg100*in_reg111 + m_reg101*in_reg110;
        out_reg110 += m_reg110*in_reg110 - m_reg111*in_reg111;
        out_reg111 += m_reg110*in_reg111 + m_reg111*in_reg110;

        M_+=32; // Jump to (column+2).
        IN0_+=4;
        IN1_+=4;
      }
      OUT0[ 0]=out_reg000; OUT0[ 1]=out_reg001;
      OUT0[ 2]=out_reg010; OUT0[ 3]=out_reg011;
      OUT1[ 0]=out_reg100; OUT1[ 1]=out_reg101;
      OUT1[ 2]=out_reg110; OUT1[ 3]=out_reg111;
      M_+=4-64*2; // Jump back to first column (row+2).
      OUT0+=4;
      OUT1+=4;
    }
  }

void rtfmm::LaplaceKernel::matmult_8x8x2_avx(double*& M_, double*& IN0, double*& IN1, double*& OUT0, double*& OUT1)
{
    __m256d out00,out01,out10,out11;
    __m256d out20,out21,out30,out31;
    double* in0__ = IN0;
    double* in1__ = IN1;
    out00 = _mm256_load_pd(OUT0);
    out01 = _mm256_load_pd(OUT1);
    out10 = _mm256_load_pd(OUT0+4);
    out11 = _mm256_load_pd(OUT1+4);
    out20 = _mm256_load_pd(OUT0+8);
    out21 = _mm256_load_pd(OUT1+8);
    out30 = _mm256_load_pd(OUT0+12);
    out31 = _mm256_load_pd(OUT1+12);
    for(int i2=0;i2<8;i2+=2){
      __m256d m00;
      __m256d ot00;
      __m256d mt0,mtt0;
      __m256d in00,in00_r,in01,in01_r;
      in00 = _mm256_broadcast_pd((const __m128d*)in0__);
      in00_r = _mm256_permute_pd(in00,5);
      in01 = _mm256_broadcast_pd((const __m128d*)in1__);
      in01_r = _mm256_permute_pd(in01,5);
      m00 = _mm256_load_pd(M_);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out00 = _mm256_add_pd(out00,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out01 = _mm256_add_pd(out01,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+4);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out10 = _mm256_add_pd(out10,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out11 = _mm256_add_pd(out11,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+8);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out20 = _mm256_add_pd(out20,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out21 = _mm256_add_pd(out21,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+12);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out30 = _mm256_add_pd(out30,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out31 = _mm256_add_pd(out31,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      in00 = _mm256_broadcast_pd((const __m128d*) (in0__+2));
      in00_r = _mm256_permute_pd(in00,5);
      in01 = _mm256_broadcast_pd((const __m128d*) (in1__+2));
      in01_r = _mm256_permute_pd(in01,5);
      m00 = _mm256_load_pd(M_+16);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out00 = _mm256_add_pd(out00,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out01 = _mm256_add_pd(out01,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+20);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out10 = _mm256_add_pd(out10,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out11 = _mm256_add_pd(out11,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+24);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out20 = _mm256_add_pd(out20,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out21 = _mm256_add_pd(out21,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+28);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out30 = _mm256_add_pd(out30,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out31 = _mm256_add_pd(out31,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      M_ += 32;
      in0__ += 4;
      in1__ += 4;
    }
    _mm256_store_pd(OUT0,out00);
    _mm256_store_pd(OUT1,out01);
    _mm256_store_pd(OUT0+4,out10);
    _mm256_store_pd(OUT1+4,out11);
    _mm256_store_pd(OUT0+8,out20);
    _mm256_store_pd(OUT1+8,out21);
    _mm256_store_pd(OUT0+12,out30);
    _mm256_store_pd(OUT1+12,out31);
  }

void rtfmm::LaplaceKernel::m2l_fft_precompute_t(int P, Cells3& cs, PeriodicInteractionMap& m2l_parent_map)
{
    std::cout<<"m2l_fft_precompute_t\n";
    std::cout<<"m2l_parent_map.size() = "<<m2l_parent_map.size()<<std::endl;
    const int NCHILD = 8;
    int num_cell = cs.size();
    int num_surf = get_surface_point_num(P);
    int N = 2 * P - 1;
    int N3 = N * N * N;
    int N_freq = N * N * (N / 2 + 1);
    std::map<int,rtfmm::vec3i> surf_conv_map = get_surface_conv_map(P);
    size_t fft_size = 2 * NCHILD * N_freq;
    typedef std::vector<real, AlignAllocator> AlignedVec;

    // setup
    tbegin(setup);
    std::map<int, int> rel2idx_map;
    int rel2idx_map_idx = 0;
    for(int z = -1; z <= 1; z++)
    {
        for(int y = -1; y <= 1; y++)
        {
            for(int x = -1; x <= 1; x++)
            {
                if(std::abs(x) >= 1 || std::abs(y) >= 1 || std::abs(z) >= 1)
                {
                    rel2idx_map[hash(vec3i(x,y,z))] = rel2idx_map_idx;
                    rel2idx_map_idx++;
                }
            }
        }
    }

    std::vector<int> tar_cell_idxs;
    std::set<int> src_cell_idxs_;
    for(auto m2l_parent : m2l_parent_map)
    {
        tar_cell_idxs.push_back(m2l_parent.first);
        //std::cout<<"tar = "<<m2l_parent.first<<std::endl;
        for(int i = 0; i < m2l_parent.second.size(); i++)
        {
            src_cell_idxs_.insert(m2l_parent.second[i].first);
        }
    }
    std::vector<int> src_cell_idxs;
    for(auto it = src_cell_idxs_.begin(); it != src_cell_idxs_.end(); it++)
    {
        src_cell_idxs.push_back(*it);
    }
    std::map<int, int> src_cell_local_map;  // convert src_cell_idx in all cells to local index
    for(int i = 0; i < src_cell_idxs.size(); i++)
    {
        int src_idx = src_cell_idxs[i];
        src_cell_local_map[src_idx] = i;
    }

    printf("m2l_setup  src_cell_idxs.size() = %ld, tar_cell_idxs.size() = %ld\n", src_cell_idxs.size(),tar_cell_idxs.size());
    std::vector<size_t> fft_offset(src_cell_idxs.size());
    std::vector<size_t> ifft_offset(tar_cell_idxs.size());
    for(int i = 0; i < src_cell_idxs.size(); i++)
    {
        fft_offset[i] = cs[src_cell_idxs[i]].crange.offset * num_surf;
    }
    for(int i = 0; i < tar_cell_idxs.size(); i++)
    {
        ifft_offset[i] = cs[tar_cell_idxs[i]].crange.offset * num_surf;
    }
    std::map<int,std::map<int,int>> m2l_parent_octant_map;//tar_idx -> [octant->src_idx]
    for(int i = 0; i < tar_cell_idxs.size(); i++)
    {
        int tar_idx = tar_cell_idxs[i];
        Cell3& tar_cell = cs[tar_idx];
        std::vector<std::pair<int, rtfmm::vec3r>> m2l_list = m2l_parent_map[tar_idx];
        for(int d = 0; d < 26; d++)
        {
            m2l_parent_octant_map[tar_idx][d] = -1;
        }
        for(int j = 0; j < m2l_list.size(); j++)
        {
            int src_idx = m2l_list[j].first;
            Cell3& src_cell = cs[src_idx];
            vec3r rv = (src_cell.x - tar_cell.x) / (tar_cell.r * 2);
            vec3i rvi(std::round(rv[0]), std::round(rv[1]), std::round(rv[2]));
            int octant = rel2idx_map[hash(rvi)];
            if(rel2idx_map.find(hash(rvi)) == rel2idx_map.end())
            {
                std::cout<<"cannot find rvi\n";
                exit(1);
            }
            if(std::abs(rvi[0]) > 1 || std::abs(rvi[1]) > 1 || std::abs(rvi[1]) > 1 || octant >= 26)
            {
                std::cout<<"error!"<<std::endl;
                exit(1);
            }
            m2l_parent_octant_map[tar_idx][octant] = src_idx;
        }
    }
    std::vector<size_t> interaction_offset_f;
    std::vector<size_t> interaction_count_offset;
    size_t interaction_count_offset_ = 0;
    for(int k = 0; k < 26; k++)
    {
        for(int i = 0; i < tar_cell_idxs.size(); i++)
        {
            int tar_idx = tar_cell_idxs[i];
            int src_idx = m2l_parent_octant_map[tar_idx][k];
            if(src_idx != -1)
            {
                if(src_cell_local_map.find(src_idx) == src_cell_local_map.end())
                {
                    std::cout<<"cannnot find srcidx\n";
                    exit(1);
                }
                int src_cell_local_idx = src_cell_local_map[src_idx];
                if(src_cell_local_idx < 0 || src_cell_local_idx > src_cell_idxs.size() - 1)
                {
                    std::cout<<"src_cell_local_idx error\n";
                    exit(1);
                }
                //printf("[%d] %d(%d) -> %d(%d)\n", k, src_cell_local_idx, src_idx, i, tar_idx);
                interaction_offset_f.push_back(src_cell_local_idx * fft_size);
                interaction_offset_f.push_back(i * fft_size);
                interaction_count_offset_++;
            }
        }
        interaction_count_offset.push_back(interaction_count_offset_);
        //std::cout<<interaction_count_offset_<<std::endl;
    }
    std::cout<<"interaction_count_offset.size = "<<interaction_count_offset.size()<<std::endl;
    tend(setup);


    // allocate memory
    std::vector<real> all_up_equiv, all_dn_equiv;
    all_up_equiv.reserve(num_cell * num_surf);   // use reserve() to avoid the overhead of calling constructor
    all_dn_equiv.reserve(num_cell * num_surf);   // use pointer instead of iterator to access elements 
    AlignedVec fft_in, fft_out;
    fft_in.reserve(fft_offset.size()*fft_size);
    fft_out.reserve(ifft_offset.size()*fft_size);
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < num_cell; i++) 
    {
        for (int j = 0; j < num_surf; j++) 
        {
            all_up_equiv[i*num_surf+j] = cs[i].q_equiv[j];
            all_dn_equiv[i*num_surf+j] = cs[i].p_check[j];
        }
    }

    // up
    tbegin(up);
    AlignedVec fftw_in(N3 * NCHILD);
    AlignedVec fftw_out(fft_size);
    int dim[3] = {N, N, N};
    fftw_plan plan_up = fftw_plan_many_dft_r2c(3, dim, NCHILD,
                                          (real*)&fftw_in[0], nullptr, 1, N3,
                                          (fftw_complex*)(&fftw_out[0]), nullptr, 1, N_freq,
                                          FFTW_ESTIMATE);

    printf("fft_offset.size() = %ld\n", fft_offset.size());
    #pragma omp parallel for
    for (size_t node_idx=0; node_idx<fft_offset.size(); node_idx++) {
      std::vector<real> buffer(fft_size, 0);
      real* up_equiv = &all_up_equiv[fft_offset[node_idx]];
      real* up_equiv_f = &fft_in[fft_size*node_idx];
      std::memset(up_equiv_f, 0, fft_size*sizeof(real));
      for (int k=0; k<num_surf; k++) 
      {
        rtfmm::vec3i idx3 = surf_conv_map[k];
        int idx = idx3[2] * N * N + idx3[1] * N + idx3[0];
        int cpidx = src_cell_idxs[node_idx];
        Cell3& cp = cs[cpidx];
        for (int i=0; i<cp.crange.number; i++)
        {
            int octant = cs[cp.crange.offset + i].octant;
            up_equiv_f[idx + octant * N3] = up_equiv[i*num_surf+k];
        }
      }
      fftw_execute_dft_r2c(plan_up, up_equiv_f, (fftw_complex*)&buffer[0]);
      for (int k=0; k<N_freq; k++) {
        for (int j=0; j<NCHILD; j++) {
          up_equiv_f[2*(NCHILD*k+j)+0] = buffer[2*(N_freq*j+k)+0];
          up_equiv_f[2*(NCHILD*k+j)+1] = buffer[2*(N_freq*j+k)+1];
        }
      }
    }
    fftw_destroy_plan(plan_up);
    tend(up);

    // hadamard
    tbegin(hadamard);
    AlignedVec zero_vec0(fft_size, 0.);
    AlignedVec zero_vec1(fft_size, 0.);
    int npos = ccGks.size();
    int nblk_inter = interaction_count_offset.size();
    int BLOCK_SIZE = CACHE_SIZE * 2 / sizeof(real);
    printf("npos = %d, BLOCK_SIZE = %d\n", npos, BLOCK_SIZE);
    std::vector<real*> IN_(nblk_inter * BLOCK_SIZE);
    std::vector<real*> OUT_(nblk_inter * BLOCK_SIZE);

    #pragma omp parallel for
    for(size_t i = 0; i < fft_out.capacity() / fft_size; i++)
    {
        std::memset(fft_out.data() + i * fft_size, 0, fft_size * sizeof(real));
    }
    #pragma omp parallel for
    for(int i = 0; i < nblk_inter; i++)
    {
        size_t offset0 = (i == 0 ? 0 : interaction_count_offset[i - 1]);
        size_t offset1 = interaction_count_offset[i];
        size_t interaction_count = offset1 - offset0;
        for(size_t j = 0; j < interaction_count; j++)
        {
            IN_[i * BLOCK_SIZE + j]  = &fft_in [interaction_offset_f[2*(offset0 + j) + 0]];
            OUT_[i * BLOCK_SIZE + j] = &fft_out[interaction_offset_f[2*(offset0 + j) + 1]];
        }
        IN_[i * BLOCK_SIZE + interaction_count]  = &zero_vec0[0];
        OUT_[i * BLOCK_SIZE + interaction_count] = &zero_vec1[0];
    }
    #pragma omp parallel for
    for(int f = 0; f < N_freq; f++)
    {
        for(int ipos = 0; ipos < npos; ipos++)
        {
            size_t offset0 = (ipos == 0 ? 0 : interaction_count_offset[ipos - 1]);
            size_t offset1 = interaction_count_offset[ipos]; 
            size_t count = offset1 - offset0;
            real** IN  = &IN_ [ipos * BLOCK_SIZE];
            real** OUT = &OUT_[ipos * BLOCK_SIZE];
            real* M = &ccGks[ipos][f * 2 * NCHILD * NCHILD];
            /*for(size_t j = 0; j < count; j+=2)
            {
                real* M_ = M;
                real* IN0  = IN [j+0] + f * NCHILD * 2;
                real* IN1  = IN [j+1] + f * NCHILD * 2;
                real* OUT0 = OUT[j+0] + f * NCHILD * 2;
                real* OUT1 = OUT[j+1] + f * NCHILD * 2;
                matmult_8x8x2(M_, IN0, IN1, OUT0, OUT1);
                //matmult_8x8x2_avx(M_, IN0, IN1, OUT0, OUT1);
            }*/

            for(size_t j = 0; j < count; j++)
            {
                real* M_ = M;
                real* IN0  = IN [j+0] + f * NCHILD * 2; //src
                real* OUT0 = OUT[j+0] + f * NCHILD * 2; //tar
                //matmult_8x8x1(M_, IN0, OUT0);
                matmult_8x8x1_naive(M_, IN0, OUT0);
            }
        }
    }
    tend(hadamard);

    tbegin(down);
    AlignedVec fftw_in2(fft_size);
    AlignedVec fftw_out2(N3 * NCHILD);
    fftw_plan plan_down = fftw_plan_many_dft_c2r(
        3, dim, NCHILD,
        (fftw_complex*)&fftw_in[0], nullptr, 1, N_freq,
        (real*)&fftw_out[0], nullptr, 1, N3,
        FFTW_ESTIMATE
    );
    #pragma omp parallel for
    for(int i = 0; i < ifft_offset.size(); i++)
    {
        std::vector<real> buffer0(fft_size, 0);
        std::vector<real> buffer1(fft_size, 0);
        real* dn_check_f = &fft_out[i * fft_size];
        real* dn_equiv = &all_dn_equiv[ifft_offset[i]];
        for(int f = 0; f < N_freq; f++)
        {
            for(int j = 0; j < NCHILD; j++)
            {
                buffer0[2 * j * N_freq + 2 * f + 0] = dn_check_f[2 * f * NCHILD + 2 * j + 0];
                buffer0[2 * j * N_freq + 2 * f + 1] = dn_check_f[2 * f * NCHILD + 2 * j + 1];
            }
        }
        fftw_execute_dft_c2r(plan_down, (fftw_complex*)&buffer0[0], (real*)&buffer1[0]);
        int tar_idx = tar_cell_idxs[i];
        //std::cout<<"tar_idx = "<<tar_idx<<std::endl;
        Cell3& ctar = cs[tar_idx];
        real scale = std::pow(2, ctar.depth + 1);
        for(int k = 0; k < num_surf; k++)
        {
            rtfmm::vec3i idx3 = surf_conv_map[k] + rtfmm::vec3i(P-1,P-1,P-1);
            int idx = idx3[2] * N * N + idx3[1] * N + idx3[0];
            for(int j = 0; j < ctar.crange.number; j++)
            {
                int octant = cs[ctar.crange.offset + j].octant;
                dn_equiv[j * num_surf + k] += buffer1[N3 * octant + idx] / N3 * scale;
            }
        }
    }
    fftw_destroy_plan(plan_down);
    tend(down);

    #pragma omp parallel for
    for(int i = 0; i < num_cell; i++)
    {
        for(int j = 0; j < num_surf; j++)
        {
            cs[i].p_check[j] = all_dn_equiv[i * num_surf + j];
        }
    }
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
    std::map<int, int> m2l_gk_relx_idx_map; // hash(rel_x) -> idx
    int m2l_gk_relx_idx_map_cnt = 0;
    for(int k = -range; k <= range; k++)
    {
        for(int j = -range; j <= range; j++)
        {
            for(int i = -range; i <= range; i++)
            {
                if(std::abs(i) > 1 || std::abs(j) > 1 || std::abs(k) > 1)
                {
                    vec3i relx(i,j,k);
                    m2l_gk_relx_idx_map[hash(relx)] = m2l_gk_relx_idx_map_cnt;
                    //std::cout<<hash(relx)<<"->"<<m2l_gk_relx_idx_map_cnt<<std::endl;
                    m2l_gk_relx_idx_map_cnt++;
                }
            }
        }
    }
    std::cout<<"m2l_gk_relx_idx_map_cnt = "<<m2l_gk_relx_idx_map_cnt<<std::endl;
    std::vector<std::vector<rtfmm::complexr>> m2l_Gks_ordered(m2l_gk_relx_idx_map_cnt);
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
                    vec3r relative_pos = vec3r(i,j,k) * r0 * 2; // relative_pos mean src's position
                    std::vector<rtfmm::vec3r> grid = get_conv_grid(N, gsize, delta, relative_pos);
                    std::vector<real> G = get_G_matrix(grid, N); 
                    std::vector<complexr> Gk(N_freq);                   
                    fftw_execute_dft_r2c(plan_G, G.data(), reinterpret_cast<fftw_complex*>(Gk.data()));
                    #pragma omp critical // std::map is not thread-safe
                    {
                        int hv = hash(relative_idx);
                        m2l_Gks[hv] = Gk;
                        m2l_Gk_idx[hv].first = Gk;
                        m2l_Gk_idx[hv].second = m2l_gk_relx_idx_map[hv];
                        m2l_Gks_ordered[m2l_gk_relx_idx_map[hv]] = Gk;
                    }
                }
            }
        }
    }
    std::cout<<"m2l_Gk_idx.size() = "<<m2l_Gk_idx.size()<<std::endl;
    fftw_destroy_plan(plan_G);


    // precompute for t
    std::vector<vec3r> src_rels;
    for(int z = -1; z <= 1; z++)
    {
        for(int y = -1; y <= 1; y++)
        {
            for(int x = -1; x <= 1; x++)
            {
                if(std::abs(x) >= 1 || std::abs(y) >= 1 || std::abs(z) >= 1)
                {
                    src_rels.push_back(vec3r(x,y,z));
                }
            }
        }
    }
    int m2l_par_rel_num = src_rels.size(); //26

    /**
     * @brief gk8x8_map[5][17] means Gk matrix when Ci's 2nd child interacts to Cj's 1st child(namely, 17=2*8+1), 
     * where Ci and Cj are two neighbour cells in pattern-5 relative position(of 26 patterns for nonperiodic m2l).
    */
    std::vector<std::vector<int>> gk8x8_map;
    gk8x8_map = std::vector<std::vector<int>>(m2l_par_rel_num, std::vector<int>(8*8));
    for(int i = 0; i < m2l_par_rel_num; i++)
    {
        vec3r relative_x = src_rels[i] * 2;
        for(int isrc = 0; isrc < 8; isrc++)
        {
            vec3r xsrc = Tree::get_child_cell_x(relative_x,1,isrc,0);
            for(int itar = 0; itar < 8; itar++)
            {
                int ccidx = isrc * 8 + itar;
                vec3r xtar = Tree::get_child_cell_x(vec3r(0,0,0),1,itar,0);
                vec3i dx = (xsrc - xtar).round();
                int hv = hash(dx);
                if(m2l_gk_relx_idx_map.find(hv) != m2l_gk_relx_idx_map.end())
                {
                    int gk_idx = m2l_gk_relx_idx_map[hash(dx)];
                    gk8x8_map[i][ccidx] = gk_idx;
                }
                else
                {
                    gk8x8_map[i][ccidx] = -1;
                }
            }
        }
    }

    ccGks.resize(m2l_par_rel_num, AlignedVec(N_freq * 2 * 8 * 8)); 
    #pragma omp parallel for
    for(int i = 0; i < m2l_par_rel_num; i++)
    {
        for(int j = 0; j < 8 * 8; j++)
        {
            int gk_idx = gk8x8_map[i][j];
            if(gk_idx != -1) // keep GK of neighbour cells to all-zero
            {
                for(int f = 0; f < N_freq; f++)
                {
                    int idx = f * 2 * 8 * 8 + j * 2;
                    ccGks[i][idx + 0] = m2l_Gks_ordered[gk_idx][f].r;
                    ccGks[i][idx + 1] = m2l_Gks_ordered[gk_idx][f].i;
                }
            }
        }
    }
}