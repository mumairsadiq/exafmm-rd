#include "fmm.h"
#include "tree.h"
#include "mathfunc.h"
#include <omp.h>
#include "surface.h"

rtfmm::LaplaceFMM::LaplaceFMM(const Bodies3& bs_, const Argument& args_) 
    : bs(bs_), args(args_)
{
    assert_exit(bs.size() == args.n, "LaplaceFMM init body size error");
}

rtfmm::Bodies3 rtfmm::LaplaceFMM::solve()
{
    /* build tree */
    tbegin(build_and_traverse);
    tbegin(build_tree);
    Tree tree;
    //tree.build(bs, args.x, args.r, args.ncrit, args.rega, Tree::TreeType::nonuniform);
    tree.build(bs, args.x, args.r, args.ncrit, args.rega, Tree::TreeType::regnonuniform);
    cs = tree.get_cells();
    tend(build_tree);

    /* traverse to get interaction list */
    tbegin(traverse);
    traverser.traverse(tree, args.cycle, args.images, args.P);
    cs = traverser.get_cells();
    tree_depth_range = get_min_max_depth(cs);
    if(verbose) std::cout<<"tree_depth_range = "<<tree_depth_range<<std::endl;
    tend(traverse);
    tend(build_and_traverse);

    tbegin(init_cell_matrix);
    init_cell_matrix(cs);
    tend(init_cell_matrix);

    std::cout<<"cell number = "<<cs.size()<<std::endl;

    tbegin(load_leaf_body);
    for(int j = 0; j < cs.size(); j++)
    {
        Cell3& c = cs[j];
        if(c.crange.number == 0)
        {
            for(int i = 0; i < c.brange.number; i++)
            {
                Body3& b = bs[c.brange.offset + i];
                real w = LaplaceKernel::get_rega_w(b.x - c.x, c.r, args.rega);
                c.bs.push_back(b);
                c.ws.push_back(w);
            }
            for(int i = 0; i < c.rbs.size(); i++)
            {
                Body3& b = c.rbs[i];
                real w = LaplaceKernel::get_rega_w(b.x - c.x, c.r, args.rega);
                c.bs.push_back(b);
                c.ws.push_back(w);
            }
        }
    }
    tend(load_leaf_body);

    if(args.use_precompute)
    {
        TIME_BEGIN(precompute);
        TIME_BEGIN(precompute_others);
        kernel.precompute(args.P, args.r, args.images);
        TIME_END(precompute_others);
        TIME_BEGIN(precompute_m2l);
        kernel.precompute_m2l(args.P, args.r, cs, traverser.get_map(OperatorType::M2L), args.images);
        TIME_END(precompute_m2l);
        if(args.timing) {TIME_END(precompute);}
    }

    tbegin(FMM_kernels);
    P2M();
    M2M();
    M2L();
    P2L();
    L2L();
    L2P();
    M2P();
    P2P();
    tend(FMM_kernels);

    // merge
    bs = sort_bodies_by_idx(bs);
    for(int j = 0; j < cs.size(); j++)
    {
        Cell3& c = cs[j];
        if(c.crange.number == 0)
        {
            for(int i = 0; i < c.bs.size(); i++)
            {
                Body3& b = c.bs[i];
                bs[b.idx].p += b.p;
                bs[b.idx].f += b.f;
            }
        }
    }

    if(args.dipole_correction)
        dipole_correction(bs, args.cycle);
        
    if(args.divide_4pi)
        scale_bodies(bs);
    //return sort_bodies_by_idx(bs);
    return bs;
}

void rtfmm::LaplaceFMM::P2M()
{
    Indices lcidx = get_leaf_cell_indices(cs);
    int leaf_number = lcidx.size();
    if(verbose) std::cout<<"P2M leaf_number = "<<leaf_number<<std::endl;
    TIME_BEGIN(P2M);
    #pragma omp parallel for
    for(int i = 0; i < leaf_number; i++)
    {
        Cell3& c = cs[lcidx[i]];
        kernel.p2m_precompute(args.P, c);
    }
    if(args.timing)
    {
        TIME_END(P2M);
    }
}

void rtfmm::LaplaceFMM::M2M()
{
    TIME_BEGIN(M2M);
    for(int depth = tree_depth_range[1]; depth > tree_depth_range[0]; depth--)
    {
        Indices nlcidx = get_nonleaf_cell_indices(cs, depth);
        int num = nlcidx.size();
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            Cell3& ci = cs[nlcidx[i]];
            if(depth >= 0)
            {
                kernel.m2m_precompute(args.P, ci, cs);
            }
            else
            {
                kernel.m2m_img_precompute(args.P, ci, cs, args.cycle * std::pow(3, -depth - 1));
            }
        }
    }
    if(args.timing)
    {
        TIME_END(M2M);
    }
}

void rtfmm::LaplaceFMM::M2L()
{
    TIME_BEGIN(M2L);
    //PeriodicInteractionPairs m2l_pairs = traverser.get_pairs(OperatorType::M2L);
    //PeriodicInteractionMap m2l_map = traverser.get_map(OperatorType::M2L);
    PeriodicM2LMap m2l_parent_map = traverser.get_M2L_parent_map();
    //PeriodicInteractionMap m2l_map2 = traverser.get_m2l_map_from_m2l_parent_map();
    //kernel.m2l_fft_precompute_advanced2(args.P, cs, m2l_map2);
    //kernel.m2l_fft_precompute_advanced2(args.P, cs, m2l_map);
    //kernel.m2l_fft_precompute_advanced3(args.P, cs, m2l_map, m2l_pairs);
    kernel.m2l_fft_precompute_t(args.P, cs, m2l_parent_map);
    if(args.timing)
    {
        TIME_END(M2L);
    }
}

void rtfmm::LaplaceFMM::L2L()
{
    TIME_BEGIN(L2L);
    for(int depth = tree_depth_range[0]; depth < tree_depth_range[1]; depth++)
    {
        Indices nlcidx = get_nonleaf_cell_indices(cs, depth);
        int num = nlcidx.size();
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            Cell3& ci = cs[nlcidx[i]];
            if(depth >= 0)
                kernel.l2l_precompute(args.P, ci, cs);
            else
                kernel.l2l_img_precompute(args.P, ci, cs);
        }
    }
    if(args.timing)
    {
        TIME_END(L2L);
    }
}

void rtfmm::LaplaceFMM::L2P()
{
    TIME_BEGIN(L2P);
    Indices lcidx = get_leaf_cell_indices(cs);
    int leaf_number = lcidx.size();
    if(verbose) std::cout<<"L2P leaf_number = "<<leaf_number<<std::endl;
    #pragma omp parallel for
    for(int i = 0; i < leaf_number; i++)
    {
        Cell3& c = cs[lcidx[i]];
        kernel.l2p_precompute(args.P, c);
    }
    if(args.timing)
    {
        TIME_END(L2P);
    }
}

void rtfmm::LaplaceFMM::P2P()
{
    TIME_BEGIN(P2P);
    PeriodicInteractionMap p2p_map = traverser.get_map(OperatorType::P2P);
    if(verbose) std::cout<<"p2p_pair.size() = "<<traverser.get_pairs(OperatorType::P2P).size()<<std::endl;
    #pragma omp parallel for
    for(int i = 0; i < p2p_map.size(); i++)
    {
        auto p2p = p2p_map.begin();
        std::advance(p2p, i);
        Cell3& ctar = cs[p2p->first];
        kernel.p2p_1toN_256(cs,p2p->second,ctar);
    }
    if(args.timing)
    {
        TIME_END(P2P);
    }
}

void rtfmm::LaplaceFMM::M2P()
{
    TIME_BEGIN(M2P);
    PeriodicInteractionPairs m2p_pairs = traverser.get_pairs(OperatorType::M2P);
    for(auto m2p : m2p_pairs)
    {
        Cell3& ctar = cs[m2p.first];
        Cell3& csrc = cs[m2p.second.first];
        kernel.m2p(args.P, bs, csrc, ctar, m2p.second.second);
    }
    if(args.timing)
    {
        TIME_END(M2P);
    }
}

void rtfmm::LaplaceFMM::P2L()
{
    TIME_BEGIN(P2L);
    PeriodicInteractionPairs p2l_pairs = traverser.get_pairs(OperatorType::P2L);
    for(auto p2l : p2l_pairs)
    {
        Cell3& ctar = cs[p2l.first];
        Cell3& csrc = cs[p2l.second.first];
        kernel.p2l(args.P, bs, csrc, ctar, p2l.second.second);
    }
    if(args.timing)
    {
        TIME_END(P2L);
    }
}

void rtfmm::LaplaceFMM::init_cell_matrix(Cells3& cells)
{
    for(int i = 0; i < cells.size(); i++)
    {
        Cell3& cell = cells[i];
        cell.q_equiv = Matrix(get_surface_point_num(args.P), 1);
        cell.p_check = Matrix(get_surface_point_num(args.P), 1);
    }
}

rtfmm::vec2i rtfmm::LaplaceFMM::get_min_max_depth(const Cells3& cells)
{
    vec2i res = {cells[0].depth, cells[0].depth};
    for(int i = 0; i < cells.size(); i++)
    {
        res[0] = std::min(res[0], cells[i].depth);
        res[1] = std::max(res[1], cells[i].depth);
    }
    return res;
}

rtfmm::Indices rtfmm::LaplaceFMM::get_leaf_cell_indices(const Cells3& cells)
{
    Indices res;
    for(int i = 0; i < cells.size(); i++)
    {
        if(cells[i].crange.number == 0)
        {
            res.push_back(i);
        }
    }
    return res;
}

rtfmm::Indices rtfmm::LaplaceFMM::get_nonleaf_cell_indices(const Cells3& cells, int depth)
{
    Indices res;
    for(int i = 0; i < cells.size(); i++)
    {
        if(cells[i].depth == depth && cells[i].crange.number > 0)
        {
            res.push_back(i);
        }
    }
    return res;
}

