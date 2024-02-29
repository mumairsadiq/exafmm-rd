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
    tree.build(bs, args.x, args.r, args.ncrit, args.rega, Tree::TreeType::reg_nonuniform);
    //tree.build(bs, args.x, args.r, args.ncrit, args.rega, Tree::TreeType::nonuniform);
    //tree.build(bs, args.x, args.r, 3, Tree::TreeType::uniform);
    if(args.check_tree)
    {
        check_tree(tree.get_cells());
    }
    tend(build_tree);

    /* traverse to get interaction list */
    tbegin(traverse);
    traverser.traverse(tree, args.cycle, args.images);
    cs = traverser.get_cells();
    tend(traverse);
    tend(build_and_traverse);

    std::cout<<"cs.size() = "<<cs.size()<<std::endl;

    if(verbose) std::cout<<"tree_depth_range = "<<traverser.tree_depth_range<<std::endl;
    if(args.check_tree)
    {
        check_traverser(traverser);
        check_cells(cs);
    }

    tbegin(init_cell_matrix);
    init_cell_bodies(bs, cs, traverser.rw_body_cell_idxs, args.rega);
    init_cell_matrix(cs);
    tend(init_cell_matrix);

    TIME_BEGIN(precompute);
    kernel.precompute(args.P, args.r, args.images);
    kernel.precompute_m2l(args.P, args.r, cs, traverser.M2L_map, args.images);
    if(args.timing) {TIME_END(precompute);}

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

    gather_cell_bodies(bs, cs, traverser.rw_body_cell_idxs);

    if(args.dipole_correction) dipole_correction(bs, args.cycle);
    if(args.divide_4pi) scale_bodies(bs);

    return sort_bodies_by_idx(bs);
}

void rtfmm::LaplaceFMM::P2M()
{
    Indices lcidx = traverser.leaf_cell_idxs;
    int leaf_number = lcidx.size();
    if(verbose) std::cout<<"P2M leaf_number = "<<leaf_number<<std::endl;
    TIME_BEGIN(P2M);
    #pragma omp parallel for
    for(int i = 0; i < leaf_number; i++)
    {
        Cell3& c = cs[lcidx[i]];
        kernel.p2m_precompute(args.P, c);
    }
    if(args.timing){TIME_END(P2M);}
}

void rtfmm::LaplaceFMM::M2M()
{
    TIME_BEGIN(M2M);
    for(int depth = traverser.tree_depth_range[1]; depth > traverser.tree_depth_range[0]; depth--)
    {
        Indices& nlcidx = traverser.nonleaf_cell_idxs[depth];
        int num = nlcidx.size();
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            Cell3& ci = cs[nlcidx[i]];
            if(depth >= 0)
                kernel.m2m_precompute(args.P, ci, cs);
            else
                kernel.m2m_img_precompute(args.P, ci, cs, args.cycle * std::pow(3, -depth - 1));
        }
    }
    if(args.timing){TIME_END(M2M);}
}

void rtfmm::LaplaceFMM::M2L()
{
    TIME_BEGIN(M2L);
    //PeriodicInteractionPairs m2l_pairs = traverser.M2L_pairs;
    PeriodicInteractionMap m2l_map = traverser.M2L_map;
    PeriodicM2LMap& m2l_parent_map = traverser.M2L_parent_map;
    //PeriodicInteractionMap m2l_map2 = traverser.get_m2l_map_from_m2l_parent_map();

    //kernel.m2l_fft_precompute_advanced2(args.P, cs, m2l_map2);
    kernel.m2l_fft_precompute_advanced2(args.P, cs, m2l_map);
    //kernel.m2l_fft_precompute_advanced3(args.P, cs, m2l_map, m2l_pairs);
    kernel.m2l_fft_precompute_t(args.P, cs, m2l_parent_map);

    if(args.timing){TIME_END(M2L);}
}

void rtfmm::LaplaceFMM::L2L()
{
    TIME_BEGIN(L2L);
    for(int depth = traverser.tree_depth_range[0] + 1; depth < traverser.tree_depth_range[1]; depth++)
    {
        printf("depth = %d\n", depth);
        Indices& nlcidx = traverser.nonleaf_cell_idxs[depth];
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
    if(args.timing){TIME_END(L2L);}
}

void rtfmm::LaplaceFMM::L2P()
{
    TIME_BEGIN(L2P);
    Indices lcidx = traverser.leaf_cell_idxs;
    int leaf_number = lcidx.size();
    if(verbose) std::cout<<"L2P leaf_number = "<<leaf_number<<std::endl;
    #pragma omp parallel for
    for(int i = 0; i < leaf_number; i++)
    {
        Cell3& c = cs[lcidx[i]];
        kernel.l2p_precompute(args.P, c);
    }
    if(args.timing){TIME_END(L2P);}
}

void rtfmm::LaplaceFMM::P2P()
{
    TIME_BEGIN(P2P);
    PeriodicInteractionMap& p2p_map = traverser.P2P_map;
    std::cout<<"p2p_map.size = "<<p2p_map.size()<<std::endl;
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
    PeriodicInteractionPairs& m2p_pairs = traverser.M2P_pairs;
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
    PeriodicInteractionPairs& p2l_pairs = traverser.P2L_pairs;
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
    #pragma omp parallel for
    for(int i = 0; i < cells.size(); i++)
    {
        Cell3& cell = cells[i];
        cell.q_equiv = Matrix(get_surface_point_num(args.P), 1);
        cell.p_check = Matrix(get_surface_point_num(args.P), 1);
    }
}

void rtfmm::LaplaceFMM::init_cell_bodies(Bodies3& bs, Cells3& cells, Indices& idxs, real rega)
{
    //#pragma omp parallel for
    for(int i = 0; i < idxs.size(); i++)
    {
        Cell3& c = cells[idxs[i]];
        c.bs.clear();
        c.ws.clear();
        if(rega == 0)
        {
            for(int i = 0; i < c.brange.number; i++)
            {
                c.bs.push_back(bs[c.brange.offset + i]);
                c.ws.push_back(1.0);
            }
        }
        else
        {
            /*for(int i = 0; i < bs.size(); i++)
            {
                Body3& b = bs[i];
                real w = LaplaceKernel::get_rega_w(b.x - c.x, c.r, rega);
                if(w > 0)
                {
                    c.bs.push_back(b);
                    c.ws.push_back(w);
                }
                if(w > 0 && w != 1)
                {
                    print_body(b);
                    printf("w = %.4f\n", w);
                }
            }*/
            for(int i = 0; i < c.brange.number; i++)
            {
                Body3& b = bs[c.brange.offset + i];
                c.bs.push_back(b);
                real w = LaplaceKernel::get_rega_w(b.x - c.x, c.r, rega);
                c.ws.push_back(w);
                /*if(w != 1)
                {
                    print_body(b);
                    std::cout<<c<<std::endl;
                    printf("w = %.4f\n", w);
                }*/
            }
            for(int i = 0; i < c.reg_bs.size(); i++)
            {
                Body3& b = c.reg_bs[i];
                c.bs.push_back(b);
                real w = LaplaceKernel::get_rega_w(b.x - c.x, c.r, rega);
                c.ws.push_back(w);
                /*if(w != 1)
                {
                    print_body(b);
                    std::cout<<c<<std::endl;
                    printf("w = %.4f\n", w);
                }*/
            }
        }
    }
}

void rtfmm::LaplaceFMM::gather_cell_bodies(Bodies3& bs, Cells3& cells, Indices& idxs)
{
    /*for(auto idx : idxs)
    {
        Cell3& c = cells[idx];
        for(int i = 0; i < c.brange.number; i++)
        {
            bs[c.brange.offset + i].p += c.bs[i].p;
            bs[c.brange.offset + i].f += c.bs[i].f;
        }
    }*/

    bs = sort_bodies_by_idx(bs);
    for(auto idx : idxs)
    {
        Cell3& c = cells[idx];
        for(int i = 0; i < c.bs.size(); i++)
        {
            Body3& b = c.bs[i];
            bs[b.idx].p += b.p;
            bs[b.idx].f += b.f;
        }
    }
}

void rtfmm::LaplaceFMM::check_tree(const Cells3& cells)
{
    printf("check tree\n");
    printf("cells.size() = %ld\n", cells.size());
    int bcnt = 0;
    for(int i = 0; i < cells.size(); i++)
    {
        Cell3 c = cells[i];
        for(int j = 0; j < c.brange.number; j++)
        {
            Body3 b = bs[c.brange.offset + j];
            vec3r dx = (b.x - c.x).abs();
            assert_exit(dx[0] <= c.r && dx[1] <= c.r && dx[2] <= c.r, "tree body range error");
        }
        if(c.crange.number == 0)
        {
            bcnt += c.brange.number;
        }
    }
    assert_exit(bcnt == bs.size(), "tree body sum error");
}

void rtfmm::LaplaceFMM::check_traverser(Traverser& traverser)
{
    PeriodicInteractionPairs& P2P_pairs = traverser.P2P_pairs;
    PeriodicInteractionPairs& M2L_pairs = traverser.M2L_pairs;
    PeriodicInteractionPairs& M2P_pairs = traverser.M2P_pairs;
    PeriodicInteractionPairs& P2L_pairs = traverser.P2L_pairs;

    int p2p_num = P2P_pairs.size();
    int m2l_num = M2L_pairs.size();
    int m2p_num = M2P_pairs.size();
    int p2l_num = P2L_pairs.size();

    std::cout<<"P2P_pair# : "<<p2p_num<<std::endl;
    std::cout<<"M2L_pair# : "<<m2l_num<<std::endl;
    std::cout<<"M2P_pair# : "<<m2p_num<<std::endl;
    std::cout<<"P2L_pair# : "<<p2l_num<<std::endl;
    std::cout<<"sum# : "<<p2p_num+m2l_num+m2p_num+p2l_num<<std::endl;

    // check pairs
    if(args.images == 0)
    {
        for(auto p2p : P2P_pairs)
        {
            int found = 0;
            for(auto p2p2 : P2P_pairs)
            {
                if(p2p.first == p2p2.second.first && p2p.second.first == p2p2.first)
                {
                    found++;
                }
            }
            assert_exit(found == 1, "treep p2p error");
        }
        std::cout<<"p2p check ok"<<std::endl;
        for(auto m2l : M2L_pairs)
        {
            int found = 0;
            for(auto m2l2 : M2L_pairs)
            {
                if(m2l.first == m2l2.second.first && m2l.second.first == m2l2.first)
                {
                    found++;
                }
            }
            assert_exit(found == 1, "treep m2l error");
        }
        std::cout<<"m2l check ok"<<std::endl;
    }

    for(auto m2p : M2P_pairs)
    {
        int found = 0;
        for(auto p2l : P2L_pairs)
        {
            if(m2p.first == p2l.second.first && m2p.second.first == p2l.first)
            {
                found++;
            }
        }
        assert_exit(found == 1, "treep p2l error");
    }
    std::cout<<"m2p check ok"<<std::endl;

    for(auto p2l : P2L_pairs)
    {
        int found = 0;
        for(auto m2p : M2P_pairs)
        {
            if(m2p.first == p2l.second.first && m2p.second.first == p2l.first)
            {
                found++;
            }
        }
        assert_exit(found == 1, "treep m2p error");
    }
    std::cout<<"p2l check ok"<<std::endl;
}

void rtfmm::LaplaceFMM::check_cells(const Cells3& cells)
{
    std::cout<<"check cells"<<std::endl;
    for(int i = 0; i < cells.size(); i++)
    {
        Cell3 c = cells[i];
        assert_exit(c.idx == i, "cell idx error");
    }
}