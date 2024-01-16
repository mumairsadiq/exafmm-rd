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
    Tree tree;
    tree.build(bs, args.x, args.r, args.ncrit, Tree::TreeType::nonuniform);
    cs = tree.get_cells();
    if(verbose && args.check_tree) check_tree(cs);

    /* traverse to get interaction list */
    traverser.traverse(tree, args.cycle, args.images);
    cs = traverser.get_cells();
    if(verbose && args.check_tree) check_traverser(traverser);
    if(verbose && args.check_tree) check_cells(cs);

    init_cell_matrix(cs);

    if(args.use_precompute)
    {
        TIME_BEGIN(precompute);
        kernel.precompute(args.P, args.r, args.images);
        kernel.precompute_m2l(args.P, args.r, cs, traverser.get_map(OperatorType::M2L), args.images);
        if(args.timing) {TIME_END(precompute);}
    }

    P2M();
    M2M();
    M2L();
    P2L();
    L2L();
    L2P();
    M2P();
    P2P();


    #ifdef TEST_TREE
    for(auto b : bs)
    {
        if(b.p != bs.size())
        {
            printf("error! b.p = %.4f\n", b.p);
            break;
        }
    }
    printf("test passed!\n");
    exit(0);
    #endif

    dipole_correction(bs, args.cycle);
        
    return sort_bodies_by_idx(bs);
}

void rtfmm::LaplaceFMM::P2M()
{
    TIME_BEGIN(P2M);
    Indices lcidx = get_leaf_cell_indices(cs);
    int leaf_number = lcidx.size();
    for(int i = 0; i < leaf_number; i++)
    {
        Cell3& c = cs[lcidx[i]];
        #ifndef TEST_TREE
        if(args.use_precompute)
            kernel.p2m_precompute(args.P, bs, c);
        else
            kernel.p2m(args.P, bs, c);
        #else
        for(int i = 0; i < c.brange.number; i++)
        {
            c.M += 1;
        }
        #endif
    }
    if(args.timing)
    {
        TIME_END(P2M);
    }
}

void rtfmm::LaplaceFMM::M2M()
{
    TIME_BEGIN(M2M);
    int min_depth = get_min_depth(cs);
    int max_depth = get_max_depth(cs);
    if(verbose) printf("min_depth = %d, max_depth = %d\n", min_depth, max_depth);
    for(int depth = max_depth; depth >= min_depth; depth--)
    {
        Indices nlcidx = get_nonleaf_cell_indices(cs, depth);
        int num = nlcidx.size();
        for(int i = 0; i < num; i++)
        {
            Cell3& ci = cs[nlcidx[i]];
            #ifndef TEST_TREE
            if(depth >= 0)
            {
                if(args.use_precompute)
                    kernel.m2m_precompute(args.P, ci, cs);
                else
                    kernel.m2m(args.P, ci, cs);
            }
            else
            {
                if(args.use_precompute)
                    kernel.m2m_img_precompute(args.P, ci, cs, args.cycle * std::pow(3, -depth - 1));
                else
                    kernel.m2m_img(args.P, ci, cs, args.cycle * std::pow(3, -depth - 1));
            }
            #else
            for(int j = 0; j < ci.crange.number; j++)
            {
                Cell3& cj = cs[ci.crange.offset + j];
                ci.M += cj.M;
            }
            #endif
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
    PeriodicInteractionMap m2l_map = traverser.get_map(OperatorType::M2L);
    PeriodicInteractionPairs m2l_pairs = traverser.get_pairs(OperatorType::M2L);
    if(args.use_precompute)
    {
        TIME_BEGIN(M2L_kernel);
        kernel.m2l_fft_precompute_advanced2(args.P, cs, m2l_map, m2l_pairs);
        TIME_END(M2L_kernel);
    }
    else
    {
        for(int i = 0; i < m2l_map.size(); i++)
        {
            auto m2l_list = m2l_map.begin();
            std::advance(m2l_list, i);
            Cell3& ctar = cs[m2l_list->first];
            for(int j = 0; j < m2l_list->second.size(); j++)
            {
                Cell3& csrc = cs[m2l_list->second[j].first];
                vec3r offset_src = m2l_list->second[j].second;
                if(args.use_fft) 
                {
                    kernel.m2l_fft(args.P, csrc, ctar, offset_src);
                }
                else
                {
                    kernel.m2l(args.P, csrc, ctar, offset_src);
                }
            }
        }
    }
    if(args.timing)
    {
        TIME_END(M2L);
    }
}

void rtfmm::LaplaceFMM::L2L()
{
    TIME_BEGIN(L2L);
    int min_depth = get_min_depth(cs);
    int max_depth = get_max_depth(cs);
    if(verbose) printf("min_depth = %d, max_depth = %d\n", min_depth, max_depth);
    for(int depth = min_depth; depth <= max_depth; depth++)
    {
        Indices nlcidx = get_nonleaf_cell_indices(cs, depth);
        int num = nlcidx.size();
        for(int i = 0; i < num; i++)
        {
            Cell3& ci = cs[nlcidx[i]];
            #ifndef TEST_TREE
            if(args.use_precompute)
            {
                if(depth >= 0)
                    kernel.l2l_precompute(args.P, ci, cs);
                else
                    kernel.l2l_img_precompute(args.P, ci, cs);
            }
            else
                kernel.l2l(args.P, ci, cs);
            #else
            for(int j = 0; j < ci.crange.number; j++)
            {
                Cell3& cj = cs[ci.crange.offset + j];
                cj.L += ci.L;
            }
            #endif
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
    for(int i = 0; i < leaf_number; i++)
    {
        Cell3& c = cs[lcidx[i]];
        #ifndef TEST_TREE
        if(args.use_precompute)
            kernel.l2p_precompute(args.P, bs, c);
        else
            kernel.l2p(args.P, bs, c);
        #else
        for(int j = 0; j < c.brange.number; j++)
        {
            bs[c.brange.offset + j].p += c.L;
        }
        #endif
    }
    if(args.timing)
    {
        TIME_END(L2P);
    }
}

void rtfmm::LaplaceFMM::P2P()
{
    TIME_BEGIN(P2P);
    PeriodicInteractionPairs p2p_pairs = traverser.get_pairs(OperatorType::P2P);
    for(auto p2p : p2p_pairs)
    {
        Cell3& ctar = cs[p2p.first];
        Cell3& csrc = cs[p2p.second.first];
        #ifndef TEST_TREE
        kernel.p2p(bs,bs,csrc,ctar,p2p.second.second);
        #else
        for(int j = 0; j < ctar.brange.number; j++)
        {
            bs[ctar.brange.offset + j].p += csrc.brange.number;
        }
        #endif
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
        #ifndef TEST_TREE
        kernel.m2p(args.P, bs, csrc, ctar, m2p.second.second);
        #else
        for(int j = 0; j < ctar.brange.number; j++)
        {
            bs[ctar.brange.offset + j].p += csrc.M;
        }
        #endif
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
        #ifndef TEST_TREE
        kernel.p2l(args.P, bs, csrc, ctar, p2l.second.second);
        #else
        ctar.L += csrc.brange.number;
        #endif
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

int rtfmm::LaplaceFMM::get_min_depth(const Cells3& cells)
{
    int res = cells[0].depth;
    for(int i = 0; i < cells.size(); i++)
    {
        res = std::min(res, cells[i].depth);
    }
    return res;
}

int rtfmm::LaplaceFMM::get_max_depth(const Cells3& cells)
{
    int res = 0;
    for(int i = 0; i < cells.size(); i++)
    {
        res = std::max(res, cells[i].depth);
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
    PeriodicInteractionPairs P2P_pairs = traverser.get_pairs(OperatorType::P2P);
    PeriodicInteractionPairs M2L_pairs = traverser.get_pairs(OperatorType::M2L);
    PeriodicInteractionPairs M2P_pairs = traverser.get_pairs(OperatorType::M2P);
    PeriodicInteractionPairs P2L_pairs = traverser.get_pairs(OperatorType::P2L);

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