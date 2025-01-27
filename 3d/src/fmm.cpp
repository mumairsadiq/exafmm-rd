#include "fmm.h"
#include "tree.h"
#include "mathfunc.h"
#include <omp.h>
#include "surface.h"

rtfmm::LaplaceFMM::LaplaceFMM(const Bodies3& bs_, const Argument& args_) 
    : bs(bs_), args(args_)
{
    assert(bs.size() == args.n);
}

rtfmm::Bodies3 rtfmm::LaplaceFMM::solve()
{
    /* build tree */
    Tree tree;
    tree.build(bs, args.x, args.r, args.ncrit, Tree::TreeType::nonuniform);
    cs = tree.get_cells();
    check_tree(cs);
    init_cell_matrix(cs);

    /* traverse to get interaction list */
    traverser.traverse(tree);
    check_traverser(traverser);

    /* upward */
    P2M();
    M2M();
    /* downward */
    M2L();
    P2L();
    L2L();
    L2P();
    M2P();
    P2P();

    // check
    /*for(auto b : bs)
    {
        assert(b.p == bs.size());
    }*/


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
        kernel.p2m(args.P, bs, c);
        /*for(int i = 0; i < c.brange.number; i++)
        {
            c.M += 1;
        }*/
    }
    if(args.timing)
    {
        TIME_END(P2M);
    }
}

void rtfmm::LaplaceFMM::M2M()
{
    TIME_BEGIN(M2M);
    int max_depth = get_max_depth(cs);
    for(int depth = max_depth; depth >= 0; depth--)
    {
        Indices nlcidx = get_nonleaf_cell_indices(cs, depth);
        int num = nlcidx.size();
        for(int i = 0; i < num; i++)
        {
            Cell3& ci = cs[nlcidx[i]];
            kernel.m2m(args.P, ci, cs);
            /*for(int j = 0; j < ci.crange.number; j++)
            {
                Cell3& cj = cs[ci.crange.offset + j];
                ci.M += cj.M;
            }*/
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
    InteractionPairs m2l_pairs = traverser.get_pairs(OperatorType::M2L);
    for(auto m2l : m2l_pairs)
    {
        kernel.m2l(args.P, cs[m2l.first], cs[m2l.second]);
        //cs[m2l.first].L += cs[m2l.second].M;
    }
    if(args.timing)
    {
        TIME_END(M2L);
    }
}

void rtfmm::LaplaceFMM::L2L()
{
    TIME_BEGIN(L2L);
    int max_depth = get_max_depth(cs);
    for(int depth = 0; depth <= max_depth; depth++)
    {
        Indices nlcidx = get_nonleaf_cell_indices(cs, depth);
        int num = nlcidx.size();
        for(int i = 0; i < num; i++)
        {
            Cell3& ci = cs[nlcidx[i]];
            kernel.l2l(args.P, ci, cs);
            /*for(int j = 0; j < ci.crange.number; j++)
            {
                Cell3& cj = cs[ci.crange.offset + j];
                cj.L += ci.L;
            }*/
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
        kernel.l2p(args.P, bs, c);
        /*for(int j = 0; j < c.brange.number; j++)
        {
            bs[c.brange.offset + j].p += c.L;
        }*/
    }
    if(args.timing)
    {
        TIME_END(L2P);
    }
}

void rtfmm::LaplaceFMM::P2P()
{
    TIME_BEGIN(P2P);
    InteractionPairs p2p_pairs = traverser.get_pairs(OperatorType::P2P);
    for(auto p2p : p2p_pairs)
    {
        Cell3& ctar = cs[p2p.first];
        Cell3& csrc = cs[p2p.second];
        kernel.p2p_matrix(bs,bs,csrc,ctar);
        
        /*for(int j = 0; j < ctar.brange.number; j++)
        {
            bs[ctar.brange.offset + j].p += csrc.brange.number;
        }*/
    }
    if(args.timing)
    {
        TIME_END(P2P);
    }
}

void rtfmm::LaplaceFMM::M2P()
{
    TIME_BEGIN(M2P);
    InteractionPairs m2p_pairs = traverser.get_pairs(OperatorType::M2P);
    for(auto m2p : m2p_pairs)
    {
        Cell3& ctar = cs[m2p.first];
        Cell3& csrc = cs[m2p.second];
        kernel.m2p(args.P, bs, csrc, ctar);
        /*for(int j = 0; j < ctar.brange.number; j++)
        {
            bs[ctar.brange.offset + j].p += csrc.M;
        }*/
    }
    if(args.timing)
    {
        TIME_END(M2P);
    }
}

void rtfmm::LaplaceFMM::P2L()
{
    TIME_BEGIN(P2L);
    InteractionPairs p2l_pairs = traverser.get_pairs(OperatorType::P2L);
    for(auto p2l : p2l_pairs)
    {
        Cell3& ctar = cs[p2l.first];
        Cell3& csrc = cs[p2l.second];
        kernel.p2l(args.P, bs, csrc, ctar);
        //ctar.L += csrc.brange.number;
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
    printf("cells.size() = %ld\n", cells.size());
    int bcnt = 0;
    for(int i = 0; i < cells.size(); i++)
    {
        Cell3 c = cells[i];
        for(int j = 0; j < c.brange.number; j++)
        {
            Body3 b = bs[c.brange.offset + j];
            vec3r dx = (b.x - c.x).abs();
            assert(dx[0] <= c.r && dx[1] <= c.r && dx[2] <= c.r);
        }
        if(c.crange.number == 0)
        {
            bcnt += c.brange.number;
        }
    }
    assert(bcnt == bs.size());
}

void rtfmm::LaplaceFMM::check_traverser(Traverser& traverser)
{
    InteractionPairs P2P_pairs = traverser.get_pairs(OperatorType::P2P);
    InteractionPairs M2L_pairs = traverser.get_pairs(OperatorType::M2L);
    InteractionPairs M2P_pairs = traverser.get_pairs(OperatorType::M2P);
    InteractionPairs P2L_pairs = traverser.get_pairs(OperatorType::P2L);

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
    for(auto m2p : M2P_pairs)
    {
        int found = 0;
        for(auto p2l : P2L_pairs)
        {
            if(m2p.first == p2l.second && m2p.second == p2l.first)
            {
                found++;
            }
        }
        assert(found == 1);
    }

    for(auto p2l : P2L_pairs)
    {
        int found = 0;
        for(auto m2p : M2P_pairs)
        {
            if(m2p.first == p2l.second && m2p.second == p2l.first)
            {
                found++;
            }
        }
        assert(found == 1);
    }
}
