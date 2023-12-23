#include "fmm.h"
#include "tree.h"
#include <omp.h>

rtfmm::LaplaceFMM::LaplaceFMM(const Bodies3& bs_, const Argument& args_) 
    : bs(bs_), args(args_)
{
    assert(bs.size() == args.n);
}

rtfmm::Bodies3 rtfmm::LaplaceFMM::solve()
{
    Tree tree;
    tree.build(bs, args.x, args.r, args.ncrit, Tree::TreeType::nonuniform);
    Cells3 cells = tree.get_cells();
    printf("cells.size() = %ld\n", cells.size());
    check_tree(cells);

    Indices leaf_cell_indices = get_leaf_cell_indices(cells);
    int leaf_number = leaf_cell_indices.size();
    printf("leaf number = %d\n", leaf_number);

    kernel.p2p(bs, bs, cells[0], cells[0]);
    //kernel.p2m(args.P, bs, cells[0]);
    //kernel.m2l(args.P, cells[0], cells[0]);
    //kernel.l2p(args.P, bs, cells[0]);
    printf("max threads = %d\n", omp_get_max_threads());
    TIME_BEGIN(P2M);
    #pragma omp parallel for
    for(int i = 0; i < leaf_number; i++)
    {
        kernel.p2m(args.P, bs, cells[leaf_cell_indices[i]]);
    }
    TIME_END(P2M);

    Bodies3 res = sort_bodies_by_idx(bs);

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

void rtfmm::LaplaceFMM::check_tree(const Cells3& cells)
{
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
    printf("bcnt = %d\n", bcnt);
    assert(bcnt == bs.size());
}
