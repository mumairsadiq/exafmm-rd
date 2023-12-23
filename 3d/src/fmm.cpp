#include "fmm.h"
#include "tree.h"

rtfmm::LaplaceFMM::LaplaceFMM(const Bodies3& bs_, const Argument& args_) 
    : bs(bs_), args(args_)
{
    assert(bs.size() == args.n);
}

rtfmm::Bodies3 rtfmm::LaplaceFMM::solve()
{
    Tree tree;
    tree.build(bs, args.x, args.r, 3, Tree::TreeType::uniform);
    cells = tree.get_cells();
    printf("cells.size() = %ld\n", cells.size());

    kernel.p2p(bs, bs, cells[0], cells[0]);
    //kernel.p2m(args.P, bs, cells[0]);
    //kernel.m2l(args.P, cells[0], cells[0]);
    //kernel.l2p(args.P, bs, cells[0]);

    Bodies3 res = sort_bodies_by_idx(bs);

    return res;
}
