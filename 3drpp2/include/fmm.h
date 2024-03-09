#pragma once
#include "type.h"
#include "body.h"
#include "tree.h"
#include "argument.h"
#include "kernel.h"
#include "traverser.h"

namespace rtfmm
{

class LaplaceFMM
{
public:
    LaplaceFMM(const Bodies3& bs_, const Argument& args_);

    Bodies3 solve();

private:
    Argument args;
    Bodies3 bs;
    Cells3 cs;
    Traverser traverser;
    LaplaceKernel kernel;
    vec2i tree_depth_range; //[min_depth,max_depth]

private:
    void P2M();
    void M2M();
    void M2L();
    void L2L();
    void L2P();
    void M2P();
    void P2L();
    void P2P();

private:
    void init_cell_matrix(Cells3& cells);
    vec2i get_min_max_depth(const Cells3& cells);
    Indices get_leaf_cell_indices(const Cells3& cells);
    Indices get_nonleaf_cell_indices(const Cells3& cells, int depth);
};

}