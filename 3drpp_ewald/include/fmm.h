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

protected:
    Argument args;
    Bodies3 bs;
    Cells3 cs;
    Traverser traverser;
    LaplaceKernel kernel;
    vec2i tree_depth_range; //[min_depth,max_depth]

protected:
    void P2M();
    void M2M();
    void M2L();
    void L2L();
    void L2P();
    void M2P();
    void P2L();
    void P2P();

protected:
    void check_tree(const Cells3& cells);
    void check_traverser(Traverser& traverser);
    void check_cells(const Cells3& cells);
    void init_cell_matrix(Cells3& cells);
    void init_reg_body(Cells3& cells);
    vec2i get_min_max_depth(const Cells3& cells);
    Indices get_leaf_cell_indices(const Cells3& cells);
    Indices get_nonleaf_cell_indices(const Cells3& cells, int depth);
};

}