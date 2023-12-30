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
    void check_tree(const Cells3& cells);
    void check_traverser(Traverser& traverser);
    void init_cell_matrix(Cells3& cells);
    int get_min_depth(const Cells3& cells);
    int get_max_depth(const Cells3& cells);
    Indices get_leaf_cell_indices(const Cells3& cells);
    Indices get_nonleaf_cell_indices(const Cells3& cells, int depth);
};

}