#pragma once
#include "type.h"
#include "body.h"
#include "tree.h"

namespace rtfmm
{

struct Cell3
{
    int idx;
    int depth;
    real r;
    vec3r x;
    Range crange;
    Range brange;
    matrix q_equiv;
    matrix p_check;
};

using Cells3 = std::vector<Cell3>;

class Tree
{
public:
    enum class TreeType
    {
        uniform,
        nonuniform
    };

    Tree();
    
    /**
     * @brief build tree from bodies
     * @warning this function will shuffle the bodies, so DON'T FORGET to re-sort the bodies in the end
     * @param bodies bodies
     * @param x center of box
     * @param r half size of box
     * @param m for uniform-tree it stands for max_depth; for nonuniform-tree it stands for max_n_per_cell
     * @param type type of tree (namely, uniform or nonuniform)
     */
    void build(Bodies3& bodies, vec3r x, real r, int m, TreeType type);
    
    Cells3 get_cells();

private:

    void build_uniform_octree(Bodies3& bodies, vec3r x, real r, int max_depth);

    void build_nonuniform_octree(Bodies3& bodies, vec3r x, real r, int max_n_per_cell);

private:

    std::vector<Cell3> cells;
};

}