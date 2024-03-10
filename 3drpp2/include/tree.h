#pragma once
#include "type.h"
#include "body.h"
#include "tree.h"

namespace rtfmm
{

struct Cell3
{
    int idx;
    /**
     * @brief octant relative to parent cell
     * @note 
     * when [0,7], it means a relative position of 2x2x2 children of the parent cell;
     * when == 13, it means the central position of a 3x3x3 child of the image parent cell, namely, it is hitorikko of its parent(cells whose depth <= 0 are hitorikko).
    */
    int octant;
    int depth;
    real r;
    vec3r x;
    Range crange;
    Range brange;
    Matrix q_equiv;
    Matrix p_check;

    Bodies3 bs;
    Bodies3 rbs;
    std::vector<real> ws;

    friend std::ostream &operator<<(std::ostream & os, const Cell3 & cell) 
    {
        os << "---[cell]--- "
        <<"octant=" << cell.octant
        <<",idx=" << cell.idx
        << ",depth=" << cell.depth
        << ",r=" << cell.r
        <<",center=" << cell.x
        << ",crange=" << cell.crange
        << ",brange=" << cell.brange;
        return os;
    }
};

using Cells3 = std::vector<Cell3>;

class Tree
{
public:
    enum class TreeType
    {
        uniform,
        nonuniform,
        regnonuniform
    };

    Tree();
    
    /**
     * @brief build tree from bodies
     * @warning this function will shuffle the bodies, so DO NOT FORGET to re-sort the bodies in the end
     * @param bodies bodies
     * @param x0 center of box
     * @param r0 half size of box
     * @param m for uniform-tree it stands for max_depth; for nonuniform-tree it stands for max_n_per_cell
     * @param type type of tree (namely, uniform or nonuniform)
     */
    void build(Bodies3& bodies, vec3r x0, real r0, int m, real rega, int images, real cycle, TreeType type);
    
    Cells3 get_cells();

    static vec3r get_child_cell_x(vec3r x_par, real r_par, int octant, int is_periodic);

private:

    void build_uniform_octree(Bodies3& bodies, vec3r x, real r, int max_depth, Cells3& cells);

    void build_nonuniform_octree(Bodies3& bodies, vec3r x, real r, int max_n_per_cell, Cells3& cells);

    void build_reg_nonuniform_octree(Bodies3& bodies, vec3r x, real r, int max_n_per_cell, real rega, Bodies3& c0_rbs, Cells3& cells);

    void build_reg_uniform_octree(vec3r x, real r, int max_depth, vec3i direction, real rega, Bodies3& c0_rbs, Cells3& cells);

    Bodies3 init_c0_reg(Bodies3& bodies, vec3r x0, real r0, real cycle, real rega);

    Bodies3 init_rc0_reg(Bodies3& bodies, vec3r x0, real r0, real cycle, real rega);

private:

    std::vector<Cell3> all_cells;
    Cells3 reg_cells[26];
};

}