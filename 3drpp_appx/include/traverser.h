#pragma once
#include "tree.h"
#include <map>

namespace rtfmm
{

using InteractionPair = std::pair<int, int>;
using PeriodicInteractionPair = std::pair<int, std::pair<int, vec3r>>;
using InteractionPairs = std::vector<InteractionPair>;
using PeriodicInteractionPairs = std::vector<PeriodicInteractionPair>;

using PeriodicInteractionMap = std::map<int, std::vector<std::pair<int, vec3r>>>;


struct PeriodicParentSource
{
    PeriodicParentSource(int idx_, vec3r offset_, int is_single_parent_) : idx(idx_), offset(offset_), is_single_parent(is_single_parent_){}
    int idx;
    vec3r offset;
    int is_single_parent; // is the cell has only one child(namely, is the cell a image cell)
};

/**
 * @brief map to store M2L parent pair
*/
using PeriodicM2LMap = std::map<int, std::vector<PeriodicParentSource>>;

InteractionPair make_pair(int tar, int src);
PeriodicInteractionPair make_pair(int tar, int src, vec3r offset);

enum class OperatorType
{
    P2P,
    M2L,
    M2P,
    P2L
};

class Traverser
{
    
public:
    Traverser();

    void traverse(Tree& tree, real cycle, int images, int P_);

    PeriodicInteractionPairs get_pairs(OperatorType type);

    PeriodicInteractionMap get_map(OperatorType type);

    PeriodicM2LMap get_M2L_parent_map();

    Cells3 get_cells() {return cells;}

    PeriodicInteractionMap get_m2l_map_from_m2l_parent_map();

private:
    
    /**
     * @brief horizontal traverse
     * @param tc target cell idx
     * @param sc source cell idx
     * @param tcp target cell's parent cell idx
     * @param scp source cell's parent cell idx
    */
    void horizontal_origin(int tc, int sc, int tcp, int scp, vec3r offset);

    void horizontal_periodic_near(real cycle);

    void horizontal_periodic_far(real cycle, int image);

    /**
     * @brief if b+offset adjacent with a
    */
    int adjacent(int a, int b, vec3r offset = vec3r(0,0,0));

    /**
     * @brief if b+offset neighbour with a
    */
    int neighbour(int a, int b, vec3r offset = vec3r(0,0,0));

    int is_leaf(int c);

    void make_M2L_parent_map();

    void make_M2L_parent_map_i1(real cycle);

    Cells3 cells;

    int P;

public:
    PeriodicInteractionPairs P2P_pairs;
    PeriodicInteractionPairs M2L_pairs;
    PeriodicInteractionPairs M2P_pairs;
    PeriodicInteractionPairs P2L_pairs;

    PeriodicInteractionMap M2L_map;
    PeriodicInteractionMap P2P_map;

    //PeriodicInteractionMap M2L_parent_map;
    PeriodicM2LMap M2L_parent_map;

    std::vector<int> leaf_cell_idx;
};

}