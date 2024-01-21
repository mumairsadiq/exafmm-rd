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

    void traverse(Tree& tree, real cycle = 2 * M_PI, int images = 0);

    PeriodicInteractionPairs get_pairs(OperatorType type);

    PeriodicInteractionMap get_map(OperatorType type);

    Cells3 get_cells() {return cells;}

private:
    
    /**
     * @brief horizontal traverse
     * @param tc target cell idx
     * @param sc source cell idx
     * @param tcp target cell's parent cell idx
     * @param scp source cell's parent cell idx
    */
    void horizontal_origin(int tc, int sc, int tcp, int scp, vec3r offset = vec3r(0,0,0));

    void horizontal_periodic_near(real cycle);

    void horizontal_periodic_far(real cycle, int image);

    int adjacent(int a, int b, vec3r offset = vec3r(0,0,0));

    int neighbour(int a, int b, vec3r offset = vec3r(0,0,0));

    int is_leaf(int c);

    Cells3 cells;

public:
    PeriodicInteractionPairs P2P_pairs;
    PeriodicInteractionPairs M2L_pairs;
    PeriodicInteractionPairs M2P_pairs;
    PeriodicInteractionPairs P2L_pairs;

    PeriodicInteractionMap M2L_map;
    PeriodicInteractionMap P2P_map;
};

}