#pragma once
#include "tree.h"

namespace rtfmm
{

using InteractionPair = std::pair<int, int>;
using InteractionPairs = std::vector<InteractionPair>;

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

    void traverse(Tree& tree);

    InteractionPairs get_pairs(OperatorType type);

private:
    
    /**
     * @brief horizontal traverse
     * @param tc target cell idx
     * @param sc source cell idx
     * @param tcp target cell's parent cell idx
     * @param scp source cell's parent cell idx
    */
    void horizontal(int tc, int sc, int tcp, int scp);

    int adjacent(int a, int b, vec3r offset = vec3r(0,0,0));

    int neighbour(int a, int b, vec3r offset = vec3r(0,0,0));

    int is_leaf(int c);

    Cells3 cells;

public:
    InteractionPairs P2P_pairs;
    InteractionPairs M2L_pairs;
    InteractionPairs M2P_pairs;
    InteractionPairs P2L_pairs;
};

}