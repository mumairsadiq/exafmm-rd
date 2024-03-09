#pragma once
#include "tree.h"
#include <map>

namespace rtfmm
{

class Traverser
{
public:
    Traverser();

    void traverse(Tree& tree);

public:
    std::vector<int> leaf_cell_xmem;
    std::vector<int> nonleaf_cell_xmem;
    std::map<int, std::vector<int>> depth_xmem_map;

    // interaction maps, xmemj <- xmemi
    std::map<int, std::vector<int>> m2l_map;
    std::map<int, std::vector<int>> p2p_map;
    std::map<int, std::vector<int>> m2p_map;
    std::map<int, std::vector<int>> p2l_map;
};

}
