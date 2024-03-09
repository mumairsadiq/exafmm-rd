#include "traverser.h"
#include <queue>
#include <set>

rtfmm::Traverser::Traverser()
{

}

static std::vector<rtfmm::LogicCoord> get_legal_cousins(rtfmm::LogicCoord& xlogic)
{
    std::vector<rtfmm::LogicCoord> res;
    rtfmm::LogicCoord parent = xlogic.parent();
    for(int z = -1; z <= 1; z++)
    {
        for(int y = -1; y <= 1; y++)
        {
            for(int x = -1; x <= 1; x++)
            {
                rtfmm::LogicCoord uncle = parent + rtfmm::LogicCoord(0, x, y, z);
                for(int octant = 0; octant < 8; octant++)
                {
                    rtfmm::LogicCoord cousin = uncle.child(octant);
                    if(cousin.legal())
                    {
                        res.push_back(cousin);
                    }
                }
            }
        }
    }
    return res;
}

void rtfmm::Traverser::traverse(Tree& tree)
{
    std::cout<<"traverse"<<std::endl;

    std::map<int, std::set<int>> m2l_set;
    std::map<int, std::set<int>> p2p_set;
    std::map<int, std::set<int>> m2p_set;
    std::map<int, std::set<int>> p2l_set;

    // leaf/nonleaf/depth
    for(int j = 0; j < tree.cs.size(); j++)
    {
        Cell& cj = tree.cs[j];
        if(cj.is_leaf)
        {
            leaf_cell_xmem.push_back(cj.xmem);
        }
        else
        {
            nonleaf_cell_xmem.push_back(cj.xmem);
        }
        depth_xmem_map[cj.xlogic[0]].push_back(cj.xmem);
    }

    // M2L map
    for(int j = 0; j < nonleaf_cell_xmem.size(); j++)
    {
        int xmem_j = nonleaf_cell_xmem[j];
        Cell& cj = tree.cs[xmem_j];
        if(!cj.is_leaf)
        {
            for(int z = -1; z <= 1; z++)
            {
                for(int y = -1; y <= 1; y++)
                {
                    for(int x = -1; x <= 1; x++)
                    {
                        if(x != 0 || y != 0 || z != 0)
                        {
                            LogicCoord xlogic_i = cj.xlogic + vec4i(0, x, y, z);
                            int xmem_i = tree.find_cidx(xlogic_i);
                            if(xmem_i != -1 && !tree.cs[xmem_i].is_leaf)
                            {
                                m2l_map[cj.xmem].push_back(xmem_i);
                                assert_exit(xmem_j == cj.xmem, "xmem error!");
                            }
                        }
                    }
                }
            }
        }
    }

    // other maps
    for(int j = 0; j < tree.cs.size(); j++)
    {
        Cell& cj = tree.cs[j];
        std::vector<LogicCoord> cousins = get_legal_cousins(cj.xlogic);
        for(int i = 0; i < cousins.size(); i++)
        {
            LogicCoord& xi = cousins[i];
            int xmem_i = tree.find_cidx(xi);
            if(xmem_i != -1) // ci cj same depth
            {
                Cell& ci = tree.cs[xmem_i];
                if(cj.is_leaf)
                {
                    if(ci.is_leaf)
                    {
                        if(LogicCoord::adjacent(cj.xlogic, ci.xlogic))
                        {
                            p2p_set[cj.xmem].insert(ci.xmem);
                        }
                        else
                        {
                            // nothing, delt by m2l
                        }
                    }
                    else
                    {
                        if(LogicCoord::adjacent(cj.xlogic, ci.xlogic))
                        {
                            std::queue<LogicCoord> big_cell;
                            big_cell.push(ci.xlogic);
                            while(!big_cell.empty())
                            {
                                LogicCoord xlogicii = big_cell.front();
                                big_cell.pop();
                                int xmemii = tree.find_cidx(xlogicii);
                                if(xmemii != -1)
                                {
                                    Cell& cii = tree.cs[xmemii];
                                    if(LogicCoord::adjacent(cj.xlogic, cii.xlogic))
                                    {
                                        if(cii.is_leaf)
                                        {
                                            p2p_set[cj.xmem].insert(cii.xmem);
                                        }
                                        else
                                        {
                                            for(int octant = 0; octant < 8; octant++)
                                            {
                                                LogicCoord child = ci.xlogic.child(octant);
                                                big_cell.push(child);
                                            }
                                        }
                                    }
                                    else
                                    {
                                        m2p_set[cj.xmem].insert(cii.xmem);
                                    }
                                }
                            }
                        }
                        else
                        {
                            // nothing, delt by m2l
                        }
                    }
                }
            }
            else // find leaf parent of i
            {
                do
                {
                    xi = xi.parent();
                    xmem_i = tree.find_cidx(xi);
                }while(xmem_i == -1);
                Cell& ci = tree.cs[xmem_i]; // ci is larger than cj, and ci is leaf
                if(cj.is_leaf)
                {
                    if(LogicCoord::adjacent(cj.xlogic, ci.xlogic))
                    {
                        p2p_set[cj.xmem].insert(ci.xmem);
                    }
                    else
                    {
                        p2l_set[cj.xmem].insert(ci.xmem);
                    }
                }
                else
                {
                    if(LogicCoord::adjacent(cj.xlogic, ci.xlogic))
                    {
                        // nothing, delt by other loop
                    }
                    else
                    {
                        p2l_set[cj.xmem].insert(ci.xmem);
                    }
                }
            }
        }
    }

    for(auto p2p : p2p_set)
    {
        p2p_map[p2p.first] = std::vector<int>(p2p.second.begin(), p2p.second.end());
    }
    for(auto p2l : p2l_set)
    {
        p2l_map[p2l.first] = std::vector<int>(p2l.second.begin(), p2l.second.end());
    }
    for(auto m2p : m2p_set)
    {
        m2p_map[m2p.first] = std::vector<int>(m2p.second.begin(), m2p.second.end());
    }
}