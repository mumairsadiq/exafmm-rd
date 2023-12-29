#include "traverser.h"

rtfmm::Traverser::Traverser()
{
    
}

void rtfmm::Traverser::traverse(Tree& tree)
{
    cells = tree.get_cells();
    horizontal(0,0,0,0);
}

void rtfmm::Traverser::horizontal(int tc, int sc, int tcp, int scp)
{
    int divide = 0; // 0: nodivide, 1:tc, 2:sc
    if(!adjacent(tc, sc))
    {
        if(adjacent(tc, scp) && is_leaf(tc) && cells[scp].depth >= cells[tc].depth)
        {
            M2P_pairs.push_back(std::make_pair(tc, sc));
        }
        else if(adjacent(tcp, sc) && is_leaf(sc) && cells[tcp].depth >= cells[sc].depth)
        {
            P2L_pairs.push_back(std::make_pair(tc, sc));
        }
        else if(neighbour(tcp, scp))
        {
            M2L_pairs.push_back(std::make_pair(tc, sc));
        }
        else
        {
            if(is_leaf(tc) && is_leaf(sc))
            {
                //M2L_pairs.push_back(std::make_pair(tc, sc));
                printf("here\n");
            }
            else if(!is_leaf(tc) && is_leaf(sc))
            {
                divide = 1;
            }
            else if(is_leaf(tc) && !is_leaf(sc))
            {
                divide = 2;
            }
            else
            {
                divide = cells[tc].r >= cells[sc].r ? 1 : 2;
            }
        }
    }
    else
    {
        if(is_leaf(tc) && is_leaf(sc))
        {
            P2P_pairs.push_back(std::make_pair(tc, sc));
        }
        else if(!is_leaf(tc) && is_leaf(sc))
        {
            divide = 1;
        }
        else if(is_leaf(tc) && !is_leaf(sc))
        {
            divide = 2;
        }
        else
        {
            divide = cells[tc].r >= cells[sc].r ? 1 : 2;
        }
    }
    if(divide == 1)
    {
        Cell3 c = cells[tc];
        for(int i = 0; i < c.crange.number; i++)
        {
            horizontal(c.crange.offset + i, sc, tc, scp);
        }
    }
    else if(divide == 2)
    {
        Cell3 c = cells[sc];
        for(int i = 0; i < c.crange.number; i++)
        {
            horizontal(tc, c.crange.offset + i, tcp, sc);
        }
    }
}

int rtfmm::Traverser::adjacent(int a, int b, vec3r offset)
{
    Cell3 ca = cells[a];
    Cell3 cb = cells[b];
    vec3r dx = (ca.x - cb.x - offset).abs();
    real dist = ca.r + cb.r;
    return dx[0] <= dist && dx[1] <= dist && dx[2] <= dist;
}

int rtfmm::Traverser::neighbour(int a, int b, vec3r offset)
{
    return adjacent(a,b,offset) && cells[a].depth == cells[b].depth;
}

int rtfmm::Traverser::is_leaf(int c)
{
    return cells[c].crange.number == 0;
}

rtfmm::InteractionPairs rtfmm::Traverser::get_pairs(OperatorType type)
{
    switch(type)
    {
        case OperatorType::P2P : return P2P_pairs; break;
        case OperatorType::M2L : return M2L_pairs; break;
        case OperatorType::M2P : return M2P_pairs; break;
        case OperatorType::P2L : return P2L_pairs; break;
    }
}