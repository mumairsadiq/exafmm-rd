#include "traverser.h"
#include "argument.h"

rtfmm::InteractionPair rtfmm::make_pair(int tar, int src)
{
    return std::make_pair(tar, src);
}

rtfmm::PeriodicInteractionPair rtfmm::make_pair(int tar, int src, vec3r offset)
{
    return std::make_pair(tar, std::make_pair(src, offset));
}

rtfmm::Traverser::Traverser()
{
    
}

void rtfmm::Traverser::traverse(Tree& tree, real cycle, int images)
{
    cells = tree.get_cells();
    if(images == 0)
    {
        if(verbose) printf("horizontal_origin\n");
        horizontal_origin(0,0,0,0);
    }
    else if(images >= 1)
    {
        horizontal_periodic_near(cycle);
        if(images >= 2)
        {
            horizontal_periodic_far(cycle,images);
        }
    }
}

void rtfmm::Traverser::horizontal_origin(int tc, int sc, int tcp, int scp, vec3r offset)
{
    int divide = 0; // 0: nodivide, 1:tc, 2:sc
    if(!adjacent(tc, sc, offset))
    {
        if(adjacent(tc, scp, offset) && is_leaf(tc) && cells[scp].depth >= cells[tc].depth)
        {
            M2P_pairs.push_back(make_pair(tc, sc, offset));
        }
        else if(adjacent(tcp, sc, offset) && is_leaf(sc) && cells[tcp].depth >= cells[sc].depth)
        {
            P2L_pairs.push_back(make_pair(tc, sc, offset));
        }
        else if(neighbour(tcp, scp, offset))
        {
            M2L_pairs.push_back(make_pair(tc, sc, offset));
            M2L_map[tc].push_back(std::make_pair(sc, offset));
        }
        else
        {
            if(is_leaf(tc) && is_leaf(sc))
            {
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
            P2P_pairs.push_back(make_pair(tc, sc, offset));
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
            horizontal_origin(c.crange.offset + i, sc, tc, scp, offset);
        }
    }
    else if(divide == 2)
    {
        Cell3 c = cells[sc];
        for(int i = 0; i < c.crange.number; i++)
        {
            horizontal_origin(tc, c.crange.offset + i, tcp, sc, offset);
        }
    }
}

void rtfmm::Traverser::horizontal_periodic_near(real cycle)
{
    if(verbose) printf("horizontal_periodic_near\n");
    /* images == 1 */
    Cell3 c;
    c.idx = cells.size();
    c.depth = -1;
    c.r = cells[0].r * 3;
    c.x = cells[0].x;
    c.crange = Range(0, 1);
    c.brange = Range(0, cells[0].brange.number);
    cells.push_back(c);
    int cidx = cells.size() - 1;
    if(verbose) printf("add cell %ld, r = %.4f, depth = %d, child = %d\n", cells.size()-1, c.r, c.depth, c.crange.offset);
    for(int pz = -1; pz <= 1; pz++)
    {
        for(int py = -1; py <= 1; py++)
        {
            for(int px = -1; px <= 1; px++)
            {
                horizontal_origin(0,0,cidx,cidx, vec3r(px,py,pz) * cycle);
            }   
        }
    }
}

void rtfmm::Traverser::horizontal_periodic_far(real cycle, int images)
{
    if(verbose) printf("horizontal_periodic_far\n");
    /* images >= 2 */
    int child_idx = cells.size() - 1;
    int idx = cells.size() - 1;
    for(int m = 1; m < images; m++)
    {
        Cell3 c;
        c.idx = cells.size();
        c.depth = -1 - m;
        c.r = cells[child_idx].r * 3;
        c.x = cells[child_idx].x;
        c.crange = Range(child_idx, 1);
        cells.push_back(c);
        child_idx = cells.size() - 1;
        if(verbose) printf("add cell %ld, r = %.4f, depth = %d, child = %d\n", cells.size()-1, c.r, c.depth, c.crange.offset);
    }
    for(int m = 0; m < images - 1; m++)
    {
        int icell_idx = cells[idx].crange.offset;
        if(verbose) printf("icell_idx = %d, cycle = %.4f\n", icell_idx, cycle);
        for(int pz = -1; pz <= 1; pz++)
        {
            for(int py = -1; py <= 1; py++)
            {
                for(int px = -1; px <= 1; px++)
                {
                    if(px == 0 && py == 0 && pz == 0) 
                        continue;
                    for(int cz = -1; cz <= 1; cz++)
                    {
                        for(int cy = -1; cy <= 1; cy++)
                        {
                            for(int cx = -1; cx <= 1; cx++)
                            {
                                int ox = px * 3 + cx;
                                int oy = py * 3 + cy;
                                int oz = pz * 3 + cz;
                                M2L_pairs.push_back(make_pair(0,icell_idx,vec3r(ox,oy,oz) * cycle));
                                M2L_map[0].push_back(std::make_pair(icell_idx, vec3r(ox,oy,oz) * cycle));
                            }
                        }
                    }
                }
            }
        }
        idx++;

        cycle *= 3;
    }
}

int rtfmm::Traverser::adjacent(int a, int b, vec3r offset)
{
    Cell3 ca = cells[a];
    Cell3 cb = cells[b];
    vec3r dx = (ca.x - cb.x - offset).abs();
    real dist = (ca.r + cb.r) * 1.001;  //warning : DO NOT ignore the float error
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

rtfmm::PeriodicInteractionPairs rtfmm::Traverser::get_pairs(OperatorType type)
{
    switch(type)
    {
        case OperatorType::P2P : return P2P_pairs; break;
        case OperatorType::M2L : return M2L_pairs; break;
        case OperatorType::M2P : return M2P_pairs; break;
        case OperatorType::P2L : return P2L_pairs; break;
        default: return PeriodicInteractionPairs();
    }
}

rtfmm::PeriodicInteractionMap rtfmm::Traverser::get_map(OperatorType type)
{
    switch(type)
    {
        case OperatorType::M2L : return M2L_map; break;
        default: return PeriodicInteractionMap();
    }
}