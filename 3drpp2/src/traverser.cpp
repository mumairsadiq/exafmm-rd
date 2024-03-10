#include "traverser.h"
#include "argument.h"
#include "surface.h"

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

void rtfmm::Traverser::traverse(Tree& tree, real cycle, int images, int P)
{
    this->P = P;
    //cells = tree.get_cells();
    tartree = tree.get_cells();
    srctree = tartree;
    if(images == 0)
    {
        horizontal_origin(0,0,0,0);
        make_M2L_parent_map(); // for t
    }
    else if(images >= 1)
    {
        tbegin(horizontal_periodic_near);
        horizontal_periodic_near(cycle);
        tend(horizontal_periodic_near);
        tbegin(make_M2L_parent_map_i1);
        make_M2L_parent_map_i1(cycle); // for t
        tend(make_M2L_parent_map_i1);
        if(images >= 2)
        {
            add_i2_cells(srctree, images);
            add_i2_cells(tartree, images);
            // for t
            tbegin(t_image2);
            int taridx = tartree.size() - 1 - images;
            int srcidx = srctree.size() - 1 - images;
            for(int m = 1; m < images; m++)
            {
                int tarcidx = taridx + m;
                Cell3& ctar = tartree[tarcidx];
                int srccidx = srcidx + m;
                Cell3& csrc = srctree[srccidx];
                for(int pz = -1; pz <= 1; pz++)
                {
                    for(int py = -1; py <= 1; py++)
                    {
                        for(int px = -1; px <= 1; px++)
                        {
                            if(px != 0 || py != 0 || pz != 0)
                            {
                                vec3r offset = vec3r(px,py,pz) * csrc.r * 2;
                                // since c's child is hitorikko, we should mark c->c as image/periodic interaction()
                                M2L_parent_map[tarcidx].push_back(PeriodicParentSource(srccidx, offset , 1)); 
                            }
                        }
                    }
                }
            }
            tend(t_image2);
        }
    }
}

void rtfmm::Traverser::horizontal_origin(int tc, int sc, int tcp, int scp, vec3r offset)
{
    int divide = 0; // 0: nodivide, 1:tc, 2:sc
    if(!adjacent(tartree[tc], srctree[sc], offset))
    {
        if(adjacent(tartree[tc], srctree[scp], offset) && is_leaf(tartree[tc]) && srctree[scp].depth >= tartree[tc].depth)
        {
            if(is_leaf(srctree[sc]) && srctree[sc].brange.number <= get_surface_point_num(P))
            {
                P2P_map[tc].push_back(std::make_pair(sc, offset));
            }
            else
            {
                M2P_pairs.push_back(make_pair(tc, sc, offset));
                std::cout<<"M2P"<<std::endl;
            }
        }
        else if(adjacent(tartree[tcp], srctree[sc], offset) && is_leaf(srctree[sc]) && tartree[tcp].depth >= srctree[sc].depth)
        {
            if(is_leaf(tartree[tc]) && tartree[tc].brange.number <= get_surface_point_num(P))
            {
                P2P_map[tc].push_back(std::make_pair(sc, offset));
            }
            else
            {
                P2L_pairs.push_back(make_pair(tc, sc, offset));
            }
        }
        else if(neighbour(tartree[tcp], srctree[scp], offset))
        {
            M2L_map[tc].push_back(std::make_pair(sc, offset));
        }
        else
        {
            if(is_leaf(tartree[tc]) && is_leaf(srctree[sc]))
            {
                printf("here\n");
            }
            else if(!is_leaf(tartree[tc]) && is_leaf(srctree[sc]))
            {
                divide = 1;
            }
            else if(is_leaf(tartree[tc]) && !is_leaf(srctree[sc]))
            {
                divide = 2;
            }
            else
            {
                divide = tartree[tc].r >= srctree[sc].r ? 1 : 2;
            }
        }
    }
    else
    {
        if(is_leaf(tartree[tc]) && is_leaf(srctree[sc]))
        {
            P2P_map[tc].push_back(std::make_pair(sc, offset));
        }
        else if(!is_leaf(tartree[tc]) && is_leaf(srctree[sc]))
        {
            divide = 1;
        }
        else if(is_leaf(tartree[tc]) && !is_leaf(srctree[sc]))
        {
            divide = 2;
        }
        else
        {
            divide = tartree[tc].r >= srctree[sc].r ? 1 : 2;
        }
    }
    if(divide == 1)
    {
        Cell3 c = tartree[tc];
        for(int i = 0; i < c.crange.number; i++)
        {
            horizontal_origin(c.crange.offset + i, sc, tc, scp, offset);
        }
    }
    else if(divide == 2)
    {
        Cell3 c = srctree[sc];
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
    {
        Cell3 c;
        c.idx = srctree.size();
        c.depth = -1;
        c.octant = 13;
        c.r = srctree[0].r * 3;
        c.x = srctree[0].x;
        c.crange = Range(0, 1);
        c.brange = Range(0, srctree[0].brange.number);
        srctree.push_back(c);
        if(verbose) printf("add srccell %ld, r = %.4f, depth = %d, child = %d\n", srctree.size() - 1, c.r, c.depth, c.crange.offset);
    }
    {
        Cell3 c;
        c.idx = tartree.size();
        c.depth = -1;
        c.octant = 13;
        c.r = tartree[0].r * 3;
        c.x = tartree[0].x;
        c.crange = Range(0, 1);
        c.brange = Range(0, tartree[0].brange.number);
        tartree.push_back(c);
        if(verbose) printf("add tarcell %ld, r = %.4f, depth = %d, child = %d\n", tartree.size() - 1, c.r, c.depth, c.crange.offset);
    }
    int srccidx = srctree.size() - 1;
    int tarcidx = srctree.size() - 1;
    for(int pz = -1; pz <= 1; pz++)
    {
        for(int py = -1; py <= 1; py++)
        {
            for(int px = -1; px <= 1; px++)
            {
                vec3r offset = vec3r(px,py,pz) * cycle;
                horizontal_origin(0, 0, srccidx, tarcidx, offset);
            }   
        }
    }
}

void rtfmm::Traverser::make_M2L_parent_map()
{
    for(int j = 0; j < tartree.size(); j++)
    {
        Cell3& ctar = tartree[j];
        for(int i = 0; i < srctree.size(); i++)
        {
            Cell3& csrc = srctree[i];
            if(i != j && !is_leaf(ctar) && !is_leaf(csrc) && neighbour(csrc,ctar))
            {
                M2L_parent_map[j].push_back(PeriodicParentSource(i, vec3r(0,0,0),0));
            }
        }
    }
}

void rtfmm::Traverser::make_M2L_parent_map_i1(real cycle)
{
    make_M2L_parent_map();
    #pragma omp parallel for
    for(int j = 0; j < tartree.size(); j++)
    {
        Cell3& ctar = tartree[j];
        if(ctar.depth == -1) continue; // the largest cell always neighbour with its 1st periodic images, so we need to exclude it
        for(int i = 0; i < srctree.size(); i++)
        {
            Cell3& csrc = srctree[i];
            for(int pz = -1; pz <= 1; pz++)
            {
                for(int py = -1; py <= 1; py++)
                {
                    for(int px = -1; px <= 1; px++)
                    {
                        if(px != 0 || py != 0 || pz != 0)
                        {
                            vec3r offset = vec3r(px,py,pz) * cycle;
                            if(!is_leaf(ctar) && !is_leaf(csrc) && neighbour(ctar,csrc,offset))
                            {
                                M2L_parent_map[j].push_back(PeriodicParentSource(i, offset, 0));
                            }
                        }
                    }
                }
            }
        }
    }
}

int rtfmm::Traverser::adjacent(Cell3& ca, Cell3& cb, vec3r offset)
{
    vec3r dx = (ca.x - cb.x - offset).abs();
    real dist = (ca.r + cb.r) * 1.001;  //warning : DO NOT ignore the float error
    return dx[0] <= dist && dx[1] <= dist && dx[2] <= dist;
}

int rtfmm::Traverser::neighbour(Cell3& ca, Cell3& cb, vec3r offset)
{
    return adjacent(ca,cb,offset) && ca.depth == cb.depth;
}

int rtfmm::Traverser::is_leaf(Cell3& c)
{
    return c.crange.number == 0;
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
        case OperatorType::P2P : return P2P_map; break;
        default: return PeriodicInteractionMap();
    }
}

rtfmm::PeriodicM2LMap rtfmm::Traverser::get_M2L_parent_map()
{
    return M2L_parent_map;
}

void rtfmm::Traverser::add_i2_cells(Cells3& cells, int images)
{
    /* images >= 2 */
    int child_idx = cells.size() - 1;
    int idx = cells.size() - 1;
    for(int m = 1; m < images; m++)
    {
        Cell3 c;
        c.idx = cells.size();
        c.depth = -1 - m;
        c.octant = 13;
        c.r = cells[child_idx].r * 3;
        c.x = cells[child_idx].x;
        c.crange = Range(child_idx, 1);
        cells.push_back(c);
        child_idx = cells.size() - 1;
        if(verbose) printf("add cell %ld, r = %.4f, depth = %d, child = %d\n", cells.size()-1, c.r, c.depth, c.crange.offset);
    }
}