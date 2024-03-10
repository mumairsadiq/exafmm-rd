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
    cells = tree.get_cells();
    if(images == 0)
    {
        if(verbose) printf("horizontal_origin\n");
        horizontal_origin(0,0,0,0);
        make_M2L_parent_map(); // for t
    }
    else if(images >= 1)
    {
        tbegin(horizontal_periodic_near);
        horizontal_periodic_near(cycle, cells, cells);
        tend(horizontal_periodic_near);
        tbegin(make_M2L_parent_map_i1);
        make_M2L_parent_map_i1(cycle); // for t
        tend(make_M2L_parent_map_i1);
        if(images >= 2)
        {
            tbegin(horizontal_periodic_far);
            horizontal_periodic_far(cycle,images);
            tend(horizontal_periodic_far);

            // for t
            tbegin(t_image2);
            int idx = cells.size() - 1 - images;
            for(int m = 1; m < images; m++)
            {
                int cidx = idx + m;
                Cell3& c = cells[cidx];
                for(int pz = -1; pz <= 1; pz++)
                {
                    for(int py = -1; py <= 1; py++)
                    {
                        for(int px = -1; px <= 1; px++)
                        {
                            if(px != 0 || py != 0 || pz != 0)
                            {
                                vec3r offset = vec3r(px,py,pz) * c.r * 2;
                                // since c's child is hitorikko, we should mark c->c as image/periodic interaction()
                                M2L_parent_map[cidx].push_back(PeriodicParentSource(cidx, offset , 1)); 
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
    if(!adjacent(cells[tc], cells[sc], offset))
    {
        if(adjacent(cells[tc], cells[scp], offset) && is_leaf(cells[tc]) && cells[scp].depth >= cells[tc].depth)
        {
            if(is_leaf(cells[sc]) && cells[sc].brange.number <= get_surface_point_num(P))
            {
                P2P_map[tc].push_back(std::make_pair(sc, offset));
            }
            else
            {
                M2P_pairs.push_back(make_pair(tc, sc, offset));
                std::cout<<"M2P"<<std::endl;
            }
        }
        else if(adjacent(cells[tcp], cells[sc], offset) && is_leaf(cells[sc]) && cells[tcp].depth >= cells[sc].depth)
        {
            if(is_leaf(cells[tc]) && cells[tc].brange.number <= get_surface_point_num(P))
            {
                P2P_map[tc].push_back(std::make_pair(sc, offset));
            }
            else
            {
                P2L_pairs.push_back(make_pair(tc, sc, offset));
            }
        }
        else if(neighbour(cells[tcp], cells[scp], offset))
        {
            M2L_map[tc].push_back(std::make_pair(sc, offset));
        }
        else
        {
            if(is_leaf(cells[tc]) && is_leaf(cells[sc]))
            {
                printf("here\n");
            }
            else if(!is_leaf(cells[tc]) && is_leaf(cells[sc]))
            {
                divide = 1;
            }
            else if(is_leaf(cells[tc]) && !is_leaf(cells[sc]))
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
        if(is_leaf(cells[tc]) && is_leaf(cells[sc]))
        {
            P2P_map[tc].push_back(std::make_pair(sc, offset));
        }
        else if(!is_leaf(cells[tc]) && is_leaf(cells[sc]))
        {
            divide = 1;
        }
        else if(is_leaf(cells[tc]) && !is_leaf(cells[sc]))
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

void rtfmm::Traverser::horizontal_periodic_near(real cycle, Cells3& srcs, Cells3& tars)
{
    if(verbose) printf("horizontal_periodic_near\n");
    /* images == 1 */
    Cell3 c;
    c.idx = cells.size();
    c.depth = -1;
    c.octant = 13;
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
                vec3r offset = vec3r(px,py,pz) * cycle;
                horizontal_origin(0,0,cidx,cidx, offset);
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
        c.octant = 13;
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
                                //M2L_pairs.push_back(make_pair(0,icell_idx,vec3r(ox,oy,oz) * cycle));
                                //M2L_map[0].push_back(std::make_pair(icell_idx, vec3r(ox,oy,oz) * cycle));
                                M2L_pairs.push_back(make_pair(icell_idx,icell_idx,vec3r(ox,oy,oz) * cycle));
                                M2L_map[icell_idx].push_back(std::make_pair(icell_idx, vec3r(ox,oy,oz) * cycle));
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

void rtfmm::Traverser::make_M2L_parent_map()
{
    for(int j = 0; j < cells.size(); j++)
    {
        Cell3& ctar = cells[j];
        for(int i = 0; i < cells.size(); i++)
        {
            Cell3& csrc = cells[i];
            if(i != j && !is_leaf(cells[j]) && !is_leaf(cells[i]) && neighbour(cells[i],cells[j]))
            {
                //M2L_parent_map[j].push_back(std::make_pair(i, vec3r(0,0,0)));
                M2L_parent_map[j].push_back(PeriodicParentSource(i, vec3r(0,0,0),0));
            }
        }
    }
}

void rtfmm::Traverser::make_M2L_parent_map_i1(real cycle)
{
    make_M2L_parent_map();
    #pragma omp parallel for
    for(int j = 0; j < cells.size(); j++)
    {
        Cell3& ctar = cells[j];
        if(ctar.depth == -1) continue; // the largest cell always neighbour with its 1st periodic images, so we need to exclude it
        for(int i = 0; i < cells.size(); i++)
        {
            Cell3& csrc = cells[i];
            for(int pz = -1; pz <= 1; pz++)
            {
                for(int py = -1; py <= 1; py++)
                {
                    for(int px = -1; px <= 1; px++)
                    {
                        if(px != 0 || py != 0 || pz != 0)
                        {
                            vec3r offset = vec3r(px,py,pz) * cycle;
                            if(!is_leaf(cells[j]) && !is_leaf(cells[i]) && neighbour(cells[j],cells[i],offset))
                            {
                                //M2L_parent_map[j].push_back(std::make_pair(i, offset));
                                M2L_parent_map[j].push_back(PeriodicParentSource(i, offset, 0));
                            }
                        }
                    }
                }
            }
        }
    }
}

rtfmm::PeriodicInteractionMap rtfmm::Traverser::get_m2l_map_from_m2l_parent_map()
{
    std::cout<<M2L_parent_map.size()<<std::endl;
    PeriodicInteractionMap res;
    for(auto m2l_parent : M2L_parent_map)
    {
        Cell3& ctar = cells[m2l_parent.first];
        for(int i = 0; i < m2l_parent.second.size(); i++)
        {
            Cell3& csrc = cells[m2l_parent.second[i].idx];
            vec3r offset = m2l_parent.second[i].offset;
            for(int cj = 0; cj < ctar.crange.number; cj++)
            {
                int cjidx = ctar.crange.offset + cj;
                Cell3& ctarc = cells[cjidx];
                for(int ci = 0; ci < csrc.crange.number; ci++)
                {
                    int ciidx = csrc.crange.offset + ci;
                    Cell3& csrcc = cells[ciidx];
                    if(!adjacent(cells[cjidx], cells[ciidx], offset))
                    {
                        res[cjidx].push_back(std::make_pair(ciidx, offset));
                    }
                }
            }
        }
    }
    return res;
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