#include "traverser.h"
#include "argument.h"
#include <set>

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
        make_M2L_parent_map(); // for t
        for(int j = 0; j < cells[0].crange.number; j++)
        {
            int jidx = cells[0].crange.offset + j;
            for(int i = 0; i < cells[0].crange.number; i++)
            {
                if(i != j)
                {
                    int iidx = cells[0].crange.offset + i;
                    if(!adjacent(jidx, iidx))
                    {
                        printf("%d <- %d\n", jidx, iidx);
                        M2L_map[jidx].push_back(std::make_pair(iidx, vec3r(0,0,0)));
                    }
                }
            }
        }
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

    build_tree_depth_range();

    build_leaf_cell_idxs();

    build_nonleaf_cell_idxs();

    build_rw_body_cell_idxs();
}

void rtfmm::Traverser::horizontal_origin(int tc, int sc, int tcp, int scp, vec3r offset)
{
    //std::cout<<cells[tc]<<"   "<<cells[sc]<<std::endl;
    int divide = 0; // 0: nodivide, 1:tc, 2:sc
    if(!adjacent(tc, sc, offset))
    {
        if(adjacent(tc, scp, offset) && is_leaf(tc) && cells[scp].depth >= cells[tc].depth)
        {
            //M2P_pairs.push_back(make_pair(tc, sc, offset));
            P2P_pairs.push_back(make_pair(tc, sc, offset));
            P2P_map[tc].push_back(std::make_pair(sc, offset));
        }
        else if(adjacent(tcp, sc, offset) && is_leaf(sc) && cells[tcp].depth >= cells[sc].depth)
        {
            //P2L_pairs.push_back(make_pair(tc, sc, offset));
            P2P_pairs.push_back(make_pair(tc, sc, offset));
            P2P_map[tc].push_back(std::make_pair(sc, offset));
        }
        else if(neighbour(tcp, scp, offset))
        {
            M2L_pairs.push_back(make_pair(tc, sc, offset));
            //M2L_map[tc].push_back(std::make_pair(sc, offset));
        }
        else
        {
            if(is_leaf(tc) && is_leaf(sc))
            {
                printf("here\n");
                std::cout<<cells[tc]<<"   "<<cells[sc]<<std::endl;
                std::cout<<cells[tcp]<<"   "<<cells[scp]<<std::endl;
                exit(0);
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
            P2P_map[tc].push_back(std::make_pair(sc, offset));
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
                                //M2L_map[icell_idx].push_back(std::make_pair(icell_idx, vec3r(ox,oy,oz) * cycle));
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
            if(i != j && !is_leaf(j) && !is_leaf(i) && neighbour(i,j))
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
                            if(!is_leaf(j) && !is_leaf(i) && neighbour(j,i,offset))
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
                    if(!adjacent(cjidx, ciidx, offset))
                    {
                        res[cjidx].push_back(std::make_pair(ciidx, offset));
                    }
                }
            }
        }
    }
    return res;
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
        case OperatorType::P2P : return P2P_map; break;
        default: return PeriodicInteractionMap();
    }
}

rtfmm::PeriodicM2LMap rtfmm::Traverser::get_M2L_parent_map()
{
    return M2L_parent_map;
}

void rtfmm::Traverser::build_leaf_cell_idxs()
{
    for(int i = 0; i < cells.size(); i++)
    {
        if(cells[i].crange.number == 0)
        {
            leaf_cell_idxs.push_back(i);
        }
    }
}

void rtfmm::Traverser::build_nonleaf_cell_idxs()
{
    for(int depth = tree_depth_range[0]; depth <= tree_depth_range[1]; depth++)
    {
        Indices res;
        for(int i = 0; i < cells.size(); i++)
        {
            if(cells[i].depth == depth && cells[i].crange.number > 0)
            {
                res.push_back(i);
            }
        }
        nonleaf_cell_idxs[depth] = res;
    }
}

void rtfmm::Traverser::build_tree_depth_range()
{
    vec2i res = {cells[0].depth, cells[0].depth};
    for(int i = 0; i < cells.size(); i++)
    {
        res[0] = std::min(res[0], cells[i].depth);
        res[1] = std::max(res[1], cells[i].depth);
    }
    tree_depth_range = res;
}

void rtfmm::Traverser::build_rw_body_cell_idxs()
{
    rw_body_cell_idxs.clear();
    std::set<int> idx_set;
    for(auto idx : leaf_cell_idxs)
    {
        idx_set.insert(idx);
    }
    for(auto p2p : P2P_map)
    {
        idx_set.insert(p2p.first);
        for(int i = 0; i < p2p.second.size(); i++)
        {
            idx_set.insert(p2p.second[i].first);
        }
    }
    for(auto idx : idx_set)
    {
        rw_body_cell_idxs.push_back(idx);
    }
}