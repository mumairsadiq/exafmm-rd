#include "tree.h"
#include <queue>
#include "mathfunc.h"
#include <algorithm>

rtfmm::Tree::Tree(Bodies& bs, const Argument& args)
{
    r0 = args.r0;
    x0 = args.x0;
    ncrit = args.ncrit;
    images = args.images;

    build_octree(bs);   // nonreg part
    build_base27tree(args.images); // reg part
    precompute_attributes();
}

void rtfmm::Tree::build_octree(Bodies& bs)
{
    if(verbose) std::cout<<"build adaptive octree"<<std::endl;

    int num_body = bs.size();
    Cell root;
    root.xlogic = vec4i(0,0,0,0);
    root.xphys = logic2phys(root.xlogic, r0, x0);
    root.xmem = 0;
    root.r = r0;
    root.key = root.xlogic.hash(images);
    root.qidx = -1;
    root.pidx = -1;
    root.brange = vec2i(0, num_body);
    root.is_leaf = 1;
    this->cs.push_back(root);

    std::queue<int> big_cells;
    big_cells.push(0);

    while(!big_cells.empty())
    {
        int branch_idx = big_cells.front();
        big_cells.pop();
        Cell branch = this->cs[branch_idx];
        
        if(branch.brange[1] > ncrit)
        {
            this->cs[branch_idx].is_leaf = 0;
            int begin = branch.brange[0];
            int num = branch.brange[1];
            int end = branch.brange[0] + branch.brange[1] - 1;
            int split_num = 8;
            //count
            std::vector<int> quad_num(split_num);
            for(int i = begin; i <= end; i++)
            {
                int idx = get_body_octant(bs[i].x, branch.xphys, branch.r);
                quad_num[idx]++;
            }
            std::vector<int> offset = exclusive_scan(quad_num);
            //sort
            std::vector<Body> bodies_buffer;
            bodies_buffer.resize(num);
            for(int i = begin; i <= end; i++)
            {
                int idx = get_body_octant(bs[i].x, branch.xphys, branch.r);
                bodies_buffer[offset[idx]] = bs[i];
                offset[idx]++;
            }
            for(int i = 0; i < num; i++)
            {
                bs[begin + i] = bodies_buffer[i];
            }
            //create leaves
            for(int octant = 0; octant < split_num; octant++)
            {
                Cell child;
                child.xlogic = branch.xlogic.child(octant);
                child.xphys = logic2phys(child.xlogic, r0, x0);
                child.xmem = this->cs.size();
                child.r = branch.r / 2;
                child.key = child.xlogic.hash(images);
                child.pidx = -1;
                child.qidx = -1;
                child.brange = vec2i(begin + offset[octant] - quad_num[octant], quad_num[octant]);
                child.is_leaf = 1;
                this->cs.push_back(child);
                big_cells.push(child.xmem);
            }
        }
    }
    for(int i = 0; i < this->cs.size(); i++)
    {
        Cell& c = this->cs[i];
        cidx_map[c.key] = c.xmem;
    }
}

void rtfmm::Tree::build_base27tree(int images)
{
    for(int depth = -1; depth >= -images; depth--)
    {
        real r = r0 * std::pow(3, -depth - 1);
        for(int z = 0; z < 3; z++)
        {
            for(int y = 0; y < 3; y++)
            {
                for(int x = 0; x < 3; x++)
                {
                    Cell cell;
                    cell.xlogic = vec4i(depth,x,y,z);
                    cell.xphys = logic2phys(cell.xlogic, r0, x0);
                    cell.xmem = this->cs.size();
                    cell.r = r;
                    cell.key = cell.xlogic.hash(images);
                    cell.pidx = -1;
                    cell.qidx = -1;
                    cell.brange = 0;
                    cell.is_leaf = 0;
                    this->cs.push_back(cell);
                    cidx_map[cell.key] = cell.xmem;
                }
            }
        }
        //this->cs[find_cidx(LogicCoord(depth, 1, 1, 1))];
    }
}

void rtfmm::Tree::precompute_attributes()
{
    vec2i res = {cs[0].xlogic[0], cs[0].xlogic[0]};
    for(int i = 0; i < cs.size(); i++)
    {
        res[0] = std::min(res[0], cs[i].xlogic[0]);
        res[1] = std::max(res[1], cs[i].xlogic[0]);
    }
    depth_range = res;
}