#include "tree.h"
#include "argument.h"
#include <queue>
#include "kernel.h"

rtfmm::Tree::Tree()
{

}

rtfmm::Bodies3 rtfmm::Tree::init_c0_reg(Bodies3& bodies, vec3r x0, real r0, real cycle, real rega)
{
    Bodies3 res;
    for(int z = -1; z <= 1; z++)
    {
        for(int y = -1; y <= 1; y++)
        {
            for(int x = -1; x <= 1; x++)
            {
                if(x != 0 || y != 0 || z != 0)
                {
                    vec3r offset = vec3r(x,y,z) * cycle;
                    for(int i = 0; i < bodies.size(); i++)
                    {
                        Body3 b = bodies[i];
                        real w = LaplaceKernel::get_rega_w(b.x + offset - x0, r0, rega);
                        if(w > 0)
                        {
                            b.idx = -1; // only for M
                            b.x += offset;
                            res.push_back(b);
                        }
                    }
                }
            }
        }
    }
    return res;
}

rtfmm::Bodies3 rtfmm::Tree::init_rc0_reg(Bodies3& bodies, vec3r x0, real r0, real cycle, real rega)
{
    Bodies3 res;
    for(int i = 0; i < bodies.size(); i++)
    {
        Body3 b = bodies[i];
        real w = LaplaceKernel::get_rega_w(b.x - x0, r0, rega);
        if(w > 0)
        {
            std::cout<<w<<" "<<b<<std::endl; 
            res.push_back(b);
        }
    }
    return res;
}

void rtfmm::Tree::build(Bodies3& bodies, vec3r x0, real r0, int m, real rega, int images, real cycle, TreeType type)
{
    if(type == TreeType::uniform)
    {
        build_uniform_octree(bodies, x0, r0, m, this->all_cells);
    }
    else if(type == TreeType::nonuniform)
    {
        build_nonuniform_octree(bodies, x0, r0, m, this->all_cells);
    }
    else if(type == TreeType::regnonuniform)
    {
        Bodies3 rbs0;
        if(images > 0)
        {
            rbs0 = init_c0_reg(bodies, x0, r0, cycle, rega);
        }
        build_reg_nonuniform_octree(bodies, x0, r0, m, rega, rbs0, this->all_cells);
        int max_depth = 0;
        for(int i = 0; i < this->all_cells.size(); i++)
        {
            max_depth = std::max(max_depth, this->all_cells[i].depth);
        }
        std::cout<<"bodies.size = "<<bodies.size()<<std::endl;
        std::cout<<"max_depth = "<<max_depth<<std::endl;
        int cnt = 0;
        for(int z = -1; z <= 1; z++)
        {
            for(int y = -1; y <= 1; y++)
            {
                for(int x = -1; x <= 1; x++)
                {
                    if(x != 0 || y != 0 || z != 0)
                    {
                        vec3i direction(x == 0 ? -1 : (-x + 1) / 2, y == 0 ? -1 : (-y + 1) / 2, z == 0 ? -1 : (-z + 1) / 2);
                        vec3r center = x0 + vec3r(x,y,z) * cycle;
                        Bodies3 rc0_rbs = init_rc0_reg(bodies, center, r0, cycle, rega);
                        build_reg_uniform_octree(center, r0, max_depth, direction, rega, rc0_rbs, reg_cells[cnt]);
                        cnt++;
                    }
                }
            }
        }
    }
}

void rtfmm::Tree::build_uniform_octree(Bodies3& bodies, vec3r x, real r, int max_depth, Cells3& cells)
{
    if(verbose) std::cout<<"build complete balanced octree"<<std::endl;
    std::queue<int> big_cells;
    int num_body = bodies.size();
    Cell3 root;
    root.idx = 0;
    root.octant = -1;
    root.depth = 0;
    root.r = r;
    root.x = x;
    root.crange = Range(0,0);
    root.brange = Range(0, bodies.size());
    cells.push_back(root);
    big_cells.push(0);

    while(!big_cells.empty())
    {
        int branch_cell_idx = big_cells.front();
        Cell3 branch_cell = cells[branch_cell_idx];
        big_cells.pop();

        if(branch_cell.depth < max_depth)
        {
            int begin = branch_cell.brange.offset;
            int num = branch_cell.brange.number;
            int end = branch_cell.brange.offset + branch_cell.brange.number - 1;

            //divide
            int quad_num[8] = {0,0,0,0,0,0,0,0};
            int offset[8] = {0,0,0,0,0,0,0,0};
            int offset_sum = 0;
            for(int i = begin; i <= end; i++)
            {
                int idx = ((bodies[i].x[0] > branch_cell.x[0]) << 2) 
                        + ((bodies[i].x[1] > branch_cell.x[1]) << 1) 
                        + ((bodies[i].x[2] > branch_cell.x[2]) << 0);
                quad_num[idx]++;
            }
            for(int i = 0; i < 8; i++)
            {
                offset[i] = offset_sum;
                offset_sum += quad_num[i];
            }
            //sort
            std::vector<Body3> bodies_buffer;
            bodies_buffer.resize(num);
            for(int i = begin; i <= end; i++)
            {
                int idx = ((bodies[i].x[0] > branch_cell.x[0]) << 2) 
                        + ((bodies[i].x[1] > branch_cell.x[1]) << 1) 
                        + ((bodies[i].x[2] > branch_cell.x[2]) << 0);
                bodies_buffer[offset[idx]] = bodies[i];
                offset[idx]++;
            }
            //store
            for(int i = 0; i < num; i++)
            {
                bodies[begin + i] = bodies_buffer[i];
            }
            //create leaves
            int insert_offset = cells.size();
            cells[branch_cell_idx].crange = Range(insert_offset,8);
            for(int i = 0; i < 8; i++)
            {
                Cell3 child_cell;
                child_cell.idx = cells.size();
                child_cell.octant = i;
                child_cell.depth = branch_cell.depth + 1;
                child_cell.r = branch_cell.r / 2;
                child_cell.crange = Range(0,0);
                child_cell.brange = Range(begin + offset[i] - quad_num[i], quad_num[i]);
                child_cell.x = get_child_cell_x(branch_cell.x, branch_cell.r, i, 0);
                cells.push_back(child_cell);
                big_cells.push(insert_offset + i);
            }
        }
    }
}

void rtfmm::Tree::build_nonuniform_octree(Bodies3& bodies, vec3r x, real r, int max_n_per_cell, Cells3& cells)
{
    if(verbose) std::cout<<"build adaptive octree"<<std::endl;

    std::queue<int> big_cells;
    int num_body = bodies.size();
    Cell3 root;
    root.idx = 0;
    root.octant = 13;
    root.depth = 0;
    root.r = r;
    root.x = x;
    root.crange = Range(0,0);
    root.brange = Range(0, bodies.size());
    cells.push_back(root);
    big_cells.push(0);

    while(!big_cells.empty())
    {
        int branch_cell_idx = big_cells.front();
        Cell3 branch_cell = cells[branch_cell_idx];
        big_cells.pop();

        if(branch_cell.brange.number > max_n_per_cell)
        {
            int begin = branch_cell.brange.offset;
            int num = branch_cell.brange.number;
            int end = branch_cell.brange.offset + branch_cell.brange.number - 1;

            //divide
            int quad_num[8] = {0,0,0,0,0,0,0,0};
            int offset[8] = {0,0,0,0,0,0,0,0};
            int offset_sum = 0;
            int num_child = 0;
            for(int i = begin; i <= end; i++)
            {
                int idx = ((bodies[i].x[0] > branch_cell.x[0]) << 2) 
                        + ((bodies[i].x[1] > branch_cell.x[1]) << 1) 
                        + ((bodies[i].x[2] > branch_cell.x[2]) << 0);
                quad_num[idx]++;
            }
            for(int i = 0; i < 8; i++)
            {
                offset[i] = offset_sum;
                offset_sum += quad_num[i];
                if(quad_num[i] > 0) num_child++;
            }
            //sort
            std::vector<Body3> bodies_buffer;
            bodies_buffer.resize(num);
            for(int i = begin; i <= end; i++)
            {
                int idx = ((bodies[i].x[0] > branch_cell.x[0]) << 2) 
                        + ((bodies[i].x[1] > branch_cell.x[1]) << 1) 
                        + ((bodies[i].x[2] > branch_cell.x[2]) << 0);
                bodies_buffer[offset[idx]] = bodies[i];
                offset[idx]++;
            }
            //store
            for(int i = 0; i < num; i++)
            {
                bodies[begin + i] = bodies_buffer[i];
            }
            //create leaves
            int cnt = 0;
            int insert_offset = cells.size();
            cells[branch_cell_idx].crange = Range(insert_offset,num_child);
            for(int i = 0; i < 8; i++)
            {
                if(quad_num[i] > 0)
                {
                    Cell3 child_cell;
                    child_cell.idx = cells.size();
                    child_cell.octant = i;
                    child_cell.depth = branch_cell.depth + 1;
                    child_cell.r = branch_cell.r / 2;
                    child_cell.crange = Range(0,0);
                    child_cell.brange = Range(begin + offset[i] - quad_num[i], quad_num[i]);
                    child_cell.x = get_child_cell_x(branch_cell.x, branch_cell.r, i, 0);
                    cells.push_back(child_cell);
                    big_cells.push(insert_offset + cnt);
                    cnt++;
                }
            }
        }
    }
}

void rtfmm::Tree::build_reg_nonuniform_octree(Bodies3& bodies, vec3r x, real r, int max_n_per_cell, real rega, Bodies3& c0_rbs, Cells3& cells)
{
    if(verbose) std::cout<<"build reg adaptive octree"<<std::endl;

    std::queue<int> big_cells;
    int num_body = bodies.size();
    Cell3 root;
    root.idx = 0;
    root.octant = 13;
    root.depth = 0;
    root.r = r;
    root.x = x;
    root.crange = Range(0,0);
    root.brange = Range(0, bodies.size());
    root.rbs = c0_rbs;
    cells.push_back(root);
    big_cells.push(0);

    while(!big_cells.empty())
    {
        int branch_cell_idx = big_cells.front();
        Cell3 branch_cell = cells[branch_cell_idx];
        big_cells.pop();

        if(branch_cell.brange.number > max_n_per_cell)
        {
            int begin = branch_cell.brange.offset;
            int num = branch_cell.brange.number;
            int end = branch_cell.brange.offset + branch_cell.brange.number - 1;

            //divide
            int quad_num[8] = {0,0,0,0,0,0,0,0};
            int offset[8] = {0,0,0,0,0,0,0,0};
            int offset_sum = 0;
            int num_child = 0;
            for(int i = begin; i <= end; i++)
            {
                int idx = ((bodies[i].x[0] > branch_cell.x[0]) << 2) 
                        + ((bodies[i].x[1] > branch_cell.x[1]) << 1) 
                        + ((bodies[i].x[2] > branch_cell.x[2]) << 0);
                quad_num[idx]++;
            }
            for(int i = 0; i < 8; i++)
            {
                offset[i] = offset_sum;
                offset_sum += quad_num[i];
                if(quad_num[i] > 0) num_child++;
            }
            //sort
            std::vector<Body3> bodies_buffer;
            bodies_buffer.resize(num);
            for(int i = begin; i <= end; i++)
            {
                int idx = ((bodies[i].x[0] > branch_cell.x[0]) << 2) 
                        + ((bodies[i].x[1] > branch_cell.x[1]) << 1) 
                        + ((bodies[i].x[2] > branch_cell.x[2]) << 0);
                bodies_buffer[offset[idx]] = bodies[i];
                offset[idx]++;
            }
            //store
            for(int i = 0; i < num; i++)
            {
                bodies[begin + i] = bodies_buffer[i];
            }

            // check reg
            Bodies3 regs[8];
            for(int octant = 0; octant < 8; octant++)
            {
                vec3r cx = get_child_cell_x(branch_cell.x, branch_cell.r, octant, 0);
                for(int i = begin; i <= end; i++)
                {
                    Body3& b = bodies[i];
                    int idx = ((b.x[0] > branch_cell.x[0]) << 2) 
                            + ((b.x[1] > branch_cell.x[1]) << 1) 
                            + ((b.x[2] > branch_cell.x[2]) << 0);
                    if(idx != octant)
                    {
                        real w = LaplaceKernel::get_rega_w(b.x - cx, branch_cell.r / 2, rega);
                        if(w > 0)
                        {
                            regs[octant].push_back(b);
                        }
                    }
                }
                for(int i = 0; i < branch_cell.rbs.size(); i++)
                {
                    Body3& b = branch_cell.rbs[i];
                    int idx = ((b.x[0] > branch_cell.x[0]) << 2) 
                            + ((b.x[1] > branch_cell.x[1]) << 1) 
                            + ((b.x[2] > branch_cell.x[2]) << 0);
                    real w = LaplaceKernel::get_rega_w(b.x - cx, branch_cell.r / 2, rega);
                    if(w > 0)
                    {
                        regs[octant].push_back(b);
                    }
                }
            }

            //create leaves
            int cnt = 0;
            int insert_offset = cells.size();
            cells[branch_cell_idx].crange = Range(insert_offset,num_child);
            for(int octant = 0; octant < 8; octant++)
            {
                if(quad_num[octant] > 0 || regs[octant].size() > 0)
                {
                    Cell3 child_cell;
                    child_cell.idx = cells.size();
                    child_cell.octant = octant;
                    child_cell.depth = branch_cell.depth + 1;
                    child_cell.r = branch_cell.r / 2;
                    child_cell.crange = Range(0,0);
                    child_cell.brange = Range(begin + offset[octant] - quad_num[octant], quad_num[octant]);
                    child_cell.x = get_child_cell_x(branch_cell.x, branch_cell.r, octant, 0);
                    child_cell.rbs = regs[octant];
                    cells.push_back(child_cell);
                    big_cells.push(insert_offset + cnt);
                    cnt++;
                }
            }
        }
    }
}

void rtfmm::Tree::build_reg_uniform_octree(vec3r x, real r, int max_depth, vec3i direction, real rega, Bodies3& c0_rbs, Cells3& cells)
{
    if(verbose) std::cout<<"build reg uniform octree"<<std::endl;
    std::queue<int> big_cells;
    Cell3 root;
    root.idx = 0;
    root.octant = -1;
    root.depth = 0;
    root.r = r;
    root.x = x;
    root.crange = Range(0,0);
    root.brange = Range(0,0);
    root.rbs = c0_rbs;
    cells.push_back(root);
    big_cells.push(0);

    while(!big_cells.empty())
    {
        int branch_cell_idx = big_cells.front();
        Cell3 branch_cell = cells[branch_cell_idx];
        big_cells.pop();
        if(branch_cell.depth < max_depth)
        {
            int insert_offset = cells.size();
            cells[branch_cell_idx].crange = Range(insert_offset,8);
            // check reg
            Bodies3 regs[8];
            for(int octant = 0; octant < 8; octant++)
            {
                if((direction[0] == -1 || ((octant >> 2) & 1) == direction[0])
                 &&(direction[1] == -1 || ((octant >> 1) & 1) == direction[1])
                 &&(direction[2] == -1 || ((octant >> 0) & 1) == direction[2]))
                {
                    vec3r cx = get_child_cell_x(branch_cell.x, branch_cell.r, octant, 0);
                    for(int i = 0; i < branch_cell.rbs.size(); i++)
                    {
                        Body3& b = branch_cell.rbs[i];
                        int idx = ((b.x[0] > branch_cell.x[0]) << 2) 
                                + ((b.x[1] > branch_cell.x[1]) << 1) 
                                + ((b.x[2] > branch_cell.x[2]) << 0);
                        real w = LaplaceKernel::get_rega_w(b.x - cx, branch_cell.r / 2, rega);
                        if(w > 0)
                        {
                            regs[octant].push_back(b);
                        }
                    }
                }
            }
            for(int i = 0; i < 8; i++)
            {
                Cell3 child_cell;
                child_cell.idx = cells.size();
                child_cell.octant = i;
                child_cell.depth = branch_cell.depth + 1;
                child_cell.r = branch_cell.r / 2;
                child_cell.crange = Range(0,0);
                child_cell.brange = Range(0,0);
                child_cell.x = get_child_cell_x(branch_cell.x, branch_cell.r, i, 0);
                child_cell.rbs = regs[i];
                cells.push_back(child_cell);
                if((direction[0] == -1 || ((i >> 2) & 1) == direction[0])
                 &&(direction[1] == -1 || ((i >> 1) & 1) == direction[1])
                 &&(direction[2] == -1 || ((i >> 0) & 1) == direction[2]))
                {
                    big_cells.push(insert_offset + i);
                }
            }
        }
    }
}

rtfmm::Cells3 rtfmm::Tree::get_cells()
{
    return all_cells;
}

rtfmm::vec3r rtfmm::Tree::get_child_cell_x(vec3r x_par, real r_par, int octant, int is_periodic)
{
    vec3r x;
    if(!is_periodic)
    {
        assert_exit(octant >= 0 && octant <= 7, "octant out of range");
        if(octant == 0)      x = x_par + vec3r(-r_par/2, -r_par/2, -r_par/2);
        else if(octant == 1) x = x_par + vec3r(-r_par/2, -r_par/2,  r_par/2);
        else if(octant == 2) x = x_par + vec3r(-r_par/2,  r_par/2, -r_par/2);
        else if(octant == 3) x = x_par + vec3r(-r_par/2,  r_par/2,  r_par/2);
        else if(octant == 4) x = x_par + vec3r( r_par/2, -r_par/2, -r_par/2);
        else if(octant == 5) x = x_par + vec3r( r_par/2, -r_par/2,  r_par/2);
        else if(octant == 6) x = x_par + vec3r( r_par/2,  r_par/2, -r_par/2);
        else if(octant == 7) x = x_par + vec3r( r_par/2,  r_par/2,  r_par/2);
    }
    else
    {
        assert_exit(octant >= 0 && octant <= 26, "periodic octant out of range");
        int k = octant / 9 - 1;
        int j = (octant % 9) / 3 - 1;
        int i = (octant % 3) - 1;
        x = x_par + vec3r(i,j,k) * (r_par * 2 / 3.0);
    }
    return x;
}