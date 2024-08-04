#include "tree.h"
#include "argument.h"
#include <queue>

rtfmm::Tree::Tree()
{

}

void rtfmm::Tree::build(Bodies3& bodies, vec3r x, real r, int m, TreeType type)
{
    if(type == TreeType::uniform)
    {
        build_uniform_octree(bodies, x, r, m);
    }
    else if(type == TreeType::nonuniform)
    {
        build_nonuniform_octree(bodies, x, r, m);
    }
}

void rtfmm::Tree::build_uniform_octree(Bodies3& bodies, vec3r x, real r, int max_depth)
{
    if(verbose) std::cout<<"build complete balanced octree"<<std::endl;
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
    this->cells.push_back(root);
    big_cells.push(0);

    while(!big_cells.empty())
    {
        int branch_cell_idx = big_cells.front();
        Cell3 branch_cell = this->cells[branch_cell_idx];
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
            int insert_offset = this->cells.size();
            this->cells[branch_cell_idx].crange = Range(insert_offset,8);
            for(int i = 0; i < 8; i++)
            {
                Cell3 child_cell;
                child_cell.idx = this->cells.size();
                child_cell.octant = i;
                child_cell.depth = branch_cell.depth + 1;
                child_cell.r = branch_cell.r / 2;
                child_cell.crange = Range(0,0);
                child_cell.brange = Range(begin + offset[i] - quad_num[i], quad_num[i]);
                child_cell.x = get_child_cell_x(branch_cell.x, branch_cell.r, i, 0);
                this->cells.push_back(child_cell);
                big_cells.push(insert_offset + i);
            }
        }
    }
}

void rtfmm::Tree::build_nonuniform_octree(Bodies3& bodies, vec3r x, real r, int max_n_per_cell)
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
    this->cells.push_back(root);
    big_cells.push(0);

    while(!big_cells.empty())
    {
        int branch_cell_idx = big_cells.front();
        Cell3 branch_cell = this->cells[branch_cell_idx];
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
            int insert_offset = this->cells.size();
            if(num_child >= 1) num_child = 8;
            this->cells[branch_cell_idx].crange = Range(insert_offset,num_child);
            for(int i = 0; i < 8; i++)
            {
                if(num_child == 8)
                {
                    Cell3 child_cell;
                    child_cell.idx = this->cells.size();
                    child_cell.octant = i;
                    child_cell.depth = branch_cell.depth + 1;
                    child_cell.r = branch_cell.r / 2;
                    child_cell.crange = Range(0,0);
                    child_cell.brange = Range(begin + offset[i] - quad_num[i], quad_num[i]);
                    child_cell.x = get_child_cell_x(branch_cell.x, branch_cell.r, i, 0);
                    this->cells.push_back(child_cell);
                    big_cells.push(insert_offset + cnt);
                    cnt++;
                }
            }
        }
    }
}

rtfmm::Cells3 rtfmm::Tree::get_cells()
{
    return cells;
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