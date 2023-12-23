#include "tree.h"
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
    std::cout<<"build complete balanced octree"<<std::endl;
    std::queue<int> big_cells;
    int num_body = bodies.size();
    Cell3 root;
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
            this->cells[branch_cell_idx].crange = OffsetAndNumber(insert_offset,8);
            for(int i = 0; i < 8; i++)
            {
                Cell3 leaf_cell;
                leaf_cell.depth = branch_cell.depth + 1;
                leaf_cell.r = branch_cell.r / 2;
                leaf_cell.crange = OffsetAndNumber(0,0);
                leaf_cell.brange = OffsetAndNumber(begin + offset[i] - quad_num[i], quad_num[i]);
                if(i == 0)      leaf_cell.x = branch_cell.x + vec3r(-branch_cell.r/2, -branch_cell.r/2, -branch_cell.r/2);
                else if(i == 1) leaf_cell.x = branch_cell.x + vec3r(-branch_cell.r/2, -branch_cell.r/2,  branch_cell.r/2);
                else if(i == 2) leaf_cell.x = branch_cell.x + vec3r(-branch_cell.r/2,  branch_cell.r/2, -branch_cell.r/2);
                else if(i == 3) leaf_cell.x = branch_cell.x + vec3r(-branch_cell.r/2,  branch_cell.r/2,  branch_cell.r/2);
                else if(i == 4) leaf_cell.x = branch_cell.x + vec3r( branch_cell.r/2, -branch_cell.r/2, -branch_cell.r/2);
                else if(i == 5) leaf_cell.x = branch_cell.x + vec3r( branch_cell.r/2, -branch_cell.r/2,  branch_cell.r/2);
                else if(i == 6) leaf_cell.x = branch_cell.x + vec3r( branch_cell.r/2,  branch_cell.r/2, -branch_cell.r/2);
                else if(i == 7) leaf_cell.x = branch_cell.x + vec3r( branch_cell.r/2,  branch_cell.r/2,  branch_cell.r/2);
                this->cells.push_back(leaf_cell);
                big_cells.push(insert_offset + i);
            }
        }
    }
}

void rtfmm::Tree::build_nonuniform_octree(Bodies3& bodies, vec3r x, real r, int max_n_per_cell)
{
    std::cout<<"build adaptive octree"<<std::endl;
}

rtfmm::Cells3 rtfmm::Tree::get_cells()
{
    return cells;
}