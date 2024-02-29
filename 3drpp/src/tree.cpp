#include "tree.h"
#include "argument.h"
#include <queue>
#include "kernel.h"
#include "mathfunc.h"

rtfmm::Tree::Tree()
{

}

void rtfmm::Tree::build(Bodies3& bodies, vec3r x, real r, int m, real rega, TreeType type)
{
    if(type == TreeType::uniform)
    {
        build_uniform_octree(bodies, x, r, m);
    }
    else if(type == TreeType::nonuniform)
    {
        build_nonuniform_octree(bodies, x, r, m);
    }
    else if(type == TreeType::reg_nonuniform)
    {
        build_nonuniform_octree2(bodies, x, r, m, rega);
        //build_reg_nonuniform_octree(bodies, x, r, m, rega);
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
    root.split_policy = SplitPolicy::SPLIT8;
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
            std::vector<int> quad_num(static_cast<int>(branch_cell.split_policy));
            std::vector<int> offset(static_cast<int>(branch_cell.split_policy));
            int offset_sum = 0;
            for(int i = begin; i <= end; i++)
            {
                int idx = get_body_octant(bodies[i].x, branch_cell.x, branch_cell.r, branch_cell.split_policy);
                quad_num[idx]++;
            }
            for(int i = 0; i < static_cast<int>(branch_cell.split_policy); i++)
            {
                offset[i] = offset_sum;
                offset_sum += quad_num[i];
            }
            //sort
            std::vector<Body3> bodies_buffer;
            bodies_buffer.resize(num);
            for(int i = begin; i <= end; i++)
            {
                int idx = get_body_octant(bodies[i].x, branch_cell.x, branch_cell.r, SplitPolicy::SPLIT8);
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
            this->cells[branch_cell_idx].crange = Range(insert_offset,static_cast<int>(branch_cell.split_policy));
            for(int i = 0; i < static_cast<int>(branch_cell.split_policy); i++)
            {
                Cell3 child_cell;
                child_cell.idx = this->cells.size();
                child_cell.octant = i;
                child_cell.depth = branch_cell.depth + 1;
                child_cell.r = branch_cell.r / 2;
                child_cell.crange = Range(0,0);
                child_cell.brange = Range(begin + offset[i] - quad_num[i], quad_num[i]);
                child_cell.x = get_child_cell_x(branch_cell.x, branch_cell.r, i, 0);
                child_cell.split_policy = SplitPolicy::SPLIT8;
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
    root.split_policy = SplitPolicy::SPLIT8;
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
            //count
            std::vector<int> quad_num(static_cast<int>(branch_cell.split_policy));
            for(int i = begin; i <= end; i++)
            {
                int idx = get_body_octant(bodies[i].x, branch_cell.x, branch_cell.r, branch_cell.split_policy);
                quad_num[idx]++;
            }
            std::vector<int> offset = exclusive_scan(quad_num);
            //sort
            std::vector<Body3> bodies_buffer;
            bodies_buffer.resize(num);
            for(int i = begin; i <= end; i++)
            {
                int idx = get_body_octant(bodies[i].x, branch_cell.x, branch_cell.r, branch_cell.split_policy);
                bodies_buffer[offset[idx]] = bodies[i];
                offset[idx]++;
            }
            for(int i = 0; i < num; i++)
            {
                bodies[begin + i] = bodies_buffer[i];
            }
            //create leaves
            int num_child = std::count_if(std::begin(quad_num), std::end(quad_num), [](int v){return v > 0;});
            int cnt = 0;
            int insert_offset = this->cells.size();
            this->cells[branch_cell_idx].crange = Range(insert_offset,num_child);
            for(int i = 0; i < static_cast<int>(branch_cell.split_policy); i++)
            {
                if(quad_num[i] > 0)
                {
                    Cell3 child_cell;
                    child_cell.idx = this->cells.size();
                    child_cell.octant = i;
                    child_cell.depth = branch_cell.depth + 1;
                    child_cell.r = branch_cell.r / (branch_cell.split_policy == SplitPolicy::SPLIT8 ? 2 : 3);
                    child_cell.crange = Range(0,0);
                    child_cell.brange = Range(begin + offset[i] - quad_num[i], quad_num[i]);
                    child_cell.x = get_child_cell_x(branch_cell.x, branch_cell.r, i, branch_cell.split_policy == SplitPolicy::SPLIT27);
                    child_cell.split_policy = SplitPolicy::SPLIT8;
                    this->cells.push_back(child_cell);
                    big_cells.push(insert_offset + cnt);
                    cnt++;
                }
            }
        }
    }
}

void rtfmm::Tree::build_nonuniform_octree2(Bodies3& bodies, vec3r x, real r, int max_n_per_cell, real rega)
{
    if(verbose) std::cout<<"build adaptive octree 2"<<std::endl;

    std::queue<int> big_cells;
    int num_body = bodies.size();
    Cell3 root;
    root.idx = 0;
    root.octant = 13;
    root.depth = -1;
    root.r = r * 3;
    root.x = x;
    root.crange = Range(0,0);
    root.brange = Range(0, bodies.size());
    root.split_policy = SplitPolicy::SPLIT27;
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
            //count
            std::vector<int> quad_num(static_cast<int>(branch_cell.split_policy));
            for(int i = begin; i <= end; i++)
            {
                int idx = get_body_octant(bodies[i].x, branch_cell.x, branch_cell.r, branch_cell.split_policy);
                quad_num[idx]++;
            }
            std::vector<int> offset = exclusive_scan(quad_num);
            //sort
            std::vector<Body3> bodies_buffer;
            bodies_buffer.resize(num);
            for(int i = begin; i <= end; i++)
            {
                int idx = get_body_octant(bodies[i].x, branch_cell.x, branch_cell.r, branch_cell.split_policy);
                bodies_buffer[offset[idx]] = bodies[i];
                offset[idx]++;
            }
            for(int i = 0; i < num; i++)
            {
                bodies[begin + i] = bodies_buffer[i];
            }
            //add reg
            std::vector<int> quad_reg_num(static_cast<int>(branch_cell.split_policy));
            real r_child = branch_cell.r / (branch_cell.split_policy == SplitPolicy::SPLIT8 ? 2 : 3);
            for(int octant = 0; octant < int(branch_cell.split_policy); octant++)
            {
                vec3r x_child = get_child_cell_x(branch_cell.x, branch_cell.r, octant, branch_cell.split_policy == SplitPolicy::SPLIT27);
                for(int i = begin; i <= end; i++)
                {
                    Body3& b = bodies[i];
                    int idx = get_body_octant(b.x, branch_cell.x, branch_cell.r, branch_cell.split_policy);
                    if(idx != octant)  // if b is included in this octant before, it will not be a reg body of this octant
                    {
                        vec3r dx = (b.x - x_child).abs() - r_child;
                        real w = LaplaceKernel::get_rega_w(b.x - x_child, r_child, rega);
                        if(w > 0 && dx[0] < rega && dx[1] < rega && dx[2] < rega)
                        {
                            quad_reg_num[octant]++;
                        }
                    }
                }
                for(int i = 0; i < branch_cell.reg_bs.size(); i++)
                {
                    Body3& b = branch_cell.reg_bs[i];
                    int idx = get_body_octant(b.x, branch_cell.x, branch_cell.r, branch_cell.split_policy);
                    vec3r dx = (b.x - x_child).abs() - r_child;
                    real w = LaplaceKernel::get_rega_w(b.x - x_child, r_child, rega);
                    if(w > 0 && dx[0] < rega && dx[1] < rega && dx[2] < rega)
                    {
                        quad_reg_num[octant]++;
                    }
                }
            }

            //create leaves
            //int num_child = std::count_if(std::begin(quad_num), std::end(quad_num), [](int v){return v > 0;});
            int num_child = 0;
            for(int i = 0; i < int(branch_cell.split_policy); i++)
            {
                if(quad_num[i] || quad_reg_num[i])
                    num_child++;
            }
            int cnt = 0;
            int insert_offset = this->cells.size();
            this->cells[branch_cell_idx].crange = Range(insert_offset,num_child);
            for(int i = 0; i < static_cast<int>(branch_cell.split_policy); i++)
            {
                if(quad_num[i] > 0
                 || quad_reg_num[i] > 0
                 )
                {
                    Cell3 child_cell;
                    child_cell.idx = this->cells.size();
                    child_cell.octant = i;
                    child_cell.depth = branch_cell.depth + 1;
                    child_cell.r = r_child;
                    child_cell.crange = Range(0,0);
                    child_cell.brange = Range(begin + offset[i] - quad_num[i], quad_num[i]);
                    child_cell.x = get_child_cell_x(branch_cell.x, branch_cell.r, i, branch_cell.split_policy == SplitPolicy::SPLIT27);
                    child_cell.split_policy = SplitPolicy::SPLIT8;
                    for(int j = begin; j <= end; j++)
                    {
                        Body3& b = bodies[j];
                        int idx = get_body_octant(b.x, branch_cell.x, branch_cell.r, branch_cell.split_policy);
                        if(idx != i)
                        {
                            vec3r dx = (b.x - child_cell.x).abs() - r_child;
                            real w = LaplaceKernel::get_rega_w(b.x - child_cell.x, r_child, rega);
                            if(w > 0)
                            {
                                child_cell.reg_bs.push_back(b);
                                if(b.idx == 17846 && w > 0 && w != 1)
                                {
                                    print_body(b);
                                    std::cout<<child_cell<<std::endl;
                                    printf("w = %.4f\n", w);
                                }
                            }
                        }
                    }
                    for(int j = 0; j < branch_cell.reg_bs.size(); j++)
                    {
                        Body3& b = branch_cell.reg_bs[j];
                        int idx = get_body_octant(b.x, branch_cell.x, branch_cell.r, branch_cell.split_policy);
                        vec3r dx = (b.x - child_cell.x).abs() - r_child;
                        real w = LaplaceKernel::get_rega_w(b.x - child_cell.x, r_child, rega);
                        if(w > 0)
                        {
                            child_cell.reg_bs.push_back(b);
                        }
                    }
                    this->cells.push_back(child_cell);
                    big_cells.push(insert_offset + cnt);
                    cnt++;
                    if(branch_cell.depth == -1)
                    std::cout<<"add : "<<child_cell<<std::endl;
                }
            }
            //exit(0);
        }
    }
}

void rtfmm::Tree::build_reg_nonuniform_octree(Bodies3& bodies, vec3r x, real r, int max_n_per_cell, real rega)
{
    if(verbose) std::cout<<"build regularized adaptive octree"<<std::endl;

    cells.clear();

    Cell3 root;
    root.idx = 0;
    root.octant = 13;
    root.depth = rega > 0 ? -1 : 0;
    root.r = rega > 0 ? r * 3 : r;
    root.x = x;
    root.crange = Range(0,0);
    for(int i = 0; i < bodies.size(); i++)
    {
        root.bs.push_back(bodies[i]);
    }
    cells.push_back(root);

    std::queue<int> big_cells;
    big_cells.push(0);

    while(!big_cells.empty())
    {
        int parent_idx = big_cells.front();
        big_cells.pop();
        int num1 = cells.size();
        if(cells[parent_idx].bs.size() > max_n_per_cell)
        {
            split_cell(cells[parent_idx], rega, (rega > 0 && parent_idx == 0) ? SplitPolicy::SPLIT27 : SplitPolicy::SPLIT8, bodies, cells);
            int num2 = cells.size();
            cells[parent_idx].crange = Range(num1, num2 - num1);
            for(int j = 0; j < num2 - num1; j++)
            {
                big_cells.push(num1 + j);
            }
        }
    }

    /*for(int i = 0; i < cells.size(); i++)
    {
        std::cout<<cells[i]<<std::endl;
    }*/
}

void rtfmm::Tree::split_cell(Cell3 parent, real rega, SplitPolicy split_policy, Bodies3& bs, Cells3& cells)
{
    int num_octant = int(split_policy);
    real r_child = split_policy == SplitPolicy::SPLIT8 ? parent.r / 2 : parent.r / 3;
    std::vector<int> num_body_each_octant(num_octant, 0);
    for(int octant = 0; octant < num_octant; octant++)
    {
        vec3r x_child = get_child_cell_x(parent.x, parent.r, octant, split_policy == SplitPolicy::SPLIT27);
        for(int i = 0; i < parent.bs.size(); i++)
        {
            Body3& b = parent.bs[i];
            vec3r dx = (b.x - x_child).abs() - r_child;
            real w = LaplaceKernel::get_rega_w(b.x - x_child, r_child, rega);
            if(w > 0 && dx[0] < rega && dx[1] < rega && dx[2] < rega)
            {
                if(num_body_each_octant[octant] == 0)
                {
                    Cell3 child;
                    child.idx = cells.size();
                    child.octant = octant;
                    child.depth = parent.depth + 1;
                    child.r = r_child;
                    child.x = x_child;
                    child.crange = Range(0,0);
                    cells.push_back(child);
                }
                cells[cells.size() - 1].bs.push_back(b);
                cells[cells.size() - 1].ws.push_back(w);
                if(w != 1)
                {
                    printf("w = %.4f\n", w);
                }
                num_body_each_octant[octant]++;
            }
        }
    }
}

int rtfmm::Tree::get_body_octant(const vec3r& bx, const vec3r& cx, const real& cr, SplitPolicy split_policy)
{
    int idx = -1;
    
    if(split_policy == SplitPolicy::SPLIT8)
    {
        idx = ((bx[0] > cx[0]) << 2) 
            + ((bx[1] > cx[1]) << 1) 
            + ((bx[2] > cx[2]) << 0);
    }
    else if(split_policy == SplitPolicy::SPLIT27)
    {
        real r6 = cr / 3;
        vec3i idx3(
            bx[0] > cx[0] + r6 ? 2 : (bx[0] > cx[0] - r6 ? 1 : 0),
            bx[1] > cx[1] + r6 ? 2 : (bx[1] > cx[1] - r6 ? 1 : 0),
            bx[2] > cx[2] + r6 ? 2 : (bx[2] > cx[2] - r6 ? 1 : 0)
        ); 
        idx = idx3[0] * 9 + idx3[1] * 3 + idx3[2];
    }

    return idx;
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