
#include "fmmtree.h"

void gmx::fmm::FMMTree::build_tree_uniform()
{
    std::queue<int> big_cells;
    int num_body = bodies_.size();
    FMMCell root;
    root.radius = box_radius_;
    root.center = box_center_;
    root.index = 0;
    root.depth = 0;
    root.crange = Range(0, 0);
    for (int i = 0; i < num_body; i++)
    {
        root.bodiesIndices.push_back(i);
    }
    fmm_cells_.push_back(root);
    big_cells.push(0);

    while (!big_cells.empty())
    {
        int branch_cell_idx = big_cells.front();
        big_cells.pop();

        if (fmm_cells_[branch_cell_idx].depth < max_depth_)
        {
            // divide
            int quad_num[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            FPIndices quad_indices[8];
            int num_child = 0;
            FPIndices &bodiesIndicesBranchCell = fmm_cells_[branch_cell_idx].bodiesIndices;
            for (size_t i = 0; i < bodiesIndicesBranchCell.size(); i++)
            {
                int body_idx = bodiesIndicesBranchCell[i];
                int idx = ((bodies_[body_idx].x[0] > fmm_cells_[branch_cell_idx].center[0]) << 2) +
                          ((bodies_[body_idx].x[1] > fmm_cells_[branch_cell_idx].center[1]) << 1) +
                          ((bodies_[body_idx].x[2] > fmm_cells_[branch_cell_idx].center[2]) << 0);
                quad_indices[idx].push_back(body_idx);
                quad_num[idx]++;
            }

            // create leaves
            int cnt = 0;
            num_child = 8;
            size_t insert_offset = fmm_cells_.size();
            fmm_cells_[branch_cell_idx].crange = Range(insert_offset, num_child);

            const int branch_level = fmm_cells_[branch_cell_idx].depth;
            const auto branch_center = fmm_cells_[branch_cell_idx].center;
            const auto branch_radius = fmm_cells_[branch_cell_idx].radius;
            for (int i = 0; i < 8; i++)
            {
                RVec current_cell_center = get_child_cell_x_(branch_center, branch_radius, i);
                FMMCell child_cell;
                child_cell.center = current_cell_center;
                child_cell.radius = branch_radius / 2;
                child_cell.bodiesIndices = quad_indices[i];
                child_cell.index = fmm_cells_.size();
                child_cell.octant = i;
                child_cell.crange = Range(0, 0);
                child_cell.depth = branch_level + 1;
                child_cell.radiusParent = branch_radius;
                child_cell.centerParent = branch_center;
                fmm_cells_.push_back(child_cell);
                big_cells.push(insert_offset + cnt);
                cnt++;
            }
        }
    }
}

gmx::RVec gmx::fmm::FMMTree::get_child_cell_x_(RVec x_par, real r_par, int octant)
{
    RVec x;
    if (!(octant >= 0 && octant <= 7))
    {
        // "octant out of range"
        exit(-1);
    }
    if (octant == 0)
        x = x_par + RVec(-r_par / 2, -r_par / 2, -r_par / 2);
    else if (octant == 1)
        x = x_par + RVec(-r_par / 2, -r_par / 2, r_par / 2);
    else if (octant == 2)
        x = x_par + RVec(-r_par / 2, r_par / 2, -r_par / 2);
    else if (octant == 3)
        x = x_par + RVec(-r_par / 2, r_par / 2, r_par / 2);
    else if (octant == 4)
        x = x_par + RVec(r_par / 2, -r_par / 2, -r_par / 2);
    else if (octant == 5)
        x = x_par + RVec(r_par / 2, -r_par / 2, r_par / 2);
    else if (octant == 6)
        x = x_par + RVec(r_par / 2, r_par / 2, -r_par / 2);
    else if (octant == 7)
        x = x_par + RVec(r_par / 2, r_par / 2, r_par / 2);
    return x;
}

gmx::fmm::FMMTree::FMMTree(const FBodies &bodies, const RVec box_center, const real box_radius, const size_t max_depth)
    : bodies_(bodies), box_center_(box_center), box_radius_(box_radius), max_depth_(max_depth)
{
    build_tree_uniform();
}

gmx::fmm::FMMCells &gmx::fmm::FMMTree::get_cells() { return this->fmm_cells_; }

gmx::fmm::FMMCells gmx::fmm::FMMTree::get_leaves()
{
    FMMCells leaf_cells;
    std::stack<int> stack;

    // Start DFS from the root cell (index 0)
    stack.push(0);

    while (!stack.empty())
    {
        int current_idx = stack.top();
        stack.pop();

        FMMCell current_cell = fmm_cells_[current_idx];

        // If the current cell has children, push its children onto the stack
        if (!current_cell.isLeaf())
        {
            int child_start = current_cell.crange.offset;

            int child_end = child_start + current_cell.crange.number;
            // Push children onto the stack in reverse order for DFS traversal
            for (int i = child_end - 1; i >= child_start; --i)
            {
                stack.push(i);
            }
        }
        else
        {
            // If the cell is a leaf (no children), store it
            leaf_cells.push_back(current_cell);
        }
    }

    return leaf_cells;
}

gmx::RVec gmx::fmm::FMMTree::get_box_center() const { return box_center_; }

gmx::real gmx::fmm::FMMTree::get_box_radius() const { return box_radius_; }