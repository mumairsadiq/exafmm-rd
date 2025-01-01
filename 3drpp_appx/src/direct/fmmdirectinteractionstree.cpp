
#include "fmmdirectinteractionstree.h"
#include <fstream>

void gmx::fmm::FMMDirectInteractionsTree::rebuild_and_reprocess_tree()
{
    fmm_cells_.clear();
    FMMTree::build_tree_uniform();
    process_tree_();
}

gmx::fmm::FMMDirectInteractionsTree::FMMDirectInteractionsTree(const FBodies &bodies, const RVec box_center, const real box_radius, const size_t max_depth) : FMMTree(bodies, box_center, box_radius, max_depth)
{
    process_tree_();
}

int gmx::fmm::FMMDirectInteractionsTree::get_neighbour_idx(RVec neighbor_center_coord)
{
    int neighbor_idx = -1;
    long neighbor_center = coord_to_long(neighbor_center_coord, box_center_, box_radius_);
    if (cells_map.find(neighbor_center) != cells_map.end())
    {
        neighbor_idx = cells_map[neighbor_center];
    }

    return neighbor_idx;
}

void gmx::fmm::FMMDirectInteractionsTree::process_tree_()
{
    FMMCells leaves = get_leaves();
    for (size_t i = 0; i < leaves.size(); i++)
    {
        leaves[i].index = i;
    }
    fmm_cells_.clear();
    fmm_cells_ = std::move(leaves);

    for (size_t i = 0; i < fmm_cells_.size(); i++)
    {
        const auto &current_cell = fmm_cells_[i];
        cells_map[coord_to_long(current_cell.center, box_center_, box_radius_)] = current_cell.index;
    }
}