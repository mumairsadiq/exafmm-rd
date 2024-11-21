
#include "fmmdirectinteractionstree.h"
#include <fstream>

gmx::fmm::FPIndices &
gmx::fmm::FMMDirectInteractionsTree::get_adjacent_cells(size_t i)
{
    if (i >= fmm_cells_.size())
    {
        exit(-1);
    }
    return this->adjacent_cells_info_[i];
}

void gmx::fmm::FMMDirectInteractionsTree::rebuild_and_reprocess_tree()
{
    fmm_cells_.clear();
    FMMTree::build_tree_non_uniform();
    find_all_adjacent_cells_();
}

inline int
gmx::fmm::FMMDirectInteractionsTree::check_adjacent_parent_(int ca_idx,
                                                            int cb_idx)
{
    const FMMCell &ca = fmm_cells_[ca_idx];
    const FMMCell &cb = fmm_cells_[cb_idx];
    if ((cb.depth - 1) < ca.depth)
    {
        return 0;
    }

    RVec dx = ca.center - cb.centerParent;
    dx[0] = fabs(dx[0]);
    dx[1] = fabs(dx[1]);
    dx[2] = fabs(dx[2]);
    real neighbour_factor = (num_neighbours - 1) * 2;
    real dist = (ca.radius + cb.radiusParent + (ca.radius * neighbour_factor)) *
                1.001; // warning : DO NOT ignore the float error
    return dx[0] <= dist && dx[1] <= dist && dx[2] <= dist;
}

void gmx::fmm::FMMDirectInteractionsTree::find_all_adjacent_cells_()
{

    std::vector<std::unordered_map<long, int>> cells_at_depth_i(max_depth_ -
                                                                min_depth_ + 1);
    for (size_t i = 0; i < fmm_cells_.size(); i++)
    {
        const auto &current_cell = fmm_cells_[i];
        if (current_cell.depth >= min_depth_)
        {
            cells_at_depth_i[current_cell.depth - min_depth_][coord_to_long(
                current_cell.center, box_center_, box_radius_)] =
                current_cell.index;
        }
    }

    std::vector<FPIndices> fmm_cells_children_leaves(fmm_cells_.size());

    for (int i = max_depth_; i >= min_depth_; i--)
    {
        for (auto &v : cells_at_depth_i[i - min_depth_])
        {
            FMMCell &cell_curr = fmm_cells_[v.second];

            if (!cell_curr.isLeaf())
            {
                size_t total_leaves_to_insert = 0;
                for (int k = 0; k < cell_curr.crange.number; k++)
                {
                    total_leaves_to_insert +=
                        fmm_cells_children_leaves[cell_curr.crange.offset + k]
                            .size();
                }
                fmm_cells_children_leaves[cell_curr.index].reserve(
                    total_leaves_to_insert);
                for (int k = 0; k < cell_curr.crange.number; k++)
                {
                    std::copy(
                        fmm_cells_children_leaves[cell_curr.crange.offset + k]
                            .begin(),
                        fmm_cells_children_leaves[cell_curr.crange.offset + k]
                            .end(),
                        std::back_inserter(
                            fmm_cells_children_leaves[cell_curr.index]));
                }
            }
            else
            {
                fmm_cells_children_leaves[cell_curr.index].push_back(
                    cell_curr.index);
            }
        }
    }

    FMMCells leaves = getLeaves();
    std::unordered_map<int, int> full_tree_to_leaves_map;
    for (size_t i = 0; i < leaves.size(); i++)
    {
        full_tree_to_leaves_map[leaves[i].index] = i;
    }

    FPIndices full_tree_to_leaves_indices(fmm_cells_.size(), 0);
    for (size_t i = 0; i < fmm_cells_.size(); i++)
    {
        const FMMCell &curr_cell = fmm_cells_[i];
        if (curr_cell.isLeaf())
        {
            full_tree_to_leaves_indices[i] =
                full_tree_to_leaves_map[curr_cell.index];
        }
    }
    full_tree_to_leaves_map.clear();

    adjacent_cells_info_.resize(leaves.size());
    for (size_t i = 0; i < leaves.size(); i++)
    {
        adjacent_cells_info_[i].push_back(i);
    }

    for (size_t i = 0; i < cells_at_depth_i.size(); i++)
    {
        for (const auto &val : cells_at_depth_i[i])
        {
            const int tar_cell_idx = val.second;
            const auto &target_cell = fmm_cells_[tar_cell_idx];
            const int tar_cell_idx_l =
                full_tree_to_leaves_indices[tar_cell_idx];

            if (target_cell.isLeaf())
            {
                real ddistance = target_cell.radius * 2;
                for (int dz = -1 * num_neighbours; dz <= num_neighbours; dz++)
                {
                    for (int dy = -1 * num_neighbours; dy <= num_neighbours;
                         dy++)
                    {
                        for (int dx = -1 * num_neighbours; dx <= num_neighbours;
                             dx++)
                        {
                            // Skip the current cell (dx == 0, dy == 0, dz == 0)
                            if (dx == 0 && dy == 0 && dz == 0)
                            {
                                continue;
                            }

                            // Compute the neighboring cell's center coordinates
                            const long neighbor_center = coord_to_long(
                                target_cell.center + RVec(dx * ddistance,
                                                          dy * ddistance,
                                                          dz * ddistance),
                                box_center_, box_radius_);

                            if (cells_at_depth_i[i].find(neighbor_center) !=
                                cells_at_depth_i[i].end())
                            {
                                int src_idx =
                                    cells_at_depth_i[i][neighbor_center];

                                FMMCell &source_cell = fmm_cells_[src_idx];

                                if (source_cell.isLeaf())
                                {

                                    const int src_cell_idx_l =
                                        full_tree_to_leaves_indices[src_idx];

                                    adjacent_cells_info_[tar_cell_idx_l]
                                        .push_back(src_cell_idx_l);

                                    if (target_cell.depth < source_cell.depth)
                                    {
                                        adjacent_cells_info_[src_cell_idx_l]
                                            .push_back(tar_cell_idx_l);
                                    }
                                }
                                else
                                {
                                    const auto &leave_cells_to_curr_source =
                                        fmm_cells_children_leaves[source_cell
                                                                      .index];

                                    for (size_t k = 0;
                                         k < leave_cells_to_curr_source.size();
                                         k++)
                                    {
                                        int csrc_idx =
                                            leave_cells_to_curr_source[k];

                                        const int csrc_cell_idx_l =
                                            full_tree_to_leaves_indices
                                                [csrc_idx];
                                        if (check_adjacent_parent_(tar_cell_idx,
                                                                   csrc_idx))
                                        {
                                            adjacent_cells_info_[tar_cell_idx_l]
                                                .push_back(csrc_cell_idx_l);

                                            adjacent_cells_info_
                                                [csrc_cell_idx_l]
                                                    .push_back(tar_cell_idx_l);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fmm_cells_children_leaves.clear();

    for (size_t i = 0; i < leaves.size(); i++)
    {
        leaves[i].index = i;
    }
    full_tree_to_leaves_indices.clear();
    fmm_cells_.clear();
    fmm_cells_ = std::move(leaves);

    // std::ofstream fout_adj("my_adjacency.txt");
    // // Dump interaction information
    // for (size_t i = 0; i < fmm_cells_.size(); i++)
    // {
    //     FMMCell &ctar = fmm_cells_[i];

    //     fout_adj << "Cell " << ctar.center << "\t"
    //              << " at " << ctar.depth << " have following adjacents"
    //              << ": \n";
    //     for (auto &adj : get_adjacent_cells(i))
    //     {
    //         FMMCell &csrc = fmm_cells_[adj];

    //         fout_adj << csrc.center << "\n";
    //     }
    //     fout_adj << "\n";
    // }

    // fout_adj.
}

gmx::fmm::FMMDirectInteractionsTree::FMMDirectInteractionsTree(
    const FBodies &bodies, const RVec box_center, const real box_radius,
    const size_t cell_limit_param, const bool is_tree_uniform)
    : FMMTree(bodies, box_center, box_radius, cell_limit_param,
              is_tree_uniform),
      num_neighbours(1)
{
    this->process_tree_();
}

void gmx::fmm::FMMDirectInteractionsTree::process_tree_()
{
    find_all_adjacent_cells_();
}