#include "fmm.h"
#include <fstream>
#include <iomanip>

// assuming coordinates are simd packed (like xxx..yyy..zzz.. as per simd width)
gmx::fmm::FMMDirectInteractions::FMMDirectInteractions(
    const std::vector<RVec> coordinates, const std::vector<real> charges,
    const RVec box_center, const real box_radius, const size_t cell_limit_param,
    const real reg_alpha)
    : bodies_all_(coordinates, charges),
      fmm_weights_eval_(box_center, box_radius, reg_alpha),
      fmm_direct_interactions_tree_(bodies_all_, box_center, box_radius,
                                    cell_limit_param, true)
{
    compute_weights_();
}

bool gmx::fmm::FMMDirectInteractions::is_point_within_radius(const RVec &point1,
                                                             const RVec &point2,
                                                             double radius)
{
    RVec dx = {point1[0] - point2[0], point1[1] - point2[1],
               point1[2] - point2[2]};
    double distance_squared = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    double radius_squared = radius * radius;
    return distance_squared <= radius_squared;
}

void gmx::fmm::FMMDirectInteractions::compute_weights_()
{
    atoms_interactions_list.clear();
    atoms_interactions_list.resize(bodies_all_.size());

    atoms_interactions_weights_src.clear();
    atoms_interactions_weights_src.resize(bodies_all_.size());

    atoms_interactions_weights_tar.clear();
    atoms_interactions_weights_tar.resize(bodies_all_.size());

    const FMMCells &fmm_cells = fmm_direct_interactions_tree_.get_cells();

    // needs optimization
    std::vector<std::unordered_set<int>> adj_second_cells(fmm_cells.size());
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        const auto &adj1_cells_indices =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);

        const std::unordered_set<int> &adj1_cells_indices_set =
            fmm_direct_interactions_tree_.get_adjacent_cells_set(k);

        for (size_t j = 0; j < adj1_cells_indices.size(); j++)
        {
            const int adj1_cell_idx = adj1_cells_indices[j];
            const FMMCell &adj1_cell = fmm_cells[adj1_cell_idx];
            if (adj1_cell.index != cell.index)
            {
                const FPIndices &adj2_cells_indices =
                    fmm_direct_interactions_tree_.get_adjacent_cells(
                        adj1_cell_idx);

                for (size_t l = 0; l < adj2_cells_indices.size(); l++)
                {
                    if (adj1_cells_indices_set.find(adj2_cells_indices[l]) ==
                        adj1_cells_indices_set.end())
                    {
                        adj_second_cells[k].insert(adj2_cells_indices[l]);
                    }
                }
            }
        }
    }

    std::vector<FPIndices> boundary_bodies_idxs(fmm_cells.size());
    std::vector<std::vector<bool>> is_part_of_reg_region(fmm_cells.size());

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        is_part_of_reg_region[k].resize(cell.bodiesIndices.size());
        int bidx = 0;
        for (const int &body_idx : cell.bodiesIndices)
        {
            const FBody &body = bodies_all_[body_idx];
            const real w = fmm_weights_eval_.compute_weight_within_cell(
                body.x, cell.center, cell.radius, false);

            if (w < 1)
            {
                boundary_bodies_idxs[k].push_back(body_idx);
                is_part_of_reg_region[k][bidx] = true;
            }
            else
            {
                is_part_of_reg_region[k][bidx] = false;
            }
            bidx++;
        }
    }

    std::vector<std::unordered_map<int, real>> bodies_weights_to_adj_cells(
        bodies_all_.size());

    // default interaction weight is one
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        const FPIndices adj1_cells_indices =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);

        for (const int body_idx_tar : cell.bodiesIndices)
        {
            for (size_t j = 0; j < adj1_cells_indices.size(); j++)
            {
                bodies_weights_to_adj_cells[body_idx_tar]
                                           [adj1_cells_indices[j]] = 1;
            }
        }
    }

    std::vector<FPIndices> reg_bodies_cell(fmm_cells.size());
    std::vector<std::vector<real>> reg_weights_cell(fmm_cells.size());
    std::vector<FPIndices> reg_bodies_from_cell_idx(fmm_cells.size());

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {

        const FMMCell &cell = fmm_cells[k];

        const auto &adj1_cell =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);

        // find out bodies from adjacent cells that are influencing this cell
        for (size_t j = 0; j < adj1_cell.size(); j++)
        {
            const int adj_cell_idx = adj1_cell[j];
            if (cell.index != adj_cell_idx)
            {
                // checkout for regularization boundary bodies from adjacent
                // cells having weight less than one in their cell
                for (const int &body_idx : boundary_bodies_idxs[adj_cell_idx])
                {

                    const FBody &body = bodies_all_[body_idx];
                    const RVec dx = body.x - cell.center;
                    const real w =
                        fmm_weights_eval_.compute_weight_outside_cell(
                            body.x, cell.center, cell.radius, false);

                    if (w > 0)
                    {
                        reg_bodies_cell[k].push_back(body_idx);
                        reg_weights_cell[k].push_back(w);
                        reg_bodies_from_cell_idx[k].push_back(adj_cell_idx);
                    }
                }
            }
        }
    }

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        const auto &adj1_cells =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);

        size_t bidxt = 0;
        // should be possible to iterate over only bodies in regularization
        // region to be seen once all things work correctly
        for (const int body_idx_tar : cell.bodiesIndices)
        {
            if (is_part_of_reg_region[k][bidxt] == true)
            {
                const FBody &body_tar = bodies_all_[body_idx_tar];
                for (size_t j = 0; j < adj1_cells.size(); j++)
                {
                    const int adj_cell_idx = adj1_cells[j];
                    const FMMCell &adj_cell = fmm_cells[adj_cell_idx];

                    // needs optimization, this is brute force
                    // approach to compute w_total, could be optimized by
                    // using geometric position of the body
                    real w_total = 0;
                    const std::unordered_set<int> &adj_cells_src_set =
                        fmm_direct_interactions_tree_.get_adjacent_cells_set(
                            adj_cell_idx);
                    for (size_t l = 0; l < adj1_cells.size(); l++)
                    {
                        if (adj1_cells[l] != cell.index)
                        {
                            if (adj_cells_src_set.find(adj1_cells[l]) !=
                                adj_cells_src_set.end())
                            {
                                const real w_temp =
                                    fmm_weights_eval_
                                        .compute_weight_outside_cell(
                                            body_tar.x,
                                            fmm_cells[adj1_cells[l]].center,
                                            fmm_cells[adj1_cells[l]].radius,
                                            false);
                                w_total += w_temp;
                            }
                        }
                        else
                        {
                            const real w_temp =
                                fmm_weights_eval_.compute_weight_within_cell(
                                    body_tar.x, cell.center, cell.radius,
                                    false);
                            w_total += w_temp;
                        }
                    }
                    // optimize till this point
                    bodies_weights_to_adj_cells[body_idx_tar][adj_cell_idx] =
                        w_total;
                }
            }
            bidxt++;
        }
    }

    // interactions within same cell for each body
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        for (const int &body_idx_tar : cell.bodiesIndices)
        {

            const FBody &body_tar = bodies_all_[body_idx_tar];
            for (const int body_idx_src : cell.bodiesIndices)
            {
                if (body_idx_tar != body_idx_src)
                {
                    const FBody &body_src = bodies_all_[body_idx_src];

                    // it is possible that target weight is split into current
                    // cell and first adjacent cell, and source may only
                    // interact from adjacent cell
                    real target_weight_within_cell =
                        fmm_weights_eval_.compute_weight_within_cell(
                            body_tar.x, cell.center, cell.radius, false);

                    real source_weight_within_cell =
                        fmm_weights_eval_.compute_weight_within_cell(
                            body_src.x, cell.center, cell.radius, false);

                    atoms_interactions_list[body_idx_tar].push_back(
                        body_idx_src);

                    atoms_interactions_weights_tar[body_idx_tar].push_back(
                        target_weight_within_cell);

                    // within cell source weight is always one
                    atoms_interactions_weights_src[body_idx_tar].push_back(1);
                }
            }
        }
    }

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        const auto &adj1_cells_indices =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);

        size_t bidxt = 0;
        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];

            for (size_t j = 0; j < adj1_cells_indices.size(); j++)
            {
                const int adj1_cell_idx = adj1_cells_indices[j];
                const FMMCell &adj1_cell = fmm_cells[adj1_cell_idx];

                if (adj1_cell.index != cell.index)
                {
                    for (const int body_idx_src : adj1_cell.bodiesIndices)
                    {
                        const FBody &body_src = bodies_all_[body_idx_src];

                        real target_weight_within_cell =
                            fmm_weights_eval_.compute_weight_within_cell(
                                body_tar.x, cell.center, cell.radius, false);

                        atoms_interactions_list[body_idx_tar].push_back(
                            body_idx_src);

                        atoms_interactions_weights_tar[body_idx_tar].push_back(
                            target_weight_within_cell);

                        atoms_interactions_weights_src[body_idx_tar].push_back(
                            bodies_weights_to_adj_cells[body_idx_src]
                                                       [cell.index]);
                    }
                }
            }
            bidxt++;
        }
    }

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        const auto &adj1_cells_indices =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);

        const std::unordered_set<int> &adj1_cells_indices_set =
            fmm_direct_interactions_tree_.get_adjacent_cells_set(k);

        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];

            for (size_t j = 0; j < adj1_cells_indices.size(); j++)
            {
                const int adj1_cell_idx = adj1_cells_indices[j];
                const FMMCell &adj1_cell = fmm_cells[adj1_cell_idx];
                if (cell.index != adj1_cell_idx)
                {
                    real target_weight_in_adj1 =
                        fmm_weights_eval_.compute_weight_outside_cell(
                            body_tar.x, adj1_cell.center, adj1_cell.radius,
                            false);

                    if (target_weight_in_adj1 > 0)
                    {
                        const auto &adj2_cells_indices =
                            fmm_direct_interactions_tree_.get_adjacent_cells(
                                adj1_cell_idx);

                        const auto &adj2_cells_indices_set =
                            fmm_direct_interactions_tree_
                                .get_adjacent_cells_set(adj1_cell_idx);

                        for (size_t l = 0; l < adj2_cells_indices.size(); l++)
                        {
                            const int adj2_cell_idx = adj2_cells_indices[l];

                            const FMMCell &adj2_cell = fmm_cells[adj2_cell_idx];

                            for (const int body_idx_src :
                                 adj2_cell.bodiesIndices)
                            {
                                if (body_idx_tar != body_idx_src)
                                {
                                    const FBody &body_src =
                                        bodies_all_[body_idx_src];

                                    atoms_interactions_list[body_idx_tar]
                                        .push_back(body_idx_src);

                                    atoms_interactions_weights_tar[body_idx_tar]
                                        .push_back(target_weight_in_adj1);

                                    atoms_interactions_weights_src[body_idx_tar]
                                        .push_back(
                                            bodies_weights_to_adj_cells
                                                [body_idx_src][adj1_cell_idx]);
                                }
                            }

                            for (size_t j2 = 0; j2 < adj1_cells_indices.size();
                                 j2++)
                            {
                                if (adj2_cells_indices_set.find(
                                        adj1_cells_indices[j2]) ==
                                    adj2_cells_indices_set.end())
                                {
                                    for (const int body_idx_src :
                                         boundary_bodies_idxs
                                             [adj1_cells_indices[j2]])
                                    {
                                        const FBody &body_src =
                                            bodies_all_[body_idx_src];

                                        const real wsrc_partial =
                                            fmm_weights_eval_
                                                .compute_weight_outside_cell(
                                                    body_src.x,
                                                    adj2_cell.center,
                                                    adj2_cell.radius, false);

                                        if (wsrc_partial != 0)
                                        {
                                            atoms_interactions_list
                                                [body_idx_tar]
                                                    .push_back(body_idx_src);

                                            atoms_interactions_weights_tar
                                                [body_idx_tar]
                                                    .push_back(
                                                        target_weight_in_adj1);

                                            atoms_interactions_weights_src
                                                [body_idx_tar]
                                                    .push_back(wsrc_partial);
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

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        const auto &adj1_cells_indices =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);

        const std::unordered_set<int> &adj1_cells_indices_set =
            fmm_direct_interactions_tree_.get_adjacent_cells_set(k);

        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];
            for (const int adj2_cell_idx : adj_second_cells[k])
            {
                const FMMCell &adj2_cell = fmm_cells[adj2_cell_idx];

                for (const int body_idx_src :
                     boundary_bodies_idxs[adj2_cell_idx])
                {
                    const FBody &body_src = bodies_all_[body_idx_src];
                    for (size_t j = 0; j < adj1_cells_indices.size(); j++)
                    {
                        const int adj1_cell_idx = adj1_cells_indices[j];
                        const FMMCell &adj1_cell = fmm_cells[adj1_cell_idx];

                        if (adj1_cell_idx != cell.index)
                        {

                            real wtar_influencing_cell =
                                fmm_weights_eval_.compute_weight_within_cell(
                                    body_tar.x, cell.center, cell.radius,
                                    false);

                            std::unordered_set<int> adj1b_set =
                                fmm_direct_interactions_tree_
                                    .get_adjacent_cells_set(adj1_cell_idx);

                            std::unordered_set<int> adj2b_set =
                                fmm_direct_interactions_tree_
                                    .get_adjacent_cells_set(adj2_cell_idx);

                            for (const int adj1b_cell_idx : adj1b_set)
                            {
                                if (adj1b_cell_idx != adj1_cell_idx &&
                                    adj1b_cell_idx != cell.index)
                                {
                                    if (adj2b_set.find(adj1b_cell_idx) ==
                                        adj2b_set.end())
                                    {
                                        wtar_influencing_cell +=
                                            fmm_weights_eval_
                                                .compute_weight_outside_cell(
                                                    body_tar.x,
                                                    fmm_cells[adj1b_cell_idx]
                                                        .center,
                                                    fmm_cells[adj1b_cell_idx]
                                                        .radius,
                                                    false);
                                    }
                                }
                            }

                            real wsrc2_adj1_partial =
                                fmm_weights_eval_.compute_weight_outside_cell(
                                    body_src.x, adj1_cell.center,
                                    adj1_cell.radius, false);

                            if (wsrc2_adj1_partial > 0)
                            {
                                atoms_interactions_list[body_idx_tar].push_back(
                                    body_idx_src);

                                atoms_interactions_weights_tar[body_idx_tar]
                                    .push_back(wtar_influencing_cell);

                                atoms_interactions_weights_src[body_idx_tar]
                                    .push_back(wsrc2_adj1_partial);
                            }
                        }
                    }
                }

                // the most unoptimal loop, rethink later
                for (size_t j = 0; j < adj1_cells_indices.size(); j++)
                {

                    const int adj1_cell_idx = adj1_cells_indices[j];
                    const FMMCell &adj1_cell = fmm_cells[adj1_cell_idx];
                    const std::unordered_set<int> &adj2b_cells_indices_set =
                        fmm_direct_interactions_tree_.get_adjacent_cells_set(
                            adj1_cell_idx);

                    if (adj1_cell_idx != cell.index)
                    {
                        real wtar_adj1_partial =
                            fmm_weights_eval_.compute_weight_outside_cell(
                                body_tar.x, adj1_cell.center, adj1_cell.radius,
                                false);

                        if (wtar_adj1_partial > 0 &&
                            adj2b_cells_indices_set.find(adj2_cell_idx) !=
                                adj2b_cells_indices_set.end())
                        {
                            size_t bxr = 0;
                            for (const int body_idx_src :
                                 reg_bodies_cell[adj2_cell_idx])
                            {
                                int corresponding_cell =
                                    reg_bodies_from_cell_idx[adj2_cell_idx]
                                                            [bxr];

                                if (adj1_cells_indices_set.find(
                                        corresponding_cell) ==
                                        adj1_cells_indices_set.end() &&
                                    adj2b_cells_indices_set.find(
                                        corresponding_cell) ==
                                        adj2b_cells_indices_set.end())
                                {
                                    const FBody &body_src =
                                        bodies_all_[body_idx_src];

                                    real wsrc_reg2 =
                                        reg_weights_cell[adj2_cell_idx][bxr];

                                    atoms_interactions_list[body_idx_tar]
                                        .push_back(body_idx_src);

                                    atoms_interactions_weights_tar[body_idx_tar]
                                        .push_back(wtar_adj1_partial);

                                    atoms_interactions_weights_src[body_idx_tar]
                                        .push_back(wsrc_reg2);
                                }
                                bxr++;
                            }
                        }
                    }
                }
            }
        }
    }

    // std::ofstream output_file("atom_interactions_dump.txt");
    // if (!output_file.is_open())
    // {
    //     throw std::runtime_error(
    //         "Failed to open file for dumping atom interactions.");
    // }

    // // Write the data
    // for (size_t target = 0; target < bodies_all_.size(); ++target)
    // {
    //     for (size_t i = 0; i < atoms_interactions_list[target].size(); ++i)
    //     {
    //         const int body_src_idx = atoms_interactions_list[target][i];

    //         output_file << bodies_all_[body_src_idx].x << "--"
    //                     << bodies_all_[target].x << "--"
    //                     << atoms_interactions_weights_src[target][i] << "--"
    //                     << atoms_interactions_weights_tar[target][i] << "\n";
    //     }
    // }

    // // Close the file
    // output_file.close();

    // std::ofstream output_file2("actual_belongings.txt");

    // // Write the data
    // for (size_t k = 0; k < fmm_cells.size(); k++)
    // {
    //     const FMMCell &cell = fmm_cells[k];

    //     for (const int &body_idx_tar : cell.bodiesIndices)
    //     {
    //         output_file2 << bodies_all_[body_idx_tar].x << "--" <<
    //         cell.center
    //                      << "\n";
    //     }
    // }

    // // Close the file
    // output_file2.close();
}

std::vector<std::pair<gmx::RVec, gmx::real>>
gmx::fmm::FMMDirectInteractions::execute_direct_kernel()
{

    std::vector<std::pair<gmx::RVec, real>> forces_and_potentials(
        bodies_all_.size());

    for (size_t i = 0; i < bodies_all_.size(); i++)
    {
        forces_and_potentials[i] = std::make_pair(RVec(0, 0, 0), 0);
    }

    for (size_t i = 0; i < bodies_all_.size(); i++)
    {
        const FBody &body_tar = bodies_all_[i];
        const real xj = body_tar.x[0];
        const real yj = body_tar.x[1];
        const real zj = body_tar.x[2];
        real pj_effective = 0.0;
        real fxj_effective = 0.0, fyj_effective = 0.0, fzj_effective = 0.0;
        int ix = 0;
        for (auto &body_src_idx : atoms_interactions_list[i])
        {
            real pj = 0.0;
            real fxj = 0.0, fyj = 0.0, fzj = 0.0;
            gmx::fmm::FBody &asrc = bodies_all_[body_src_idx];
            const real wsrc = atoms_interactions_weights_src[i][ix];
            const real wtar = atoms_interactions_weights_tar[i][ix];

            // Compute distance differences
            const real dx = xj - asrc.x[0];
            const real dy = yj - asrc.x[1];
            const real dz = zj - asrc.x[2];

            // Compute squared distance
            real invr = dx * dx + dy * dy + dz * dz;

            invr = 1.0 / std::sqrt(invr); // Compute inverse distance

            const real qi = asrc.q * wsrc;

            real qinvr = qi * invr;
            pj = qinvr;
            qinvr = qinvr * invr * invr;

            fxj = qinvr * dx;
            fyj = qinvr * dy;
            fzj = qinvr * dz;

            pj_effective += pj * wtar;
            fxj_effective += fxj * wtar;
            fyj_effective += fyj * wtar;
            fzj_effective += fzj * wtar;

            ix++;
        }

        // Apply accumulated forces and potential to target bodies
        forces_and_potentials[i].second += pj_effective;
        forces_and_potentials[i].first[0] -= fxj_effective;
        forces_and_potentials[i].first[1] -= fyj_effective;
        forces_and_potentials[i].first[2] -= fzj_effective;
    }

    return forces_and_potentials;
}

// std::vector<std::pair<gmx::RVec, gmx::real>>
// gmx::fmm::FMMDirectInteractions:: execute_direct_kernel()
// {

//     std::vector<std::pair<gmx::RVec, real>> foreces_and_potentials(
//         bodies_all_.size());

//     for (size_t i = 0; i < bodies_all_.size(); i++)
//     {
//         foreces_and_potentials[i] = std::make_pair(RVec(0, 0, 0), 0);
//     }

//     FMMCells &fmm_cells = fmm_direct_interactions_tree_.get_cells();
//     for (size_t i = 0; i < fmm_cells.size(); i++)
//     {
//         gmx::fmm::FMMCell &cell_tar = fmm_cells[i];
//         auto &atomsTarAll = cell_tar.bodiesIndicesRegularized;
//         auto &weightsTarAll = cell_tar.weightsBodiesAll;

//         for (size_t j = 0; j < atomsTarAll.size(); j++)
//         {
//             const int atom_idx_tar = atomsTarAll[j];

//             gmx::fmm::FBody &atar = bodies_all_[atom_idx_tar];
//             const real wtar = weightsTarAll[j];
//             const real xj = atar.x[0];
//             const real yj = atar.x[1];
//             const real zj = atar.x[2];
//             real pj = 0.0;
//             real fxj = 0.0, fyj = 0.0, fzj = 0.0;

//             for (const auto &adjacent_cell_idx :
//                  fmm_direct_interactions_tree_.get_adjacent_cells(i))
//             {

//                 gmx::fmm::FMMCell &cell_src =
//                 fmm_cells[adjacent_cell_idx]; auto &atomsSrcAll =
//                 cell_src.bodiesIndicesRegularized; auto &weightsSrcAll =
//                 cell_src.weightsBodiesAll;

//                 for (size_t k = 0; k < atomsSrcAll.size(); k++)
//                 {
//                     int atom_idx_src = atomsSrcAll[k];
//                     gmx::fmm::FBody &asrc = bodies_all_[atom_idx_src];
//                     const real wsrc = weightsSrcAll[k];

//                     // Compute distance differences
//                     const real dx = xj - asrc.x[0];
//                     const real dy = yj - asrc.x[1];
//                     const real dz = zj - asrc.x[2];

//                     // Compute squared distance
//                     real invr = dx * dx + dy * dy + dz * dz;
//                     if (invr > 0) // Avoid division by zero
//                     {
//                         invr =
//                             1.0 / std::sqrt(invr); // Compute inverse
//                             distance

//                         const real qi = asrc.q * wsrc; // weight source

//                         real qinvr = qi * invr;
//                         pj += qinvr;
//                         qinvr = qinvr * invr * invr;

//                         fxj += qinvr * dx;
//                         fyj += qinvr * dy;
//                         fzj += qinvr * dz;
//                     }
//                 }
//             }

//             // Apply accumulated forces and potential to target bodies
//             foreces_and_potentials[atom_idx_tar].second += wtar * pj;
//             foreces_and_potentials[atom_idx_tar].first[0] -= wtar * fxj;
//             foreces_and_potentials[atom_idx_tar].first[1] -= wtar * fyj;
//             foreces_and_potentials[atom_idx_tar].first[2] -= wtar * fzj;
//         }
//     }

//     return foreces_and_potentials;
// }

void gmx::fmm::FMMDirectInteractions::recompute_weights()
{
    compute_weights_();
}

void gmx::fmm::FMMDirectInteractions::rebuild_and_reprocess_tree()
{
    fmm_direct_interactions_tree_.rebuild_and_reprocess_tree();
}