#include "fmm.h"
#include <fstream>

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
                    atoms_interactions_list[body_idx_tar].push_back(
                        body_idx_src);
                    atoms_interactions_weights_src[body_idx_tar].push_back(1);
                    atoms_interactions_weights_tar[body_idx_tar].push_back(1);
                }
            }
        }
    }

    std::vector<FPIndices> boundary_bodies_idxs(fmm_cells.size());
    std::vector<std::vector<real>> boundary_bodies_weights(fmm_cells.size());
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
                boundary_bodies_weights[k].push_back(w);
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

        const auto &aCells =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);
        std::vector<int> reg_body_idxs;
        std::vector<int> reg_body_frm_cell_idxs;

        // find out bodies from adjacent cells that are influencing this cell
        for (size_t j = 0; j < aCells.size(); j++)
        {
            const int adj_cell_idx = aCells[j];
            if (cell.index != adj_cell_idx)
            {
                // checkout for regularization boundary bodies from adjacent
                // cells having weight less than one in their cell
                for (const int &body_idx : boundary_bodies_idxs[adj_cell_idx])
                {

                    const FBody &body = bodies_all_[body_idx];
                    const RVec dx = body.x - cell.center;
                    const RVec w_inter =
                        fmm_weights_eval_.compute_w_xyz(dx, cell.radius);
                    const real w = w_inter[0] * w_inter[1] * w_inter[2];
                    if (w > 0)
                    {
                        reg_body_frm_cell_idxs.push_back(adj_cell_idx);
                        reg_body_idxs.push_back(body_idx);
                    }
                }
            }
        }

        // Compute weights for bodies from adjacent cells
        for (size_t i = 0; i < reg_body_idxs.size(); i++)
        {
            const int &body_idx = reg_body_idxs[i];
            const FBody &body = bodies_all_[body_idx];

            const real w = fmm_weights_eval_.compute_weight_outside_cell(
                body.x, cell.center, cell.radius, false);

            if (w > 0)
            {
                reg_bodies_cell[k].push_back(body_idx);
                reg_weights_cell[k].push_back(w);
                reg_bodies_from_cell_idx[k].push_back(
                    reg_body_frm_cell_idxs[i]);
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

                    const real w2 =
                        fmm_weights_eval_.compute_weight_outside_cell(
                            body_tar.x, adj_cell.center, adj_cell.radius,
                            false);
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

                        atoms_interactions_list[body_idx_tar].push_back(
                            body_idx_src);

                        atoms_interactions_weights_tar[body_idx_tar].push_back(
                            bodies_weights_to_adj_cells[body_idx_tar]
                                                       [adj1_cell_idx]);

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

                if (adj1_cell_idx != cell.index)
                {
                    const real wtar_in_original =
                        fmm_weights_eval_.compute_weight_within_cell(
                            body_tar.x, cell.center, cell.radius, false);

                    const real wtar_partial_in_adjacent =
                        fmm_weights_eval_.compute_weight_outside_cell(
                            body_tar.x, adj1_cell.center, adj1_cell.radius,
                            false);

                    const real wtar_full_in_adjacent =
                        bodies_weights_to_adj_cells[body_idx_tar]
                                                   [adj1_cell_idx];

                    if (wtar_full_in_adjacent > 0)
                    {
                        const auto &adj2_cells_indices =
                            fmm_direct_interactions_tree_.get_adjacent_cells(
                                adj1_cell_idx);
                        real wtar_rem =
                            wtar_full_in_adjacent - wtar_partial_in_adjacent;

                        for (size_t l = 0; l < adj2_cells_indices.size(); l++)
                        {
                            const int adj2_cell_idx = adj2_cells_indices[l];

                            if (adj1_cells_indices_set.find(adj2_cell_idx) ==
                                adj1_cells_indices_set.end())
                            {
                                const FMMCell &adj2_cell =
                                    fmm_cells[adj2_cell_idx];

                                size_t bidxs2 = 0;
                                for (const int body_idx_src :
                                     adj2_cell.bodiesIndices)
                                {
                                    const FBody &body_src =
                                        bodies_all_[body_idx_src];

                                    real wsrc_full_in_adj1 =
                                        bodies_weights_to_adj_cells
                                            [body_idx_src][adj1_cell_idx];

                                    const real wsrc_partial_in_adjacent =
                                        fmm_weights_eval_
                                            .compute_weight_outside_cell(
                                                body_src.x, adj1_cell.center,
                                                adj1_cell.radius, false);

                                    if (wsrc_full_in_adj1 > 0 &&
                                        wtar_partial_in_adjacent > 0)
                                    {
                                        atoms_interactions_list[body_idx_tar]
                                            .push_back(body_idx_src);

                                        atoms_interactions_weights_tar
                                            [body_idx_tar]
                                                .push_back(
                                                    wtar_partial_in_adjacent);

                                        atoms_interactions_weights_src
                                            [body_idx_tar]
                                                .push_back(wsrc_full_in_adj1);
                                    }
                                    
                                    if (is_part_of_reg_region[adj2_cell_idx]
                                                             [bidxs2])
                                    {

                                        if (wtar_rem > 0 &&
                                            wsrc_partial_in_adjacent > 0)
                                        {
                                            atoms_interactions_list
                                                [body_idx_tar]
                                                    .push_back(body_idx_src);

                                            atoms_interactions_weights_tar
                                                [body_idx_tar]
                                                    .push_back(wtar_rem);

                                            atoms_interactions_weights_src
                                                [body_idx_tar]
                                                    .push_back(
                                                        wsrc_partial_in_adjacent);
                                        }
                                    }

                                    bidxs2++;
                                }
                            }

                            if (wtar_partial_in_adjacent > 0)
                            {
                                int bx = 0;
                                for (const int body_idx_src :
                                     reg_bodies_cell[adj2_cell_idx])
                                {
                                    int corresponding_cell =
                                        reg_bodies_from_cell_idx[adj2_cell_idx]
                                                                [bx];

                                    if (adj1_cells_indices_set.find(
                                            corresponding_cell) ==
                                        adj1_cells_indices_set.end())
                                    {
                                        const FBody &body_src =
                                            bodies_all_[body_idx_src];

                                        real wsrc =
                                            reg_weights_cell[adj2_cell_idx][bx];

                                        atoms_interactions_list[body_idx_tar]
                                            .push_back(body_idx_src);

                                        atoms_interactions_weights_tar
                                            [body_idx_tar]
                                                .push_back(
                                                    wtar_partial_in_adjacent);

                                        atoms_interactions_weights_src
                                            [body_idx_tar]
                                                .push_back(wsrc);
                                    }

                                    bx++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::ofstream output_file("atom_interactions_dump.txt");
    if (!output_file.is_open())
    {
        throw std::runtime_error(
            "Failed to open file for dumping atom interactions.");
    }

    // Write the data
    for (size_t target = 0; target < bodies_all_.size(); ++target)
    {

        // Write the target body index once
        // output_file << "TargetBody: " << bodies_all_[target].x << "\n";

        // Write all its source bodies and weights
        for (size_t i = 0; i < atoms_interactions_list[target].size(); ++i)
        {
            // output_file << "\tSourceBody: "
            //             <<
            //             bodies_all_[atoms_interactions_list[target][i]].x
            //             << ", Weight: " <<
            //             atoms_interactions_weights[target][i]
            //             << "\n";

            // if (!(atoms_interactions_weights_src[target][i] == 1 &&
            //       atoms_interactions_weights_tar[target][i] == 1))
            {
                const int body_src_idx = atoms_interactions_list[target][i];
                output_file << bodies_all_[body_src_idx].x << "--"
                            << bodies_all_[target].x << "--"
                            << atoms_interactions_weights_src[target][i] << "--"
                            << atoms_interactions_weights_tar[target][i]
                            << "\n";
            }
        }
    }

    // Close the file
    output_file.close();

    std::ofstream output_file2("actual_belongings.txt");

    // Write the data
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];

        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            output_file2 << bodies_all_[body_idx_tar].x << "--" << cell.center
                         << "\n";
        }
    }

    // Close the file
    output_file2.close();
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