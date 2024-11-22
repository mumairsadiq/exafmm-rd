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
    std::vector<std::vector<bool>> is_part_of_reg_region(fmm_cells.size());

    // find weights for existing boundary bodies within the same cell
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        is_part_of_reg_region[k].resize(cell.bodiesIndices.size());
        int bidx = 0;
        for (const int &body_idx : cell.bodiesIndices)
        {
            const FBody &body = bodies_all_[body_idx];
            const real w = fmm_weights_eval_.compute_weight(body.x, cell.center,
                                                            cell.radius, false);

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

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        const auto &adj1_cells =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);

        // rethink later if it useful to create this massive data structure for
        // precomputation of weights third dimension can be added in that too to
        // totally precompute everything but such a massive space requirment
        // could be bottleneck so better to do everything on the fly for time
        // being alternatively w_total if that can be computed somehow on the
        // fly in faster way that would solve this

        /*
        real body_tar_weights_by_src_cell[adj1_cells.size()]
                                         [cell.bodiesIndices.size()];
        std::vector<std::vector<real>> body_tar_weights_by_src_cell(
            adj1_cells.size(), std::vector<real>(cell.bodiesIndices.size()));

        for (size_t j = 0; j < adj1_cells.size(); j++)
        {
            size_t bidxt = 0;
            const int adj_cell_idx = adj1_cells[j];
            const FMMCell &adj_cell = fmm_cells[adj_cell_idx];
            for (const int &body_idx_tar : cell.bodiesIndices)
            {
                const FBody &body_tar = bodies_all_[body_idx_tar];
                body_tar_weights_by_src_cell[j][bidxt] = 1;
                if (is_part_of_reg_region[k][bidxt] == true)
                {
                    const real w2 = fmm_weights_eval_.compute_weight(
                        body_tar.x, adj_cell.center, adj_cell.radius, false);
                    if (w2 == 0)
                    {
                        // needs optimization, this is brute force
                        // approach to compute w_total
                        real w_total = 0;
                        const std::unordered_set<int> &adj_cells_src_set =
                            fmm_direct_interactions_tree_
                                .get_adjacent_cells_set(adj_cell_idx);
                        for (size_t l = 0; l < adj1_cells.size(); l++)
                        {
                            if (adj_cells_src_set.find(adj1_cells[l]) !=
                                adj_cells_src_set.end())
                            {
                                const real w_temp =
                                    fmm_weights_eval_.compute_weight(
                                        body_tar.x,
                                        fmm_cells[adj1_cells[l]].center,
                                        fmm_cells[adj1_cells[l]].radius,
                                        false);
                                w_total += w_temp;
                            }
                        }
                        body_tar_weights_by_src_cell[j][bidxt] = w_total;
                    }
                }
            }
        }
        */

        int bidxt = 0;
        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];

            for (size_t j = 0; j < adj1_cells.size(); j++)
            {
                const int adj1_cell_idx = adj1_cells[j];
                const FMMCell &adj1_cell = fmm_cells[adj1_cell_idx];

                if (adj1_cell.index != cell.index)
                {

                    int bidxs = 0;
                    for (const int body_idx_src : adj1_cell.bodiesIndices)
                    {
                        const FBody &body_src = bodies_all_[body_idx_src];
                        atoms_interactions_list[body_idx_tar].push_back(
                            body_idx_src);

                        if (is_part_of_reg_region[k][bidxt] == true)
                        {
                            const real w2 = fmm_weights_eval_.compute_weight(
                                body_tar.x, adj1_cell.center, adj1_cell.radius,
                                false);

                            if (w2 != 0)
                            {
                                atoms_interactions_weights_tar[body_idx_tar]
                                    .push_back(1);
                            }
                            else
                            {
                                // needs optimization, this is brute force
                                // approach to compute w_total
                                real w_total = 0;
                                const std::unordered_set<int>
                                    &adj2_cells_indices_set =
                                        fmm_direct_interactions_tree_
                                            .get_adjacent_cells_set(
                                                adj1_cell_idx);
                                for (size_t l = 0; l < adj1_cells.size(); l++)
                                {
                                    if (adj2_cells_indices_set.find(
                                            adj1_cells[l]) !=
                                        adj2_cells_indices_set.end())
                                    {
                                        const real w_temp =
                                            fmm_weights_eval_.compute_weight(
                                                body_tar.x,
                                                fmm_cells[adj1_cells[l]].center,
                                                fmm_cells[adj1_cells[l]].radius,
                                                false);
                                        w_total += w_temp;
                                    }
                                }
                                atoms_interactions_weights_tar[body_idx_tar]
                                    .push_back(w_total);
                            }
                        }
                        else
                        {
                            atoms_interactions_weights_tar[body_idx_tar]
                                .push_back(1);
                        }

                        if (is_part_of_reg_region[adj1_cell_idx][bidxs] == true)
                        {
                            const real w1 = fmm_weights_eval_.compute_weight(
                                body_src.x, cell.center, cell.radius, false);

                            if (w1 != 0)
                            {
                                atoms_interactions_weights_src[body_idx_tar]
                                    .push_back(1);
                            }
                            else
                            {
                                // needs optimization, this is brute force
                                // approach to compute w_total
                                real w_total = 0;
                                const std::unordered_set<int>
                                    &adj1_cells_indices_set =
                                        fmm_direct_interactions_tree_
                                            .get_adjacent_cells_set(k);

                                const FPIndices &adj2_cells =
                                    fmm_direct_interactions_tree_
                                        .get_adjacent_cells(adj1_cell_idx);

                                for (size_t l = 0; l < adj2_cells.size(); l++)
                                {
                                    if (adj1_cells_indices_set.find(
                                            adj2_cells[l]) !=
                                        adj1_cells_indices_set.end())
                                    {
                                        const real w_temp =
                                            fmm_weights_eval_.compute_weight(
                                                body_src.x,
                                                fmm_cells[adj2_cells[l]].center,
                                                fmm_cells[adj2_cells[l]].radius,
                                                false);
                                        w_total += w_temp;
                                    }
                                }
                                atoms_interactions_weights_src[body_idx_tar]
                                    .push_back(w_total);
                            }
                        }
                        else
                        {
                            atoms_interactions_weights_src[body_idx_tar]
                                .push_back(1);
                        }

                        bidxs++;
                    }
                }
            }
            bidxt++;
        }
    }

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        const auto &adj1_cells =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);
        const std::unordered_set<int> &adj1_cells_set =
            fmm_direct_interactions_tree_.get_adjacent_cells_set(k);

        int bidxt = 0;
        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];

            for (size_t j = 0; j < adj1_cells.size(); j++)
            {
                const int adj1_cell_idx = adj1_cells[j];
                const FMMCell &adj1_cell = fmm_cells[adj1_cell_idx];

                if (adj1_cell.index != cell.index)
                {

                    real wtar = 1;
                    if (is_part_of_reg_region[k][bidxt] == true)
                    {
                        // needs optimization, this is brute force
                        // approach to compute w_total
                        real w_total = 0;
                        const std::unordered_set<int> &adj2_cells_set =
                            fmm_direct_interactions_tree_
                                .get_adjacent_cells_set(adj1_cell_idx);
                        for (size_t l = 0; l < adj1_cells.size(); l++)
                        {
                            if (adj2_cells_set.find(adj1_cells[l]) !=
                                adj2_cells_set.end())
                            {
                                const real w_temp =
                                    fmm_weights_eval_.compute_weight(
                                        body_tar.x,
                                        fmm_cells[adj1_cells[l]].center,
                                        fmm_cells[adj1_cells[l]].radius, false);
                                w_total += w_temp;
                            }
                        }
                        wtar = w_total;
                    }

                    const FPIndices &adj2_cells_indices =
                        fmm_direct_interactions_tree_.get_adjacent_cells(
                            adj1_cell_idx);
                    for (size_t l = 0; l < adj2_cells_indices.size(); l++)
                    {
                        if (adj1_cells_set.find(adj2_cells_indices[l]) ==
                            adj1_cells_set.end())
                        {
                            const FMMCell &adj2_cell =
                                fmm_cells[adj2_cells_indices[l]];

                            for (size_t v = 0;
                                 v <
                                 boundary_bodies_idxs[adj2_cell.index].size();
                                 v++)
                            {
                                const FBody &body_2_away = bodies_all_
                                    [boundary_bodies_idxs[adj2_cell.index][v]];
                                const real w = fmm_weights_eval_.compute_weight(
                                    body_2_away.x, adj1_cell.center,
                                    adj1_cell.radius, false);

                                if (w > 0)
                                {
                                    atoms_interactions_list[body_idx_tar]
                                        .push_back(body_2_away.idx);
                                    atoms_interactions_weights_tar[body_idx_tar]
                                        .push_back(wtar);
                                    atoms_interactions_weights_src[body_idx_tar]
                                        .push_back(w);

                                    real wtar_rev =
                                        fmm_weights_eval_.compute_weight(
                                            body_tar.x, cell.center,
                                            cell.radius, false);

                                    atoms_interactions_list[body_2_away.idx]
                                        .push_back(body_idx_tar);
                                    atoms_interactions_weights_tar[body_2_away
                                                                       .idx]
                                        .push_back(w);
                                    atoms_interactions_weights_src[body_2_away
                                                                       .idx]
                                        .push_back(wtar_rev);
                                }
                            }
                        }
                    }
                }
            }
            bidxt++;
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
            //             << bodies_all_[atoms_interactions_list[target][i]].x
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

//                 gmx::fmm::FMMCell &cell_src = fmm_cells[adjacent_cell_idx];
//                 auto &atomsSrcAll = cell_src.bodiesIndicesRegularized;
//                 auto &weightsSrcAll = cell_src.weightsBodiesAll;

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