#include "fmm.h"
#include <fstream>

// assuming coordinates are simd packed (like xxx..yyy..zzz.. as per simd width)
gmx::fmm::FMMDirectInteractions::FMMDirectInteractions(
    const std::vector<RVec> coordinates, const std::vector<real> charges,
    const RVec box_center, const real box_radius,
    const size_t max_particles_per_cell, const real reg_alpha)
    : bodies_all_(coordinates, charges), fmm_weights_eval_(reg_alpha),
      fmm_direct_interactions_tree_(bodies_all_, box_center, box_radius,
                                    max_particles_per_cell)
{
    compute_weights_();
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

    std::vector<FPIndices> boundary_bodies_idxs(fmm_cells.size());

    // find weights for existing boundary bodies within the same cell
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        for (const int &body_idx : cell.bodiesIndices)
        {
            const FBody &body = bodies_all_[body_idx];
            const RVec dx = body.x - cell.center;

            RVec dx_simcenter_inter =
                body.x - fmm_direct_interactions_tree_.get_box_center();
            dx_simcenter_inter[0] = fabs(dx_simcenter_inter[0]);
            dx_simcenter_inter[1] = fabs(dx_simcenter_inter[1]);
            dx_simcenter_inter[2] = fabs(dx_simcenter_inter[2]);
            const RVec dx_simcenter(
                dx_simcenter_inter[0] + fmm_weights_eval_.getRegAlpha(),
                dx_simcenter_inter[1] + fmm_weights_eval_.getRegAlpha(),
                dx_simcenter_inter[2] + fmm_weights_eval_.getRegAlpha());

            RVec ws = fmm_weights_eval_.compute_w_xyz(dx, cell.radius);

            for (int d = 0; d <= 2; d++)
            {
                if (dx_simcenter[d] >=
                    fmm_direct_interactions_tree_.get_box_radius())
                {
                    ws[d] = 1;
                }
            }
            const real w = ws[0] * ws[1] * ws[2];

            if (w < 1)
            {
                boundary_bodies_idxs[k].push_back(body_idx);
            }
        }
    }

    std::vector<std::vector<int>> bodies_regularized(fmm_cells.size());
    std::vector<std::vector<real>> weights(fmm_cells.size());

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {

        const FMMCell &cell = fmm_cells[k];

        const auto &aCells =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);
        std::vector<int> reg_body_idxs;

        // find out bodies from adjacent cells that are influencing this cell
        for (size_t j = 0; j < aCells.size(); j++)
        {
            const int adj_cell_idx = aCells[j];
            if (cell.index != fmm_cells[adj_cell_idx].index)
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
            RVec body_x = body.x;
            const RVec dx = body_x - cell.center;
            RVec dx_simcenter_inter =
                body_x - fmm_direct_interactions_tree_.get_box_center();
            dx_simcenter_inter[0] = fabs(dx_simcenter_inter[0]);
            dx_simcenter_inter[1] = fabs(dx_simcenter_inter[1]);
            dx_simcenter_inter[2] = fabs(dx_simcenter_inter[2]);
            const RVec dx_simcenter(
                dx_simcenter_inter[0] + fmm_weights_eval_.getRegAlpha(),
                dx_simcenter_inter[1] + fmm_weights_eval_.getRegAlpha(),
                dx_simcenter_inter[2] + fmm_weights_eval_.getRegAlpha());
            RVec ws = fmm_weights_eval_.compute_w_xyz(dx, cell.radius);

            for (int d = 0; d <= 2; d++)
            {
                if (dx_simcenter[d] >=
                    fmm_direct_interactions_tree_.get_box_radius())
                {
                    ws[d] = 1;
                }
            }
            const real w = ws[0] * ws[1] * ws[2];
            if (w > 0)
            {
                bodies_regularized[k].push_back(body_idx);
                weights[k].push_back(w);
            }
        }
    }
    // find weights for existing bodies within the same cell
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];
            const auto &adj_cells =
                fmm_direct_interactions_tree_.get_adjacent_cells(k);

            for (size_t j = 0; j < adj_cells.size(); j++)
            {
                const int adj_cell_idx = adj_cells[j];
                const FMMCell &adj_cell = fmm_cells[adj_cell_idx];

                for (const int body_idx_src : adj_cell.bodiesIndices)
                {
                    if (body_idx_tar != body_idx_src)
                    {
                        atoms_interactions_list[body_idx_tar].push_back(
                            body_idx_src);
                        atoms_interactions_weights_src[body_idx_tar].push_back(
                            1);
                        atoms_interactions_weights_tar[body_idx_tar].push_back(
                            1);
                    }
                }
            }
        }
    }

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];
            const auto &adj_cells =
                fmm_direct_interactions_tree_.get_adjacent_cells(k);
            for (size_t j = 0; j < adj_cells.size(); j++)
            {
                const int adj_cell_idx = adj_cells[j];
                const FMMCell &adj_cell = fmm_cells[adj_cell_idx];
                if (cell.index != adj_cell.index)
                {
                    int ix = 0;
                    for (const int body_idx_src :
                         bodies_regularized[adj_cell_idx])
                    {
                        if (body_idx_tar != body_idx_src)
                        {
                            bool is_interaction_by_w_one = false;
                            for (int tt = 0;
                                 tt <
                                 atoms_interactions_list[body_idx_tar].size();
                                 tt++)
                            {
                                if (body_idx_src ==
                                        atoms_interactions_list[body_idx_tar]
                                                               [tt] &&
                                    atoms_interactions_weights_src[body_idx_tar]
                                                                  [tt] == 1 &&
                                    atoms_interactions_weights_tar[body_idx_tar]
                                                                  [tt] == 1)
                                {
                                    is_interaction_by_w_one = true;
                                    break;
                                }
                            }

                            if (is_interaction_by_w_one == false)
                            {
                                atoms_interactions_list[body_idx_tar].push_back(
                                    body_idx_src);
                                atoms_interactions_weights_src[body_idx_tar]
                                    .push_back(weights[adj_cell_idx][ix]);
                                atoms_interactions_weights_tar[body_idx_tar]
                                    .push_back(1);

                                atoms_interactions_list[body_idx_src].push_back(
                                    body_idx_tar);
                                atoms_interactions_weights_tar[body_idx_src]
                                    .push_back(weights[adj_cell_idx][ix]);
                                atoms_interactions_weights_src[body_idx_src]
                                    .push_back(1);
                            }
                        }
                        ix++;
                    }
                }
            }
        }
    }

    weights.clear();
    bodies_regularized.clear();

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

            if (!(atoms_interactions_weights_src[target][i] == 1 &&
                atoms_interactions_weights_tar[target][i] == 1))
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
