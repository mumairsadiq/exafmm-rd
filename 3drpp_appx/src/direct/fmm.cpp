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
    TIME_BEGIN(compute_weights_func);

    pair_list.clear();
    pair_list.resize(bodies_all_.size());

    pair_list_b_src.clear();
    pair_list_b_src.resize(bodies_all_.size());

    pair_list_b_tar.clear();
    pair_list_b_tar.resize(bodies_all_.size());

    w_per_atom.clear();
    w_per_atom.resize(bodies_all_.size());

    FMMCells &fmm_cells = fmm_direct_interactions_tree_.get_cells();
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];
            const RVec ws = fmm_weights_eval_.compute_weight_in_cell(
                body_tar.x, cell.center, cell.radius, false);

            w_per_atom[body_idx_tar] = ws;

            bool bxt = ws[0] != 1 ? 1 : 0;
            bool byt = ws[0] != 1 ? 1 : 0;
            bool bzt = ws[0] != 1 ? 1 : 0;
            pair_list_b_tar[body_idx_tar] = {bxt, byt, bzt};
        }
    }

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        const auto &adj1_cell =
            fmm_direct_interactions_tree_.get_adjacent_cells(k);
        size_t bidxt = 0;
        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];

            // it is possible that target weight is split into current
            // cell and first adjacent cell, and source may only
            // interact from adjacent cell
            // real target_weight_within_cell =
            //     fmm_weights_eval_.compute_weight_within_cell(
            //         body_tar.x, cell.center, cell.radius, false);

            // all bodies within same cell

            for (const int body_idx_src : cell.bodiesIndices)
            {
                if (body_idx_tar != body_idx_src)
                {
                    const FBody &body_src = bodies_all_[body_idx_src];
                    // pair_list[body_idx_tar].push_back(body_idx_src);
                    // pair_list_b_src[body_idx_tar].push_back({1, 1, 1});
                }
            }
            bidxt++;

            for (size_t j = 0; j < adj1_cell.size(); j++)
            {
                const int adj_cell_idx = adj1_cell[j];
                if (cell.index != adj_cell_idx)
                {
                    const FMMCell &adj_cell = fmm_cells[adj_cell_idx];
                    // checkout for regularization boundary bodies from adjacent
                    // cells having weight less than one in their cell
                    for (const int &body_idx_src : adj_cell.bodiesIndices)
                    {

                        const RVec ws = w_per_atom[body_idx_src];
                        bool bx = ws[0] < 1 ? 0 : 1;
                        bool by = ws[1] < 1 ? 0 : 1;
                        bool bz = ws[2] < 1 ? 0 : 1;
                        pair_list[body_idx_tar].push_back(body_idx_src);
                        pair_list_b_src[body_idx_tar].push_back({bx, by, bz});
                    }
                }
            }
        }
    }
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
        for (auto &body_src_idx : pair_list[i])
        {
            const RVec wtar_ws = w_per_atom[i];
            const BVec wtar_b = pair_list_b_tar[i];

            real wtar = wtar_ws[0] * wtar_ws[1] * wtar_ws[2];
            gmx::fmm::FBody &asrc = bodies_all_[body_src_idx];

            const BVec wsrc_b = pair_list_b_src[i][ix];
            const RVec wsrc_ws = w_per_atom[body_src_idx];

            // real wsrc_own = wsrc_ws[0] * wsrc_ws[1] * wsrc_ws[2];

            const real dx = xj - asrc.x[0];
            const real dy = yj - asrc.x[1];
            const real dz = zj - asrc.x[2];

            real wsrc_x =
                wsrc_b[0] == 1 ? 1 : (dx < 0 ? wsrc_ws[0] : 1 - wsrc_ws[0]);
            real wsrc_y =
                wsrc_b[1] == 1 ? 1 : (dy < 0 ?  wsrc_ws[1] :1 - wsrc_ws[1]);
            real wsrc_z =
                wsrc_b[2] == 1 ? 1 : (dz < 0 ? wsrc_ws[2] : 1 - wsrc_ws[2]);

            real pj = 0.0;
            real fxj = 0.0, fyj = 0.0, fzj = 0.0;

            real wsrc = wsrc_x * wsrc_y * wsrc_z;

            std::cout << asrc.x << "--" << body_tar.x << "--" << wsrc << "--"
                      << wtar << "\n";

            // const real wsrcx = pair_list_b_src[i][ix][0] != 0
            //                        ? w_per_atom[body_src_idx][0]
            //                        : 1;
            // const real wsrcy = pair_list_b_src[i][ix][1] != 0
            //                        ? w_per_atom[body_src_idx][1]
            //                        : 1;
            // const real wsrcz = pair_list_b_src[i][ix][2] != 0
            //                        ? w_per_atom[body_src_idx][2]
            //                        : 1;
            // const real wsrc = wsrcx * wsrcy * wsrcz;

            // Compute distance differences

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