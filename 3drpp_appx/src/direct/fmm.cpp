#include "fmm.h"
#include <fstream>
#include <iomanip>

// assuming coordinates are simd packed (like xxx..yyy..zzz.. as per simd width)
gmx::fmm::FMMDirectInteractions::FMMDirectInteractions(const std::vector<RVec> coordinates, const std::vector<real> charges, const RVec box_center,
                                                       const real box_radius, const size_t max_depth, const real reg_alpha)
    : bodies_all_(coordinates, charges), fmm_weights_eval_(box_center, box_radius, reg_alpha),
      fmm_direct_interactions_tree_(bodies_all_, box_center, box_radius, max_depth)
{
    compute_weights_();
}

bool gmx::fmm::FMMDirectInteractions::is_point_within_radius(const RVec &point1, const RVec &point2, double radius)
{
    RVec dx = {point1[0] - point2[0], point1[1] - point2[1], point1[2] - point2[2]};
    double distance_squared = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    double radius_squared = radius * radius;
    return distance_squared <= radius_squared;
}

void gmx::fmm::FMMDirectInteractions::compute_weights_()
{
    TIME_BEGIN(compute_weights_func);

    pair_list.clear();
    pair_list.resize(bodies_all_.size());

    pair_list_bxyz_src.clear();
    pair_list_bxyz_src.resize(bodies_all_.size());

    pair_list_bxyz_tar.clear();
    pair_list_bxyz_tar.resize(bodies_all_.size());

    pair_list_w_tar.clear();
    pair_list_w_tar.resize(bodies_all_.size());

    w_per_atom.clear();
    w_per_atom.resize(bodies_all_.size());

    FMMCells &fmm_cells = fmm_direct_interactions_tree_.get_cells();

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];
            const RVec ws = fmm_weights_eval_.compute_weight_in_cell(body_tar.x, cell.center, cell.radius, false);
            w_per_atom[body_idx_tar] = ws;
        }
    }

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        size_t bidxt = 0;
        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];

            for (const int body_idx_src : cell.bodiesIndices)
            {
                if (body_idx_tar != body_idx_src)
                {
                    // const FBody &body_src = bodies_all_[body_idx_src];
                    // pair_list[body_idx_tar].push_back(body_idx_src);
                    // pair_list_w_tar[body_idx_tar].push_back(cell.weights[bidxt]);
                    // pair_list_b_src[body_idx_tar].push_back({1, 1, 1});
                }
            }

            const real interaction_region_tcell = cell.radius * 3; // one and half cell
            const real width_of_tcell = cell.radius * 2;
            const RVec wst = w_per_atom[body_idx_tar];

            for (int dz = -2; dz <= 2; dz++)
            {
                for (int dy = -2; dy <= 2; dy++)
                {
                    for (int dx = -2; dx <= 2; dx++)
                    {
                        // Skip the current cell
                        if (dx == 0 && dy == 0 && dz == 0)
                        {
                            continue;
                        }
                        const RVec neighbor_center = cell.center + RVec(dx * width_of_tcell, dy * width_of_tcell, dz * width_of_tcell);

                        const int adj_cell_idx = fmm_direct_interactions_tree_.get_neighbour_idx(neighbor_center);

                        if (adj_cell_idx != -1)
                        {

                            short num_away = 1; // Default to "one away"
                            if (std::abs(dx) == 2 || std::abs(dy) == 2 || std::abs(dz) == 2)
                            {
                                num_away = 2; // If any direction is "two away"
                            }

                            if (num_away == 1)
                            {
                                const FMMCell &adj_cell = fmm_cells[adj_cell_idx];
                                const real interaction_region_scell = adj_cell.radius * 3; // one and half cell
                                bool bxt = wst[0] != 1 ? (fabs(body_tar.x[0] - adj_cell.center[0]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scell  ? 1 : 0) : 1;

                                bool byt = wst[1] != 1 ? (fabs(body_tar.x[1] - adj_cell.center[1]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scell  ? 1 : 0) : 1;

                                bool bzt = wst[2] != 1 ? (fabs(body_tar.x[2] - adj_cell.center[2]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scell  ? 1 : 0) : 1;

                                for (const int &body_idx_src : adj_cell.bodiesIndices)
                                {
                                    const FBody &body_src = bodies_all_[body_idx_src];
                                    if (body_idx_tar != body_idx_src)
                                    {
                                        const RVec ws = w_per_atom[body_idx_src];
                                        bool bx =
                                            ws[0] < 1
                                                ? (fabs(body_src.x[0] - cell.center[0]) + fmm_weights_eval_.getRegAlpha() >= interaction_region_tcell ? 0 : 1)
                                                : 1;

                                        bool by =
                                            ws[1] < 1
                                                ? (fabs(body_src.x[1] - cell.center[1]) + fmm_weights_eval_.getRegAlpha() >= interaction_region_tcell ? 0 : 1)
                                                : 1;

                                        bool bz =
                                            ws[2] < 1
                                                ? (fabs(body_src.x[2] - cell.center[2]) + fmm_weights_eval_.getRegAlpha() >= interaction_region_tcell ? 0 : 1)
                                                : 1;

                                        pair_list[body_idx_tar].push_back(body_idx_src);
                                        pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                        pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                    }
                                }
                            }
                            else
                            {
                                // for (const int body_idx_src : boundary_bodies_idxs[adj_cell_idx])
                                // {
                                //     const FBody &body_src = bodies_all_[body_idx_src];

                                //     if (fabs(body_src.x[0] - cell.center[0]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha() &&
                                //         fabs(body_src.x[1] - cell.center[1]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha() &&
                                //         fabs(body_src.x[2] - cell.center[2]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha())
                                //     {

                                //         const RVec ws = w_per_atom[body_idx_src];
                                //         bool bx = ws[0] < 1 ? 0 : 1;
                                //         bool by = ws[1] < 1 ? 0 : 1;
                                //         bool bz = ws[2] < 1 ? 0 : 1;

                                //         pair_list_w_tar[body_idx_tar].push_back(cell.weights[bidxt]);
                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                //         pair_list_b_src[body_idx_tar].push_back({bx, by, bz});
                                //     }
                                // }
                            }
                        }
                    }
                }
            }
            bidxt++;
        }
    }
}

std::vector<std::pair<gmx::RVec, gmx::real>> gmx::fmm::FMMDirectInteractions::execute_direct_kernel()
{

    std::vector<std::pair<gmx::RVec, real>> forces_and_potentials(bodies_all_.size());

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
        const RVec wtar_ws = w_per_atom[i];

        real pj_effective = 0.0;
        real fxj_effective = 0.0, fyj_effective = 0.0, fzj_effective = 0.0;
        int ix = 0;
        for (auto &body_src_idx : pair_list[i])
        {

            gmx::fmm::FBody &asrc = bodies_all_[body_src_idx];

            const BVec bxyz_src = pair_list_bxyz_src[i][ix];
            const BVec bxyz_tar = pair_list_bxyz_tar[i][ix];

            const RVec wsrc_ws = w_per_atom[body_src_idx];

            const real dx = xj - asrc.x[0];
            const real dy = yj - asrc.x[1];
            const real dz = zj - asrc.x[2];

            // switch between w and 1-w
            real wsrc_x = bxyz_src[0] == 1 ? 1 : (dx < 0 ? 1 - wsrc_ws[0] : wsrc_ws[0]);
            real wsrc_y = bxyz_src[1] == 1 ? 1 : (dy < 0 ? 1 - wsrc_ws[1] : wsrc_ws[1]);
            real wsrc_z = bxyz_src[2] == 1 ? 1 : (dz < 0 ? 1 - wsrc_ws[2] : wsrc_ws[2]);
            real wsrc = wsrc_x * wsrc_y * wsrc_z;

            real wtar_x = bxyz_tar[0] == 1 ? 1 : (dx > 0 ? 1 - wtar_ws[0] : wtar_ws[0]);
            real wtar_y = bxyz_tar[1] == 1 ? 1 : (dy > 0 ? 1 - wtar_ws[1] : wtar_ws[1]);
            real wtar_z = bxyz_tar[2] == 1 ? 1 : (dz > 0 ? 1 - wtar_ws[2] : wtar_ws[2]);
            const real wtar = wtar_x * wtar_y * wtar_z;

            std::cout << asrc.x << "--" << body_tar.x << "--" << wsrc << "--" << wtar <<  "--" << wtar_ws[0] << "," <<wtar_ws[1] << "," << wtar_ws[2] << std::endl;

            real pj = 0.0;
            real fxj = 0.0, fyj = 0.0, fzj = 0.0;

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

void gmx::fmm::FMMDirectInteractions::recompute_weights() { compute_weights_(); }

void gmx::fmm::FMMDirectInteractions::rebuild_and_reprocess_tree() { fmm_direct_interactions_tree_.rebuild_and_reprocess_tree(); }
