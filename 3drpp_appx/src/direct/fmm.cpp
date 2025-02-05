#include "fmm.h"
#include <array>
#include <fstream>
#include <iomanip>

gmx::fmm::FMMDirectInteractions::FMMDirectInteractions(const std::vector<RVec> coordinates, const std::vector<real> charges, const RVec box_center, const real box_radius,
                                                       const size_t max_depth, const real reg_alpha)
    : bodies_all_(coordinates, charges), fmm_weights_eval_(box_center, box_radius, reg_alpha), fmm_direct_interactions_tree_(bodies_all_, box_center, box_radius, max_depth)
{
    TIME_BEGIN(weights_time_new);
    compute_weights_();
    TIME_END(weights_time_new);
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
    const real reg_alpha = fmm_weights_eval_.getRegAlpha();

    pair_list.clear();
    pair_list.resize(bodies_all_.size());

    w_per_atom.clear();
    w_per_atom.resize(bodies_all_.size());

    FMMCells &fmm_cells = fmm_direct_interactions_tree_.get_cells();

    std::vector<std::unordered_map<int, FixedPairListMap>> pair_list_aux(bodies_all_.size());
    std::vector<FPIndices> boundary_bodies_idxs(fmm_cells.size());
    std::vector<bool> is_reg_body(bodies_all_.size(), false);

    std::vector<FPIndices> bodiesIndicesReg(fmm_cells.size());
    std::vector<std::vector<WeightFlags>> w_flags(fmm_cells.size());
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        for (const int &body_idx : cell.bodiesIndices)
        {
            const FBody &body = bodies_all_[body_idx];
            const RVec ws = fmm_weights_eval_.compute_weight_in_cell(body.x, cell.center, cell.radius, false);
            const real w = ws[0] * ws[1] * ws[2];
            w_per_atom[body_idx] = ws;

            if (w < 1)
            {
                boundary_bodies_idxs[k].push_back(body_idx);
                is_reg_body[body_idx] = true;
            }
            const bool is_x_fully_in = ws[0] == 1;
            const bool is_y_fully_in = ws[1] == 1;
            const bool is_z_fully_in = ws[2] == 1;
            bodiesIndicesReg[k].push_back(body_idx);
            w_flags[k].emplace_back(is_x_fully_in, is_y_fully_in, is_z_fully_in, 1, 1, 1);
        }
    }

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        FMMCell &cell = fmm_cells[k];
        const real width_of_tcell = cell.radius * 2;

        for (int dz = -1; dz <= 1; dz++)
        {
            for (int dy = -1; dy <= 1; dy++)
            {
                for (int dx = -1; dx <= 1; dx++)
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
                        for (const int &body_idx : boundary_bodies_idxs[adj_cell_idx])
                        {
                            const RVec ws = w_per_atom[body_idx];
                            const bool is_x_fully_in = ws[0] == 1;
                            const bool is_y_fully_in = ws[1] == 1;
                            const bool is_z_fully_in = ws[2] == 1;

                            const FBody &body = bodies_all_[body_idx];
                            const RVec dx = body.x - cell.center;

                            const real dist_x = fabs(dx[0]);
                            const real dist_y = fabs(dx[1]);
                            const real dist_z = fabs(dx[2]);

                            const bool is_dist_x_in_range = dist_x <= cell.radius + reg_alpha;
                            const bool is_dist_y_in_range = dist_y <= cell.radius + reg_alpha;
                            const bool is_dist_z_in_range = dist_z <= cell.radius + reg_alpha;
                            if (is_dist_x_in_range && is_dist_y_in_range && is_dist_z_in_range)
                            {
                                bodiesIndicesReg[k].push_back(body_idx);
                                w_flags[k].emplace_back(is_x_fully_in, is_y_fully_in, is_z_fully_in, dist_x <= cell.radius, dist_y <= cell.radius, dist_z <= cell.radius);
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
        const real region_of_tcell = cell.radius * 3;
        const real red_region_tcell = region_of_tcell - reg_alpha;
        const real width_of_tcell = cell.radius * 2;
        size_t bidxt = 0;
        for (const int &body_idx_tar : bodiesIndicesReg[k])
        {
            const WeightFlags entry_tar_btar(w_flags[k][bidxt].bx, w_flags[k][bidxt].by, w_flags[k][bidxt].bz, w_flags[k][bidxt].x_in, w_flags[k][bidxt].y_in, w_flags[k][bidxt].z_in);
            for (const int body_idx_src : cell.bodiesIndices)
            {
                WeightFlags entry_tar = entry_tar_btar;                                    
                if (body_idx_tar != body_idx_src)
                {
                    WeightFlags entry_src(true, true, true, true, true, true);

                    auto &entry_map = pair_list_aux[body_idx_tar][body_idx_src];
                    WeightFlags *entry_pre = entry_map.find(entry_src);
                    if (entry_pre)
                    {
                        entry_pre->bx |= (entry_tar.x_in ^ entry_pre->x_in);
                        entry_pre->by |= (entry_tar.y_in ^ entry_pre->y_in);
                        entry_pre->bz |= (entry_tar.z_in ^ entry_pre->z_in);
                        entry_tar = *entry_pre;
                    }
                    entry_map.insert(entry_src, entry_tar);
                }
            }
            bidxt++;
        }

        

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
                        const FMMCell &adj_cell = fmm_cells[adj_cell_idx];
                        

                        const real dist_x = fabs(adj_cell.center[0] - cell.center[0]);
                        const real dist_y = fabs(adj_cell.center[1] - cell.center[1]);
                        const real dist_z = fabs(adj_cell.center[2] - cell.center[2]);

                        const bool is_distx_largest = dist_x >= dist_y && dist_x >= dist_z;
                        const bool is_disty_largest = dist_y >= dist_x && dist_y >= dist_z;
                        const bool is_distz_largest = dist_z >= dist_x && dist_z >= dist_y;

                        const real interaction_region_x = dist_x - adj_cell.radius;
                        const real interaction_region_y = dist_y - adj_cell.radius;
                        const real interaction_region_z = dist_z - adj_cell.radius;

                        const real ext_regionx_scell = interaction_region_x + reg_alpha;
                        const real ext_regiony_scell = interaction_region_y + reg_alpha;
                        const real ext_regionz_scell = interaction_region_z + reg_alpha;

                        short num_away = 1; // Default to "one away"

                        if (std::abs(dx) == 2 || std::abs(dy) == 2 || std::abs(dz) == 2)
                        {
                            num_away = 2; // If any direction is "two away"
                        }

                        bidxt = 0;
                        for (const int &body_idx_tar : bodiesIndicesReg[k])
                        {
                            const WeightFlags entry_tar_btar(w_flags[k][bidxt].bx, w_flags[k][bidxt].by, w_flags[k][bidxt].bz, w_flags[k][bidxt].x_in, w_flags[k][bidxt].y_in, w_flags[k][bidxt].z_in);
                            for (const int &body_idx_src : adj_cell.bodiesIndices)
                            {
                                WeightFlags entry_tar = entry_tar_btar;
                                const FBody &body_src = bodies_all_[body_idx_src];
                                
                                const RVec ws = w_per_atom[body_idx_src];

                                const real dist_x_bd = fabs(body_src.x[0] - cell.center[0]);
                                const real dist_y_bd = fabs(body_src.x[1] - cell.center[1]);
                                const real dist_z_bd = fabs(body_src.x[2] - cell.center[2]);
                                
                                
                                const bool bx = (ws[0] == 1) || (dist_x_bd <= red_region_tcell);
                                const bool by = (ws[1] == 1) || (dist_y_bd <= red_region_tcell);
                                const bool bz = (ws[2] == 1) || (dist_z_bd <= red_region_tcell);

                                const bool dist_x_bd_in_region = dist_x_bd <= ext_regionx_scell;
                                const bool dist_y_bd_in_region = dist_y_bd <= ext_regiony_scell;
                                const bool dist_z_bd_in_region = dist_z_bd <= ext_regionz_scell;

                                if (num_away == 1)
                                {
                                    if (body_idx_tar != body_idx_src)
                                    {
                                        WeightFlags entry_src(bx, by, bz, true, true, true);
                                        auto &entry_map = pair_list_aux[body_idx_tar][body_idx_src];
                                        WeightFlags *entry_pre = entry_map.find(entry_src);
                                        if (entry_pre)
                                        {
                                            entry_pre->bx |= (entry_tar.x_in ^ entry_pre->x_in);
                                            entry_pre->by |= (entry_tar.y_in ^ entry_pre->y_in);
                                            entry_pre->bz |= (entry_tar.z_in ^ entry_pre->z_in);
                                            entry_tar = *entry_pre;
                                        }
                                        entry_map.insert(entry_src, entry_tar);
                                    }
                                }
                                else if (num_away == 2)
                                {

                                    const bool sx_within = !(is_distx_largest && dist_x_bd_in_region);
                                    const bool sy_within = !(is_disty_largest && dist_y_bd_in_region);
                                    const bool sz_within = !(is_distz_largest && dist_z_bd_in_region);
                                    const bool is_valid_interaction =
                                        !((is_distx_largest && !dist_x_bd_in_region) || (is_disty_largest && !dist_y_bd_in_region) || (is_distz_largest && !dist_z_bd_in_region));

                                    if (is_valid_interaction && (!sx_within || !sy_within || !sz_within))
                                    {
                                        WeightFlags entry_src(bx, by, bz, sx_within, sy_within, sz_within);

                                        auto &entry_map = pair_list_aux[body_idx_tar][body_idx_src];
                                        WeightFlags *entry_pre = entry_map.find(entry_src);
                                        if (entry_pre)
                                        {
                                            entry_pre->bx |= (entry_tar.x_in ^ entry_pre->x_in);
                                            entry_pre->by |= (entry_tar.y_in ^ entry_pre->y_in);
                                            entry_pre->bz |= (entry_tar.z_in ^ entry_pre->z_in);
                                            entry_tar = *entry_pre;
                                        }
                                        entry_map.insert(entry_src, entry_tar);
                                    }
                                }
                            }
                            bidxt++;
                        }
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < bodies_all_.size(); i++)
    {
        for (const auto &[body_idx_src, pairListMap] : pair_list_aux[i])
        {
            for (const auto &[srcFlags, tarFlags] : pairListMap)
            {
                pair_list[i].emplace_back(body_idx_src, srcFlags.bx, srcFlags.by, srcFlags.bz, srcFlags.x_in, srcFlags.y_in, srcFlags.z_in,
                                          tarFlags.bx, tarFlags.by, tarFlags.bz, tarFlags.x_in, tarFlags.y_in, tarFlags.z_in);
            }
        }
    }
}

void gmx::fmm::FMMDirectInteractions::execute_direct_kernel(real *forces_and_potentials)
{
    size_t fp_idx = 0;
    for (size_t i = 0; i < bodies_all_.size(); i++)
    {
        const FBody &body_tar = bodies_all_[i];
        const real xt = body_tar.x[0];
        const real yt = body_tar.x[1];
        const real zt = body_tar.x[2];
        const RVec wtar_ws = w_per_atom[i];

        real pj_effective = 0.0;
        real fxj_effective = 0.0, fyj_effective = 0.0, fzj_effective = 0.0;

        for (auto &ent : pair_list[i])
        {

            gmx::fmm::FBody &asrc = bodies_all_[ent.body_idx_src];

            const BVec bxyz_src = {ent.bx_src, ent.by_src, ent.bz_src};
            const BVec bxyz_tar = {ent.bx_tar, ent.by_tar, ent.bz_tar};
            const BVec is_wihin_src = {ent.sx_within, ent.sy_within, ent.sz_within};
            const BVec is_within_tar = {ent.tx_within, ent.ty_within, ent.tz_within};

            const RVec wsrc_ws = w_per_atom[ent.body_idx_src];

            const real xs = asrc.x[0];
            const real ys = asrc.x[1];
            const real zs = asrc.x[2];

            const real dx = xt - xs;
            const real dy = yt - ys;
            const real dz = zt - zs;

            real wsrc_x = (bxyz_src[0] == 1) + (bxyz_src[0] != 1) * ((is_wihin_src[0] == 1) * wsrc_ws[0] + (is_wihin_src[0] != 1) * (1 - wsrc_ws[0]));
            real wsrc_y = (bxyz_src[1] == 1) + (bxyz_src[1] != 1) * ((is_wihin_src[1] == 1) * wsrc_ws[1] + (is_wihin_src[1] != 1) * (1 - wsrc_ws[1]));
            real wsrc_z = (bxyz_src[2] == 1) + (bxyz_src[2] != 1) * ((is_wihin_src[2] == 1) * wsrc_ws[2] + (is_wihin_src[2] != 1) * (1 - wsrc_ws[2]));
            real wsrc = wsrc_x * wsrc_y * wsrc_z;

            real wtar_x = (bxyz_tar[0] == 1) + (bxyz_tar[0] != 1) * ((is_within_tar[0] == 1) * wtar_ws[0] + (is_within_tar[0] != 1) * (1 - wtar_ws[0]));
            real wtar_y = (bxyz_tar[1] == 1) + (bxyz_tar[1] != 1) * ((is_within_tar[1] == 1) * wtar_ws[1] + (is_within_tar[1] != 1) * (1 - wtar_ws[1]));
            real wtar_z = (bxyz_tar[2] == 1) + (bxyz_tar[2] != 1) * ((is_within_tar[2] == 1) * wtar_ws[2] + (is_within_tar[2] != 1) * (1 - wtar_ws[2]));
            const real wtar = wtar_x * wtar_y * wtar_z;

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
        }
        // Apply accumulated forces and potential to target bodies
        forces_and_potentials[fp_idx++] = -fxj_effective;
        forces_and_potentials[fp_idx++] = -fyj_effective;
        forces_and_potentials[fp_idx++] = -fzj_effective;
        forces_and_potentials[fp_idx++] = pj_effective;
    }
}

void gmx::fmm::FMMDirectInteractions::recompute_weights() { compute_weights_(); }

void gmx::fmm::FMMDirectInteractions::rebuild_and_reprocess_tree() { fmm_direct_interactions_tree_.rebuild_and_reprocess_tree(); }
