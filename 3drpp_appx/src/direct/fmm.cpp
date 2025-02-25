#include "fmm.h"
#include <array>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <numeric>

gmx::fmm::FMMDirectInteractions::FMMDirectInteractions(const std::vector<RVec> coordinates, const std::vector<real> charges, const RVec box_center, const real box_radius,
                                                       const size_t max_depth, const real reg_alpha)
    : bodies_all_(coordinates, charges), fmm_weights_eval_(box_center, box_radius, reg_alpha), fmm_direct_interactions_tree_(bodies_all_, box_center, box_radius, max_depth)
{
    TIME_BEGIN(weights_time_new);
    compute_weights_();
    TIME_END(weights_time_new);
}

void gmx::fmm::FMMDirectInteractions::compute_weights_()
{
    const real reg_alpha = fmm_weights_eval_.getRegAlpha();

    pair_list_bits_.clear();

    FMMCells &fmm_cells = fmm_direct_interactions_tree_.get_cells();
    std::vector<FPIndices> boundary_bodies_idxs(fmm_cells.size());

    std::vector<bool> is_reg_body(bodies_all_.size(), false);
    std::vector<FPIndices> bodies_indices_ext(fmm_cells.size());
    std::vector<std::vector<WeightFlags>> w_flags(fmm_cells.size());

    bodies_cells.resize(bodies_all_.size(), -1);

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        for (const int &body_idx : cell.bodiesIndices)
        {

            FBody &body = bodies_all_[body_idx];
            bodies_cells[body_idx] = k;
            const RVec ws = fmm_weights_eval_.compute_weight_in_cell(body.x, cell.center, cell.radius, false);
            const real w = ws[0] * ws[1] * ws[2];
            body.w = ws;

            if (w < 1)
            {
                boundary_bodies_idxs[k].push_back(body_idx);
                is_reg_body[body_idx] = true;
            }
            bodies_indices_ext[k].push_back(body_idx);
            const bool is_x_fully_in = ws[0] == 1;
            const bool is_y_fully_in = ws[1] == 1;
            const bool is_z_fully_in = ws[2] == 1;
            w_flags[k].emplace_back(is_x_fully_in, is_y_fully_in, is_z_fully_in, 1, 1, 1);
        }
    }

    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        FMMCell &cell = fmm_cells[k];
        const real width_of_tcell = cell.radius * 2;

        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            int belongs_to_cells[7] = {-1, -1, -1, -1, -1, -1, -1};
            size_t btc_idx = 0;
            FBody &body_tar = bodies_all_[body_idx_tar];

            if (is_reg_body[body_idx_tar])
            {
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
                                const FMMCell &adj_cell = fmm_cells[adj_cell_idx];
                                const RVec dx = body_tar.x - adj_cell.center;

                                const real dist_x = fabs(dx[0]);
                                const real dist_y = fabs(dx[1]);
                                const real dist_z = fabs(dx[2]);

                                const bool is_dist_x_in_range = dist_x <= adj_cell.radius + reg_alpha;
                                const bool is_dist_y_in_range = dist_y <= adj_cell.radius + reg_alpha;
                                const bool is_dist_z_in_range = dist_z <= adj_cell.radius + reg_alpha;

                                if (is_dist_x_in_range && is_dist_y_in_range && is_dist_z_in_range)
                                {
                                    belongs_to_cells[btc_idx++] = adj_cell_idx;
                                }
                            }
                        }
                    }
                }
            }
            u_int32_t gid = get_group_id(cell.index, belongs_to_cells, btc_idx);
            body_tar.gid = gid;
        }
    }
    group_map.clear();
    std::vector<bool> should_compute_w(bodies_all_.size(), false);

    const u_int32_t num_groups = get_num_groups();
    std::vector<int> group_bodies;
    group_bodies.resize(num_groups, -1);
    for (const FBody &body : bodies_all_)
    {
        int group_id = body.gid;
        int pre_group_bidx = group_bodies[group_id];
        if (pre_group_bidx == -1)
        {
            group_bodies[group_id] = body.idx;
            should_compute_w[body.idx] = true;
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
                            const RVec ws = bodies_all_[body_idx].w;
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
                                bodies_indices_ext[k].push_back(body_idx);
                                w_flags[k].emplace_back(is_x_fully_in, is_y_fully_in, is_z_fully_in, dist_x <= cell.radius, dist_y <= cell.radius, dist_z <= cell.radius);
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<std::vector<WeightFlags>> w_flags_filt(fmm_cells.size());
    std::vector<FPIndices> bodies_idxs_to_comp(fmm_cells.size());
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        size_t ix = 0;
        for (const int &body_idx : bodies_indices_ext[k])
        {
            if (should_compute_w[body_idx])
            {
                bodies_idxs_to_comp[k].push_back(body_idx);
                w_flags_filt[k].emplace_back(w_flags[k][ix]);
            }
            ix++;
        }
    }

    w_flags.clear();
    bodies_indices_ext.clear();
    is_reg_body.clear();
    boundary_bodies_idxs.clear();
    w_flags = std::move(w_flags_filt);

    std::vector<std::map<int, FixedPairListMap, std::greater<int>>> pair_list_aux(num_groups);
    
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        const real region_of_tcell = cell.radius * 3;
        const real ext_region_tcell = region_of_tcell + reg_alpha;
        const real rdu_region_tcell = region_of_tcell - reg_alpha;
        const real width_of_tcell = cell.radius * 2;
        size_t bidxt = 0;
        for (const int &body_idx_tar : bodies_idxs_to_comp[k])
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];
            const WeightFlags entry_tar_btar(w_flags[k][bidxt].bx, w_flags[k][bidxt].by, w_flags[k][bidxt].bz, w_flags[k][bidxt].x_in, w_flags[k][bidxt].y_in,
                                             w_flags[k][bidxt].z_in);

            for (const int body_idx_src : cell.bodiesIndices)
            {
                const FBody &body_src = bodies_all_[body_idx_src];
                WeightFlags entry_tar = entry_tar_btar;
                if (body_tar.gid != body_src.gid && body_idx_tar < body_idx_src)
                {
                    WeightFlags entry_src(true, true, true, true, true, true);

                    auto &entry_map = pair_list_aux[body_tar.gid][body_idx_src];
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

                        const bool is_distx_largest = ((dist_x >= dist_y || fabs(dist_x - dist_y) <= tol) && (dist_x >= dist_z || fabs(dist_x - dist_z) <= tol));

                        const bool is_disty_largest = ((dist_y >= dist_x || fabs(dist_y - dist_x) <= tol) && (dist_y >= dist_z || fabs(dist_y - dist_z) <= tol));

                        const bool is_distz_largest = ((dist_z >= dist_x || fabs(dist_z - dist_x) <= tol) && (dist_z >= dist_y || fabs(dist_z - dist_y) <= tol));

                        short num_away = 1; // Default to "one away"

                        if (std::abs(dx) == 2 || std::abs(dy) == 2 || std::abs(dz) == 2)
                        {
                            num_away = 2; // If any direction is "two away"
                        }

                        size_t bidxt_in = 0;
                        for (const int &body_idx_tar : bodies_idxs_to_comp[k])
                        {
                            const FBody &body_tar = bodies_all_[body_idx_tar];
                            const WeightFlags entry_tar_btar(w_flags[k][bidxt_in].bx, w_flags[k][bidxt_in].by, w_flags[k][bidxt_in].bz, w_flags[k][bidxt_in].x_in,
                                                             w_flags[k][bidxt_in].y_in, w_flags[k][bidxt_in].z_in);
                            for (const int &body_idx_src : adj_cell.bodiesIndices)
                            {
                                const FBody &body_src = bodies_all_[body_idx_src];
                                if (body_src.gid != body_tar.gid && body_idx_tar < body_idx_src)
                                {
                                    WeightFlags entry_tar = entry_tar_btar;

                                    const RVec ws = body_src.w;
                                    const real dist_x_bd = fabs(body_src.x[0] - cell.center[0]);
                                    const real dist_y_bd = fabs(body_src.x[1] - cell.center[1]);
                                    const real dist_z_bd = fabs(body_src.x[2] - cell.center[2]);

                                    const bool bx = (ws[0] == 1) || (dist_x_bd <= rdu_region_tcell);
                                    const bool by = (ws[1] == 1) || (dist_y_bd <= rdu_region_tcell);
                                    const bool bz = (ws[2] == 1) || (dist_z_bd <= rdu_region_tcell);

                                    const bool dist_x_bd_in_region = dist_x_bd <= ext_region_tcell;
                                    const bool dist_y_bd_in_region = dist_y_bd <= ext_region_tcell;
                                    const bool dist_z_bd_in_region = dist_z_bd <= ext_region_tcell;

                                    if (num_away == 1)
                                    {

                                        WeightFlags entry_src(bx, by, bz, true, true, true);
                                        auto &entry_map = pair_list_aux[body_tar.gid][body_idx_src];
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
                                    else if (num_away == 2)
                                    {

                                        const bool sx_within = !(is_distx_largest && dist_x_bd_in_region);
                                        const bool sy_within = !(is_disty_largest && dist_y_bd_in_region);
                                        const bool sz_within = !(is_distz_largest && dist_z_bd_in_region);
                                        const bool is_valid_interaction = dist_x_bd_in_region && dist_y_bd_in_region && dist_z_bd_in_region;

                                        if (is_valid_interaction && (!sx_within || !sy_within || !sz_within))
                                        {
                                            WeightFlags entry_src(bx, by, bz, sx_within, sy_within, sz_within);

                                            auto &entry_map = pair_list_aux[body_tar.gid][body_idx_src];
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
                            }
                            bidxt_in++;
                        }
                    }
                }
            }
        }
    }
    pair_list_bits_.resize(num_groups);
    pair_list_bidx_srcs_.resize(num_groups);

    for (size_t i = 0; i < num_groups; i++)
    {
        for (const auto &[body_idx_src, pairListMap] : pair_list_aux[i])
        {
            for (const auto &[srcFlags, tarFlags] : pairListMap)
            {
                pair_list_bits_[i].emplace_back(srcFlags.bx, srcFlags.by, srcFlags.bz, srcFlags.x_in, srcFlags.y_in, srcFlags.z_in, tarFlags.bx, tarFlags.by, tarFlags.bz,
                                                tarFlags.x_in, tarFlags.y_in, tarFlags.z_in);

                pair_list_bidx_srcs_[i].push_back(body_idx_src);
            }
        }
    }
    pair_list_aux.clear();
}

void gmx::fmm::FMMDirectInteractions::compute_group_interactions_(const std::vector<FBody> &sbodies, const std::vector<size_t> &group_counts,
                                                                  const std::vector<size_t> &group_prefix_sum, real *forces_potentials)
{

    const u_int32_t num_groups = get_num_groups();
    for (size_t g_num = 0; g_num < num_groups; ++g_num)
    {
        size_t k = group_prefix_sum[g_num];

        const size_t num_particles_in_grp = group_counts[g_num];
        const size_t num_particles_in_grp_m1 = num_particles_in_grp - 1;

        std::vector<real> p_opp(num_particles_in_grp, 0);
        std::vector<real> fx_opp(num_particles_in_grp, 0);
        std::vector<real> fy_opp(num_particles_in_grp, 0);
        std::vector<real> fz_opp(num_particles_in_grp, 0);

        for (size_t kit = 0; kit < num_particles_in_grp_m1; ++kit)
        {
            const FBody &btar = sbodies[k + kit];
            const real qt = btar.q;

            real pji = 0.0;
            real fxji = 0.0, fyji = 0.0, fzji = 0.0;

            for (size_t ki = kit + 1; ki < num_particles_in_grp; ++ki)
            {
                const FBody &body_src = sbodies[k + ki];

                const real dx = btar.x[0] - body_src.x[0];
                const real dy = btar.x[1] - body_src.x[1];
                const real dz = btar.x[2] - body_src.x[2];

                const real r2 = dx * dx + dy * dy + dz * dz;
                const real invr = 1.0 / std::sqrt(r2);

                const real qsinvr = body_src.q * invr;
                const real qsinvr3 = qsinvr * invr * invr;

                pji += qsinvr;
                fxji += qsinvr3 * dx;
                fyji += qsinvr3 * dy;
                fzji += qsinvr3 * dz;

                const real qtinvr = qt * invr;
                const real qtinvr3 = qtinvr * invr * invr;

                const real fxij = qtinvr3 * -dx;
                const real fyij = qtinvr3 * -dy;
                const real fzij = qtinvr3 * -dz;

                p_opp[ki] += qtinvr;
                fx_opp[ki] += fxij;
                fy_opp[ki] += fyij;
                fz_opp[ki] += fzij;
            }

            const size_t btidx = btar.idx * 4;
            forces_potentials[btidx] -= fxji;
            forces_potentials[btidx + 1] -= fyji;
            forces_potentials[btidx + 2] -= fzji;
            forces_potentials[btidx + 3] += pji;
        }

        for (size_t ki = 0; ki < num_particles_in_grp; ++ki)
        {
            const size_t bsidx = sbodies[k + ki].idx * 4;
            forces_potentials[bsidx] -= fx_opp[ki];
            forces_potentials[bsidx + 1] -= fy_opp[ki];
            forces_potentials[bsidx + 2] -= fz_opp[ki];
            forces_potentials[bsidx + 3] += p_opp[ki];
        }
    }
}

u_int32_t gmx::fmm::FMMDirectInteractions::get_group_id(int ocell_idx, int a_cells_idxs[], size_t validSize)
{
    std::sort(a_cells_idxs, a_cells_idxs + validSize);

    std::ostringstream oss;
    oss << ocell_idx << ",";
    for (size_t i = 0; i < validSize; i++)
    {
        oss << ",";
        oss << a_cells_idxs[i];
    }
    std::string key = oss.str();

    if (group_map.find(key) != group_map.end())
    {
        return group_map[key];
    }

    group_map[key] = next_group_id++;
    return group_map[key];
}

void gmx::fmm::FMMDirectInteractions::execute_direct_kernel(real *forces_and_potentials)
{
    for (size_t i = 0, btidx = 0; i < bodies_all_.size(); i++, btidx += 4)
    {
        const FBody &body_tar = bodies_all_[i];

        const real xt = body_tar.x[0];
        const real yt = body_tar.x[1];
        const real zt = body_tar.x[2];

        const RVec wtar_ws = body_tar.w;
        const int gidt = body_tar.gid;
        const real qt = body_tar.q;

        real pji = 0.0;
        real fxji = 0.0, fyji = 0.0, fzji = 0.0;

        const auto &bidx_srcs = pair_list_bidx_srcs_[gidt];
        const auto &bits = pair_list_bits_[gidt];

        for (size_t ix = 0; ix < bidx_srcs.size() && body_tar.idx < bidx_srcs[ix]; ++ix)
        {
            const int body_src_idx = bidx_srcs[ix];
            FBody &body_src = bodies_all_[body_src_idx];
            const PairListEntry &ent = bits[ix];

            const BVec bxyz_src = {static_cast<bool>((ent.packed_flags >> PairListEntry::shiftBXS) & 1), static_cast<bool>((ent.packed_flags >> PairListEntry::shiftBYS) & 1),
                                   static_cast<bool>((ent.packed_flags >> PairListEntry::shiftBZS) & 1)};

            const BVec bxyz_tar = {static_cast<bool>((ent.packed_flags >> PairListEntry::shiftBXT) & 1), static_cast<bool>((ent.packed_flags >> PairListEntry::shiftBYT) & 1),
                                   static_cast<bool>((ent.packed_flags >> PairListEntry::shiftBZT) & 1)};

            const BVec is_within_src = {static_cast<bool>((ent.packed_flags >> PairListEntry::shiftSXW) & 1), static_cast<bool>((ent.packed_flags >> PairListEntry::shiftSYW) & 1),
                                        static_cast<bool>((ent.packed_flags >> PairListEntry::shiftSZW) & 1)};

            const BVec is_within_tar = {static_cast<bool>((ent.packed_flags >> PairListEntry::shiftTXW) & 1), static_cast<bool>((ent.packed_flags >> PairListEntry::shiftTYW) & 1),
                                        static_cast<bool>((ent.packed_flags >> PairListEntry::shiftTZW) & 1)};

            const RVec wsrc_ws = body_src.w;

            const real xs = body_src.x[0];
            const real ys = body_src.x[1];
            const real zs = body_src.x[2];

            const real dx = xt - xs;
            const real dy = yt - ys;
            const real dz = zt - zs;

            const real wsrc_x = bxyz_src[0] == 1 ? 1 : (is_within_src[0] ? wsrc_ws[0] : 1 - wsrc_ws[0]);
            const real wsrc_y = bxyz_src[1] == 1 ? 1 : (is_within_src[1] ? wsrc_ws[1] : 1 - wsrc_ws[1]);
            const real wsrc_z = bxyz_src[2] == 1 ? 1 : (is_within_src[2] ? wsrc_ws[2] : 1 - wsrc_ws[2]);
            const real wsrc = wsrc_x * wsrc_y * wsrc_z;

            const real wtar_x = bxyz_tar[0] == 1 ? 1 : (is_within_tar[0] ? wtar_ws[0] : 1 - wtar_ws[0]);
            const real wtar_y = bxyz_tar[1] == 1 ? 1 : (is_within_tar[1] ? wtar_ws[1] : 1 - wtar_ws[1]);
            const real wtar_z = bxyz_tar[2] == 1 ? 1 : (is_within_tar[2] ? wtar_ws[2] : 1 - wtar_ws[2]);

            const real wtar = wtar_x * wtar_y * wtar_z;

            const real r2 = dx * dx + dy * dy + dz * dz;
            const real invr = 1.0 / std::sqrt(r2);

            const real qs = body_src.q * wsrc;
            const real qsinvr = qs * invr;
            const real qsinvr3 = qsinvr * invr * invr;

            pji += qsinvr * wtar;
            fxji += qsinvr3 * dx * wtar;
            fyji += qsinvr3 * dy * wtar;
            fzji += qsinvr3 * dz * wtar;

            const real qtinvr = qt * invr * wtar;
            const real qtinvr3 = qtinvr * invr * invr;

            const size_t bsidx = body_src.idx * 4;
            forces_and_potentials[bsidx] -= qtinvr3 * -dx * wsrc;
            forces_and_potentials[bsidx + 1] -= qtinvr3 * -dy * wsrc;
            forces_and_potentials[bsidx + 2] -= qtinvr3 * -dz * wsrc;
            forces_and_potentials[bsidx + 3] += qtinvr * wsrc;
        }

        forces_and_potentials[btidx] -= fxji;
        forces_and_potentials[btidx + 1] -= fyji;
        forces_and_potentials[btidx + 2] -= fzji;
        forces_and_potentials[btidx + 3] += pji;
    }

    std::vector<FBody> sbodies = bodies_all_;

    std::sort(sbodies.begin(), sbodies.end(), [](const FBody &a, const FBody &b) { return a.gid < b.gid; });

    const uint32_t num_groups = get_num_groups();

    std::vector<size_t> group_counts(num_groups, 0);
    std::vector<size_t> group_prefix_sum(num_groups, 0);

    for (size_t i = 0, g = 0; i < sbodies.size(); ++i)
    {
        if (i > 0 && sbodies[i].gid != sbodies[i - 1].gid)
        {

            ++g;
        }

        group_counts[g]++;
    }

    std::exclusive_scan(group_counts.begin(), group_counts.end(), group_prefix_sum.begin(), 0);

    compute_group_interactions_(sbodies, group_counts, group_prefix_sum, forces_and_potentials);
}

void gmx::fmm::FMMDirectInteractions::recompute_weights() { compute_weights_(); }

void gmx::fmm::FMMDirectInteractions::rebuild_and_reprocess_tree() { fmm_direct_interactions_tree_.rebuild_and_reprocess_tree(); }

u_int32_t gmx::fmm::FMMDirectInteractions::get_num_groups() { return next_group_id; }
