#include "fmm.h"
#include <array>
#include <fstream>
#include <iomanip>

gmx::fmm::FMMDirectInteractions::FMMDirectInteractions(const std::vector<RVec> coordinates, const std::vector<real> charges, const RVec box_center, const real box_radius, const size_t max_depth,
                                                       const real reg_alpha)
    : bodies_all_(coordinates, charges), fmm_weights_eval_(box_center, box_radius, reg_alpha), fmm_direct_interactions_tree_(bodies_all_, box_center, box_radius, max_depth)
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

    w_per_atom.clear();
    w_per_atom.resize(bodies_all_.size());

    FMMCells &fmm_cells = fmm_direct_interactions_tree_.get_cells();

    std::vector<FPIndices> boundary_bodies_idxs(fmm_cells.size());
    std::vector<bool> is_reg_body(bodies_all_.size(), false);
    for (size_t k = 0; k < fmm_cells.size(); k++)
    {
        const FMMCell &cell = fmm_cells[k];
        for (const int &body_idx : cell.bodiesIndices)
        {
            const FBody &body = bodies_all_[body_idx];
            const real w = fmm_weights_eval_.compute_weight_within_cell(body.x, cell.center, cell.radius, false);

            if (w < 1)
            {
                boundary_bodies_idxs[k].push_back(body_idx);
                is_reg_body[body_idx] = true;
            }
        }
    }

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
        for (const int &body_idx_tar : cell.bodiesIndices)
        {
            const FBody &body_tar = bodies_all_[body_idx_tar];

            for (const int body_idx_src : cell.bodiesIndices)
            {
                if (body_idx_tar != body_idx_src)
                {
                    const FBody &body_src = bodies_all_[body_idx_src];

                    const RVec w3t = w_per_atom[body_idx_tar];
                    const RVec w3s = w_per_atom[body_idx_src];
                    real wtar_in_cell = w3t[0] * w3t[1] * w3t[2];
                    real wsrc_in_cell = w3s[0] * w3s[1] * w3s[2];

                    if (wtar_in_cell == 1 || wsrc_in_cell == 1)
                    {
                        PairListEntry entry(body_idx_src);
                        entry.set_tar_flags(1, 1, 1);
                        entry.set_trw_flags(1, 1, 1);
                        entry.set_src_flags(1, 1, 1);
                        entry.set_scw_flags(1, 1, 1);
                        pair_list[body_idx_tar].push_back(entry);
                    }
                    else
                    {

                        bool x_dir_flag = true;
                        bool y_dir_flag = true;
                        bool z_dir_flag = true;

                        real r1x = body_tar.x[0] - cell.center[0];
                        real r2x = body_src.x[0] - cell.center[0];

                        // both particles are in opposite direction of x within
                        // the cell
                        if (r1x * r2x < 0)
                        {
                            x_dir_flag = false;
                        }

                        // both particles are in opposite direction of y within
                        // the cell
                        real r1y = body_tar.x[1] - cell.center[1];
                        real r2y = body_src.x[1] - cell.center[1];
                        if (r1y * r2y < 0)
                        {
                            y_dir_flag = false;
                        }

                        // both particles are in opposite direction of z within
                        // the cell
                        real r1z = body_tar.x[2] - cell.center[2];
                        real r2z = body_src.x[2] - cell.center[2];
                        if (r1z * r2z < 0)
                        {
                            z_dir_flag = false;
                        }

                        // code below needs double checking

                        bool bxts = w3t[0] != 1 && w3s[0] != 1 ? x_dir_flag : 1;
                        bool byts = w3t[1] != 1 && w3s[1] != 1 ? y_dir_flag : 1;
                        bool bzts = w3t[2] != 1 && w3s[2] != 1 ? z_dir_flag : 1;
                        pair_list[body_idx_tar].push_back(body_idx_src);

                        PairListEntry entry(body_idx_src);
                        entry.set_tar_flags(bxts, byts, bzts);
                        entry.set_trw_flags(1, 1, 1);
                        entry.set_src_flags(1, 1, 1);
                        entry.set_scw_flags(1, 1, 1);
                        pair_list[body_idx_tar].push_back(entry);

                        // indirect weights outside the boundary of the cell
                        if (bxts == false || byts == false || bzts == false)
                        {
                            PairListEntry entry(body_idx_src);
                            entry.set_tar_flags(bxts, byts, bzts);
                            entry.set_trw_flags(1, 1, 1);
                            entry.set_src_flags(bxts, byts, bzts);
                            entry.set_scw_flags(0, 0, 0);
                            pair_list[body_idx_tar].push_back(entry);
                        }
                    }
                }
            }
        }

        const real width_of_tcell = cell.radius * 2;

        for (int dz = -3; dz <= 3; dz++)
        {
            for (int dy = -3; dy <= 3; dy++)
            {
                for (int dx = -3; dx <= 3; dx++)
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

                        const real dist_1cell = cell.radius * 2;
                        const real dist_2cells = cell.radius * 4;
                        const real dist_3cells = cell.radius * 6;

                        short num_away = 1; // Default to "one away"
                        if (std::abs(dx) == 3 || std::abs(dy) == 3 || std::abs(dz) == 3)
                        {
                            num_away = 3; // If any direction is "three away"
                        }
                        else if (std::abs(dx) == 2 || std::abs(dy) == 2 || std::abs(dz) == 2)
                        {
                            num_away = 2; // If any direction is "two away"
                        }

                        const FPIndices &tar_bodies_idxs = (num_away == 3) ? boundary_bodies_idxs[k] : cell.bodiesIndices;
                        const FPIndices &src_bodies_idxs = (num_away == 3) ? boundary_bodies_idxs[adj_cell_idx] : adj_cell.bodiesIndices;

                        for (const int &body_idx_tar : tar_bodies_idxs)
                        {
                            const FBody &body_tar = bodies_all_[body_idx_tar];
                            const RVec wst = w_per_atom[body_idx_tar];

                            const bool in_regx_tar = wst[0] != 1;
                            const bool in_regy_tar = wst[1] != 1;
                            const bool in_regz_tar = wst[2] != 1;

                            const real rtx_tc = body_tar.x[0] - cell.center[0];
                            const real rty_tc = body_tar.x[1] - cell.center[1];
                            const real rtz_tc = body_tar.x[2] - cell.center[2];

                            const real rtx_sc = body_tar.x[0] - adj_cell.center[0];
                            const real rty_sc = body_tar.x[1] - adj_cell.center[1];
                            const real rtz_sc = body_tar.x[2] - adj_cell.center[2];

                            // only valid if dist_x != 0, dist_y != 0 and in
                            // case of dist_z != 0
                            const bool is_xtar_bw_cells = rtx_tc * rtx_sc < 0;
                            const bool is_ytar_bw_cells = rty_tc * rty_sc < 0;
                            const bool is_ztar_bw_cells = rtz_tc * rtz_sc < 0;

                            // one and half cell
                            for (const int &body_idx_src : adj_cell.bodiesIndices)
                            {

                                const RVec ws = w_per_atom[body_idx_src];
                                const FBody &body_src = bodies_all_[body_idx_src];

                                const bool in_regx_src = ws[0] != 1;
                                const bool in_regy_src = ws[1] != 1;
                                const bool in_regz_src = ws[2] != 1;

                                const real wtar_in_cell = wst[0] * wst[1] * wst[2];
                                const real wsrc_in_cell = ws[0] * ws[1] * ws[2];

                                const real rsx_tc = body_src.x[0] - cell.center[0];
                                const real rsx_sc = body_src.x[0] - adj_cell.center[0];

                                const real rsy_tc = body_src.x[1] - cell.center[1];
                                const real rsy_sc = body_src.x[1] - adj_cell.center[1];

                                const real rsz_tc = body_src.x[2] - cell.center[2];
                                const real rsz_sc = body_src.x[2] - adj_cell.center[2];

                                // only valid if dist_x != 0, dist_y != 0 and in
                                // case of dist_z != 0
                                const bool is_xsrc_bw_cells = rsx_tc * rsx_sc < 0;
                                const bool is_ysrc_bw_cells = rsy_tc * rsy_sc < 0;
                                const bool is_zsrc_bw_cells = rsz_tc * rsz_sc < 0;

                                // only valid if dist_x == 0, dist_y == 0 and in case of dist_z == 0
                                const bool x_same_dir_ts = !(rtx_tc * rsx_tc < 0);
                                const bool y_same_dir_ts = !(rty_tc * rsy_tc < 0);
                                const bool z_same_dir_ts = !(rtz_tc * rsz_tc < 0);

                                // same as above, avoiding redundancy
                                // const bool x_same_dir_src = !(rsx_sc * rtx_sc < 0);
                                // const bool y_same_dir_src = !(rsy_sc * rty_sc < 0);
                                // const bool z_same_dir_src = !(rsz_sc * rtz_sc < 0);

                                const bool in_rel_regx_tar = in_regx_tar && ((dist_x == 0 && x_same_dir_ts) || (dist_x != 0 && is_xtar_bw_cells));
                                const bool in_rel_regy_tar = in_regy_tar && ((dist_y == 0 && y_same_dir_ts) || (dist_y != 0 && is_ytar_bw_cells));
                                const bool in_rel_regz_tar = in_regz_tar && ((dist_z == 0 && z_same_dir_ts) || (dist_z != 0 && is_ztar_bw_cells));

                                const bool in_rel_regx_src = in_regx_src && ((dist_x == 0 && x_same_dir_ts) || (dist_x != 0 && is_xsrc_bw_cells));
                                const bool in_rel_regy_src = in_regy_src && ((dist_y == 0 && y_same_dir_ts) || (dist_y != 0 && is_ysrc_bw_cells));
                                const bool in_rel_regz_src = in_regz_src && ((dist_z == 0 && z_same_dir_ts) || (dist_z != 0 && is_zsrc_bw_cells));

                                const bool bxt1 = !in_regx_tar;
                                const bool byt1 = !in_regy_tar;
                                const bool bzt1 = !in_regz_tar;

                                const bool bxs1 = !in_regx_src;
                                const bool bys1 = !in_regy_src;
                                const bool bzs1 = !in_regz_src;

                                const bool is_xt_fr_tscell = bxt1 || in_rel_regx_tar;
                                const bool is_yt_fr_tscell = byt1 || in_rel_regy_tar;
                                const bool is_zt_fr_tscell = bzt1 || in_rel_regz_tar;

                                const bool is_xs_fr_scell = bxs1 || in_rel_regx_src;
                                const bool is_ys_fr_scell = bys1 || in_rel_regy_src;
                                const bool is_zs_fr_scell = bzs1 || in_rel_regz_src;

                                const bool is_x_fr_cells = is_xt_fr_tscell && is_xs_fr_scell;
                                const bool is_y_fr_cells = is_yt_fr_tscell && is_ys_fr_scell;
                                const bool is_z_fr_cells = is_zt_fr_tscell && is_zs_fr_scell;

                                if (num_away == 1)
                                {
                                    if ((wtar_in_cell == 1 && wsrc_in_cell == 1) || (is_x_fr_cells && is_y_fr_cells && is_z_fr_cells))
                                    {

                                        PairListEntry entry(body_idx_src);
                                        entry.set_src_flags(1, 1, 1);
                                        entry.set_tar_flags(1, 1, 1);
                                        entry.set_scw_flags(1, 1, 1);
                                        entry.set_trw_flags(1, 1, 1);
                                        pair_list[body_idx_tar].push_back(entry);
                                    }
                                    else
                                    {

                                        // need to be fully sure of
                                        // (!in_regx_src && dist_x == 0) etc.,
                                        const bool bxt_c1 = (!in_regx_src && dist_x == 0) || is_xt_fr_tscell;
                                        const bool byt_c1 = (!in_regy_src && dist_y == 0) || is_yt_fr_tscell;
                                        const bool bzt_c1 = (!in_regz_src && dist_z == 0) || is_zt_fr_tscell;

                                        const bool bxs_c1 = !in_regx_src ? 1 : ((dist_x == 0) || (dist_x != 0 && is_xsrc_bw_cells));
                                        const bool bys_c1 = !in_regy_src ? 1 : ((dist_y == 0) || (dist_y != 0 && is_ysrc_bw_cells));
                                        const bool bzs_c1 = !in_regz_src ? 1 : ((dist_z == 0) || (dist_z != 0 && is_zsrc_bw_cells));

                                        PairListEntry entry(body_idx_src);
                                        entry.set_tar_flags(bxt_c1, byt_c1, bzt_c1);
                                        entry.set_trw_flags(1, 1, 1);
                                        entry.set_src_flags(bxs_c1, bys_c1, bzs_c1);
                                        entry.set_scw_flags(1, 1, 1);
                                        pair_list[body_idx_tar].push_back(entry);

                                        // these flags represent to take full weight if particle is not in regualization region
                                        // in a specific dimension
                                        const bool bxs_c1b = !in_regx_src;
                                        const bool bys_c1b = !in_regy_src;
                                        const bool bzs_c1b = !in_regz_src;

                                        const bool bxt_c1b = !in_regx_tar;
                                        const bool byt_c1b = !in_regy_tar;
                                        const bool bzt_c1b = !in_regz_tar;

                                        // otherwise we need to decide based on the relative position of the particles
                                        // now, a weight outside is only valid for further interaction if source and target particles are in regularized region
                                        // from a specific dimension, and they are pointing in opposite direction
                                        const bool sif_x = bxs_c1b || !(in_regx_tar && in_rel_regx_tar && !in_rel_regx_src);
                                        const bool sif_y = bys_c1b || !(in_regy_tar && in_rel_regy_tar && !in_rel_regy_src);
                                        const bool sif_z = bzs_c1b || !(in_regz_tar && in_rel_regz_tar && !in_rel_regz_src);

                                        const bool is_yaxis1_relevant = dist_y == dist_1cell && (in_rel_regy_tar || in_rel_regy_src);
                                        const bool is_zaxis1_relevant = dist_z == dist_1cell && (in_rel_regz_tar || in_rel_regz_src);
                                        const bool is_xaxis1_relevant = dist_x == dist_1cell && (in_rel_regx_tar || in_rel_regx_src);

                                        if (sif_x == false || sif_y == false || sif_z == false)
                                        {
                                            const bool tif_x = bxs_c1b || !in_rel_regx_tar;
                                            const bool tif_y = bys_c1b || !in_rel_regy_tar;
                                            const bool tif_z = bzs_c1b || !in_rel_regz_tar;

                                            const bool bxs_c1bb = (bxs_c1b || (!sif_y && (dist_x == 0 || is_xaxis1_relevant)) || (!sif_z && (dist_x == 0 || is_xaxis1_relevant)));

                                            const bool bys_c1bb = (bys_c1b || (!sif_x && (dist_y == 0 || is_yaxis1_relevant)) || (!sif_z && (dist_y == 0 || is_yaxis1_relevant)));

                                            const bool bzs_c1bb = (bzs_c1b || (!sif_x && (dist_z == 0 || is_zaxis1_relevant)) || (!sif_y && (dist_z == 0 || is_zaxis1_relevant)));

                                            const bool bxt_c1bb = (bxt_c1b || (!tif_y && (dist_x == 0 || is_xaxis1_relevant)) || (!tif_z && (dist_x == 0 || is_xaxis1_relevant)));

                                            const bool byt_c1bb = (byt_c1b || (!tif_x && (dist_y == 0 || is_yaxis1_relevant)) || (!tif_z && (dist_y == 0 || is_yaxis1_relevant)));

                                            const bool bzt_c1bb = (bzt_c1b || (!tif_x && (dist_z == 0 || is_zaxis1_relevant)) || (!tif_y && (dist_z == 0 || is_zaxis1_relevant)));

                                            PairListEntry entry(body_idx_src);
                                            entry.set_tar_flags(bxt_c1bb, byt_c1bb, bzt_c1bb);
                                            entry.set_trw_flags(tif_x, tif_y, tif_z);
                                            entry.set_src_flags(bxs_c1bb, bys_c1bb, bzs_c1bb);
                                            entry.set_scw_flags(sif_x, sif_y, sif_z);
                                            pair_list[body_idx_tar].push_back(entry);
                                        }

                                        // need more testing
                                        const bool sif_x0 = bxs_c1b || !(dist_x == 0 && in_regx_tar && !in_rel_regx_tar && !in_rel_regx_src);
                                        const bool sif_y0 = bys_c1b || !(dist_y == 0 && in_regy_tar && !in_rel_regy_tar && !in_rel_regy_src);
                                        const bool sif_z0 = bzs_c1b || !(dist_z == 0 && in_regz_tar && !in_rel_regz_tar && !in_rel_regz_src);

                                        if (sif_x0 == false || sif_y0 == false || sif_z0 == false)
                                        {
                                            const bool tif_x0 = !(dist_x == 0 && in_regx_tar && in_regx_src && !x_same_dir_ts);
                                            const bool tif_y0 = !(dist_y == 0 && in_regy_tar && in_regy_src && !y_same_dir_ts);
                                            const bool tif_z0 = !(dist_z == 0 && in_regz_tar && in_regz_src && !z_same_dir_ts);

                                            const bool bxs_c1bb = (bxs_c1b || (!sif_y0 && (dist_x == 0 || is_xaxis1_relevant)) || (!sif_z0 && (dist_x == 0 || is_xaxis1_relevant)));
                                            const bool bys_c1bb = (bys_c1b || (!sif_x0 && (dist_y == 0 || is_yaxis1_relevant)) || (!sif_z0 && (dist_y == 0 || is_yaxis1_relevant)));
                                            const bool bzs_c1bb = (bzs_c1b || (!sif_x0 && (dist_z == 0 || is_zaxis1_relevant)) || (!sif_y0 && (dist_z == 0 || is_zaxis1_relevant)));
                                            PairListEntry entry(body_idx_src);
                                            entry.set_tar_flags(bxt_c1, byt_c1, bzt_c1);
                                            entry.set_trw_flags(tif_x0, tif_y0, tif_z0);
                                            entry.set_src_flags(bxs_c1bb, bys_c1bb, bzs_c1bb);
                                            entry.set_scw_flags(sif_x, sif_y, sif_z);
                                            pair_list[body_idx_tar].push_back(entry);
                                        }

                                        // a weight outside is only valid for further interaction if source and target particles are in regularized region
                                        // from a specific dimension, and they are pointing in opposite direction
                                        const bool tif_x = bxt_c1b || !(in_regx_src && in_rel_regx_src && !in_rel_regx_tar);
                                        const bool tif_y = byt_c1b || !(in_regy_src && in_rel_regy_src && !in_rel_regy_tar);
                                        const bool tif_z = bzt_c1b || !(in_regz_src && in_rel_regz_src && !in_rel_regz_tar);

                                        if (tif_x == false || tif_y == false || tif_z == false)
                                        {
                                            const bool sif_x = bxt_c1b ? 1 : !in_rel_regx_src;
                                            const bool sif_y = byt_c1b ? 1 : !in_rel_regy_src;
                                            const bool sif_z = bzt_c1b ? 1 : !in_rel_regz_src;

                                            const bool bxs_c1bb = (bxs_c1b || (!sif_y && (dist_x == 0 || is_xaxis1_relevant)) || (!sif_z && (dist_x == 0 || is_xaxis1_relevant)));

                                            const bool bys_c1bb = (bys_c1b || (!sif_x && (dist_y == 0 || is_yaxis1_relevant)) || (!sif_z && (dist_y == 0 || is_yaxis1_relevant)));

                                            const bool bzs_c1bb = (bzs_c1b || (!sif_x && (dist_z == 0 || is_zaxis1_relevant)) || (!sif_y && (dist_z == 0 || is_zaxis1_relevant)));

                                            const bool bxt_c1bb = (bxt_c1b || (!tif_y && (dist_x == 0 || is_xaxis1_relevant)) || (!tif_z && (dist_x == 0 || is_xaxis1_relevant)));

                                            const bool byt_c1bb = (byt_c1b || (!tif_x && (dist_y == 0 || is_yaxis1_relevant)) || (!tif_z && (dist_y == 0 || is_yaxis1_relevant)));

                                            const bool bzt_c1bb = (bzt_c1b || (!tif_x && (dist_z == 0 || is_zaxis1_relevant)) || (!tif_y && (dist_z == 0 || is_zaxis1_relevant)));

                                            PairListEntry entry(body_idx_src);
                                            entry.set_tar_flags(bxt_c1bb, byt_c1bb, bzt_c1bb);
                                            entry.set_trw_flags(tif_x, tif_y, tif_z);
                                            entry.set_src_flags(bxs_c1bb, bys_c1bb, bzs_c1bb);
                                            entry.set_scw_flags(sif_x, sif_y, sif_z);
                                            pair_list[body_idx_tar].push_back(entry);
                                        }
                                    }
                                }
                                else if (num_away == 2) // num_away == 2
                                {
                                    const real interaction_region_tcell = cell.radius * 3;     // one and half cell
                                    const real interaction_region_scell = adj_cell.radius * 3; // one and half cell

                                    const bool bxt = !in_regx_tar || (fabs(body_tar.x[0] - adj_cell.center[0]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scell);
                                    const bool byt = !in_regy_tar || (fabs(body_tar.x[1] - adj_cell.center[1]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scell);
                                    const bool bzt = !in_regz_tar || (fabs(body_tar.x[2] - adj_cell.center[2]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scell);

                                    const bool bx = !in_regx_src || (fabs(body_src.x[0] - cell.center[0]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_tcell);
                                    const bool by = !in_regy_src || (fabs(body_src.x[1] - cell.center[1]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_tcell);
                                    const bool bz = !in_regz_src || (fabs(body_src.x[2] - cell.center[2]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_tcell);

                                    const bool sif_x = !(dist_x == dist_2cells && in_rel_regx_src);
                                    const bool sif_y = !(dist_y == dist_2cells && in_rel_regy_src);
                                    const bool sif_z = !(dist_z == dist_2cells && in_rel_regz_src);

                                    const bool tif_x = !(dist_x == dist_2cells && in_rel_regx_tar);
                                    const bool tif_y = !(dist_y == dist_2cells && in_rel_regy_tar);
                                    const bool tif_z = !(dist_z == dist_2cells && in_rel_regz_tar);

                                    // Check if the source particle is within
                                    // the x-range of the target cell
                                    const bool srcp_in_tarc_x = fabs(body_src.x[0] - cell.center[0]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the source particle is within
                                    // the y-range of the target cell
                                    const bool srcp_in_tarc_y = fabs(body_src.x[1] - cell.center[1]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the source particle is within
                                    // the z-range of the target cell
                                    const bool srcp_in_tarc_z = fabs(body_src.x[2] - cell.center[2]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the target particle is within
                                    // the x-range of the source cell
                                    const bool tarp_in_srcc_x = fabs(body_tar.x[0] - adj_cell.center[0]) <= interaction_region_scell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the target particle is within
                                    // the y-range of the source cell
                                    const bool tarp_in_srcc_y = fabs(body_tar.x[1] - adj_cell.center[1]) <= interaction_region_scell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the target particle is within
                                    // the z-range of the source cell
                                    const bool tarp_in_srcc_z = fabs(body_tar.x[2] - adj_cell.center[2]) <= interaction_region_scell + fmm_weights_eval_.getRegAlpha();

                                    // Determine if the source particle is
                                    // entirely within the range of the target
                                    // cell
                                    const bool srcp_in_tarc = srcp_in_tarc_x && srcp_in_tarc_y && srcp_in_tarc_z;

                                    // Determine if the target particle is
                                    // entirely within the range of the source
                                    // cell
                                    const bool tarp_in_srcc = tarp_in_srcc_x && tarp_in_srcc_y && tarp_in_srcc_z;

                                    if (!is_reg_body[body_idx_tar])
                                    {
                                        if (is_reg_body[body_idx_src])
                                        {
                                            if (srcp_in_tarc && (sif_x == false || sif_y == false || sif_z == false))
                                            {
                                                PairListEntry entry(body_idx_src);
                                                entry.set_tar_flags(bxt, byt, bzt);
                                                entry.set_trw_flags(1, 1, 1);
                                                entry.set_src_flags(bx, by, bz);
                                                entry.set_scw_flags(sif_x, sif_y, sif_z);
                                                pair_list[body_idx_tar].push_back(entry);
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (!is_reg_body[body_idx_src])
                                        {
                                            if (tarp_in_srcc && (tif_x == false || tif_y == false || tif_z == false))
                                            {
                                                PairListEntry entry(body_idx_src);
                                                entry.set_tar_flags(bxt, byt, bzt);
                                                entry.set_trw_flags(tif_x, tif_y, tif_z);
                                                entry.set_src_flags(bx, by, bz);
                                                entry.set_scw_flags(1, 1, 1);
                                                pair_list[body_idx_tar].push_back(entry);
                                            }
                                        }
                                        else if (!srcp_in_tarc && !tarp_in_srcc)
                                        {
                                            if (tif_x == false || tif_y == false || tif_z == false)
                                            {

                                                bool sx_c = 1;
                                                bool sy_c = 1;
                                                bool sz_c = 1;

                                                const bool is_zaxis_related = dist_z < dist_2cells || (in_rel_regz_src || in_rel_regz_tar);
                                                const bool is_yaxis_related = dist_y < dist_2cells || (in_rel_regy_src || in_rel_regy_tar);
                                                const bool is_xaxis_related = dist_x < dist_2cells || (in_rel_regx_src || in_rel_regx_tar);

                                                bool is_x_interacting_dim = false, is_y_interacting_dim = false, is_z_interacting_dim = false;
                                                bool multi_inf = false;

                                                // Interaction logic for each dimension
                                                if (srcp_in_tarc_x && !sif_x && ((!tif_y && is_zaxis_related) || (!tif_z && is_yaxis_related)))
                                                {
                                                    sx_c = 0;
                                                    is_x_interacting_dim = true;
                                                    is_y_interacting_dim = !tif_y && is_zaxis_related;
                                                    is_z_interacting_dim = !tif_z && is_yaxis_related;
                                                }

                                                if (srcp_in_tarc_y && !sif_y && ((!tif_x && is_zaxis_related) || (!tif_z && is_xaxis_related)))
                                                {
                                                    sy_c = 0;
                                                    is_y_interacting_dim = true;
                                                    is_x_interacting_dim = !tif_x && is_zaxis_related;
                                                    is_z_interacting_dim = !tif_z && is_xaxis_related;
                                                }

                                                if (srcp_in_tarc_z && !sif_z && ((!tif_x && is_yaxis_related) || (!tif_y && is_xaxis_related)))
                                                {
                                                    sz_c = 0;
                                                    is_z_interacting_dim = true;
                                                    is_x_interacting_dim = !tif_x && is_yaxis_related;
                                                    is_y_interacting_dim = !tif_y && is_xaxis_related;
                                                }

                                                // Multi-influence logic
                                                multi_inf = ((!is_x_interacting_dim && dist_x == 0 && in_regx_tar && in_regx_src && !in_rel_regx_tar && !in_rel_regx_src) ||
                                                             (!is_y_interacting_dim && dist_y == 0 && in_regy_tar && in_regy_src && !in_rel_regy_tar && !in_rel_regy_src) ||
                                                             (!is_z_interacting_dim && dist_z == 0 && in_regz_tar && in_regz_src && !in_rel_regz_tar && !in_rel_regz_src));

                                                // combining more than two could be problematic
                                                if (sx_c == false || sy_c == false || sz_c == false)
                                                {
                                                    if (!multi_inf)
                                                    {
                                                        PairListEntry entry(body_idx_src);
                                                        entry.set_tar_flags(bxt, byt, bzt);
                                                        entry.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry.set_src_flags(bx, by, bz);
                                                        entry.set_scw_flags(sx_c, sy_c, sz_c);
                                                        pair_list[body_idx_tar].push_back(entry);
                                                    }
                                                    else
                                                    {
                                                        if (!is_x_interacting_dim)
                                                        {
                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, byt, bzt);
                                                            entry1.set_trw_flags(1, tif_y, tif_z);
                                                            entry1.set_src_flags(bx, by, bz);
                                                            entry1.set_scw_flags(sx_c, sy_c, sz_c);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, byt, bzt);
                                                            entry2.set_trw_flags(0, tif_y, tif_z);
                                                            entry2.set_src_flags(0, by, bz);
                                                            entry2.set_scw_flags(1, sy_c, sz_c);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                        else if (!is_y_interacting_dim)
                                                        {
                                                            PairListEntry entry(body_idx_src);
                                                            entry.set_tar_flags(bxt, 0, bzt);
                                                            entry.set_trw_flags(tif_x, 1, tif_z);
                                                            entry.set_src_flags(bx, by, bz);
                                                            entry.set_scw_flags(sx_c, sy_c, sz_c);

                                                            pair_list[body_idx_tar].push_back(entry);
                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, 0, bzt);
                                                            entry2.set_trw_flags(tif_x, 0, tif_z);
                                                            entry2.set_src_flags(bx, 0, bz);
                                                            entry2.set_scw_flags(sx_c, 1, sz_c);

                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry(body_idx_src);
                                                            entry.set_tar_flags(bxt, byt, 0);
                                                            entry.set_trw_flags(tif_x, tif_y, 1);
                                                            entry.set_src_flags(bx, by, bz);
                                                            entry.set_scw_flags(sx_c, sy_c, sz_c);

                                                            pair_list[body_idx_tar].push_back(entry);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, byt, 0);
                                                            entry2.set_trw_flags(tif_x, tif_y, 0);
                                                            entry2.set_src_flags(bx, by, 0);
                                                            entry2.set_scw_flags(sx_c, sy_c, 1);

                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        else
                                        {
                                            const bool case1_4interactions_zaxis = in_regz_tar && in_regz_src && dist_z == 0 && !in_rel_regz_src;
                                            const bool case2_4interactions_zaxis = in_regz_tar && in_regz_src && dist_z == dist_1cell && (in_rel_regz_tar != in_rel_regz_src);
                                            const bool case3_4interactions_zaxis = in_regz_tar && in_regz_src && dist_z == dist_2cells && is_zs_fr_scell;

                                            const bool case1_4interactions_yaxis = in_regy_tar && in_regy_src && dist_y == 0 && !in_rel_regy_src;
                                            const bool case2_4interactions_yaxis = in_regy_tar && in_regy_src && dist_y == dist_1cell && (in_rel_regy_tar != in_rel_regy_src);
                                            const bool case3_4interactions_yaxis = in_regy_tar && in_regy_src && dist_y == dist_2cells && is_ys_fr_scell;

                                            const bool case1_4interactions_xaxis = in_regx_tar && in_regx_src && dist_x == 0 && !in_rel_regx_src;
                                            const bool case2_4interactions_xaxis = in_regx_tar && in_regx_src && dist_x == dist_1cell && (in_rel_regx_tar != in_rel_regx_src);
                                            const bool case3_4interactions_xaxis = in_regx_tar && in_regx_src && dist_x == dist_2cells && is_xs_fr_scell;

                                            if (dist_x == dist_2cells)
                                            {
                                                // bx or bxt is always 0
                                                if (in_rel_regx_src && in_rel_regx_tar)
                                                {
                                                    if (case1_4interactions_zaxis)
                                                    {
                                                        if (in_rel_regz_tar)
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, byt, bzt);
                                                            entry1.set_trw_flags(0, tif_y, tif_z);
                                                            entry1.set_src_flags(1, by, 0);
                                                            entry1.set_scw_flags(1, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, byt, bzt);
                                                            entry2.set_trw_flags(1, tif_y, tif_z);
                                                            entry2.set_src_flags(0, by, 0);
                                                            entry2.set_scw_flags(0, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(0, byt, 0);
                                                            entry3.set_trw_flags(0, tif_y, 1);
                                                            entry3.set_src_flags(1, by, 0);
                                                            entry3.set_scw_flags(1, sif_y, 0);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(0, byt, 0);
                                                            entry4.set_trw_flags(1, tif_y, 1);
                                                            entry4.set_src_flags(0, by, 0);
                                                            entry4.set_scw_flags(0, sif_y, 0);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, byt, 0);
                                                            entry1.set_trw_flags(0, tif_y, 1);
                                                            entry1.set_src_flags(1, by, bz);
                                                            entry1.set_scw_flags(1, sif_y, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, byt, 0);
                                                            entry2.set_trw_flags(1, tif_y, 1);
                                                            entry2.set_src_flags(0, by, bz);
                                                            entry2.set_scw_flags(0, sif_y, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(0, byt, 0);
                                                            entry3.set_trw_flags(0, tif_y, 0);
                                                            entry3.set_src_flags(1, by, 0);
                                                            entry3.set_scw_flags(1, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(0, byt, 0);
                                                            entry4.set_trw_flags(1, tif_y, 0);
                                                            entry4.set_src_flags(0, by, 0);
                                                            entry4.set_scw_flags(0, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                    }
                                                    else if (case2_4interactions_zaxis)
                                                    {
                                                        if (in_rel_regz_tar)
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, byt, bzt);
                                                            entry1.set_trw_flags(0, tif_y, tif_z);
                                                            entry1.set_src_flags(1, by, 0);
                                                            entry1.set_scw_flags(1, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, byt, bzt);
                                                            entry2.set_trw_flags(1, tif_y, tif_z);
                                                            entry2.set_src_flags(0, by, 0);
                                                            entry2.set_scw_flags(0, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, byt, 0);
                                                            entry1.set_trw_flags(0, tif_y, 1);
                                                            entry1.set_src_flags(1, by, bz);
                                                            entry1.set_scw_flags(1, sif_y, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, byt, 0);
                                                            entry2.set_trw_flags(1, tif_y, 1);
                                                            entry2.set_src_flags(0, by, bz);
                                                            entry2.set_scw_flags(0, sif_y, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(0, byt, 0);
                                                        entry1.set_trw_flags(0, tif_y, 0);
                                                        entry1.set_src_flags(1, by, 0);
                                                        entry1.set_scw_flags(1, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, byt, 0);
                                                        entry2.set_trw_flags(1, tif_y, 0);
                                                        entry2.set_src_flags(0, by, 0);
                                                        entry2.set_scw_flags(0, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (case3_4interactions_zaxis)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(0, byt, 0);
                                                        entry1.set_trw_flags(0, tif_y, 0);
                                                        entry1.set_src_flags(1, by, 1);
                                                        entry1.set_scw_flags(1, sif_y, 1);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, byt, 0);
                                                        entry2.set_trw_flags(0, tif_y, 1);
                                                        entry2.set_src_flags(1, by, 0);
                                                        entry2.set_scw_flags(1, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry2);

                                                        PairListEntry entry3(body_idx_src);
                                                        entry3.set_tar_flags(0, byt, 0);
                                                        entry3.set_trw_flags(1, tif_y, 0);
                                                        entry3.set_src_flags(0, by, 1);
                                                        entry3.set_scw_flags(0, sif_y, 1);
                                                        pair_list[body_idx_tar].push_back(entry3);

                                                        PairListEntry entry4(body_idx_src);
                                                        entry4.set_tar_flags(0, byt, 0);
                                                        entry4.set_trw_flags(1, tif_y, 1);
                                                        entry4.set_src_flags(0, by, 0);
                                                        entry4.set_scw_flags(0, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry4);
                                                    }
                                                    else if (case1_4interactions_yaxis)
                                                    {
                                                        if (in_rel_regy_tar)
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, byt, bzt);
                                                            entry1.set_trw_flags(0, tif_y, tif_z);
                                                            entry1.set_src_flags(1, 0, bz);
                                                            entry1.set_scw_flags(1, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, byt, bzt);
                                                            entry2.set_trw_flags(1, tif_y, tif_z);
                                                            entry2.set_src_flags(0, 0, bz);
                                                            entry2.set_scw_flags(0, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(0, 0, bzt);
                                                            entry3.set_trw_flags(0, 1, tif_z);
                                                            entry3.set_src_flags(1, 0, bz);
                                                            entry3.set_scw_flags(1, 0, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(0, 0, bzt);
                                                            entry4.set_trw_flags(1, 1, tif_z);
                                                            entry4.set_src_flags(0, 0, bz);
                                                            entry4.set_scw_flags(0, 0, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, 0, bzt);
                                                            entry1.set_trw_flags(0, 1, tif_z);
                                                            entry1.set_src_flags(1, by, bz);
                                                            entry1.set_scw_flags(1, sif_y, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);
                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, 0, bzt);
                                                            entry2.set_trw_flags(1, 1, tif_z);
                                                            entry2.set_src_flags(0, by, bz);
                                                            entry2.set_scw_flags(0, sif_y, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(0, 0, bzt);
                                                            entry3.set_trw_flags(0, 0, tif_z);
                                                            entry3.set_src_flags(1, 0, bz);
                                                            entry3.set_scw_flags(1, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(0, 0, bzt);
                                                            entry4.set_trw_flags(1, 0, tif_z);
                                                            entry4.set_src_flags(0, 0, bz);
                                                            entry4.set_scw_flags(0, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                    }
                                                    else if (case2_4interactions_yaxis)
                                                    {
                                                        if (in_rel_regy_tar)
                                                        {
                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, byt, bzt);
                                                            entry1.set_trw_flags(0, tif_y, tif_z);
                                                            entry1.set_src_flags(1, 0, bz);
                                                            entry1.set_scw_flags(1, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);
                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, byt, bzt);
                                                            entry2.set_trw_flags(1, tif_y, tif_z);
                                                            entry2.set_src_flags(0, 0, bz);
                                                            entry2.set_scw_flags(0, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, 0, bzt);
                                                            entry1.set_trw_flags(0, 1, tif_z);
                                                            entry1.set_src_flags(1, by, bz);
                                                            entry1.set_scw_flags(1, sif_y, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, 0, bzt);
                                                            entry2.set_trw_flags(1, 1, tif_z);
                                                            entry2.set_src_flags(0, by, bz);
                                                            entry2.set_scw_flags(0, sif_y, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(0, 0, bzt);
                                                        entry1.set_trw_flags(0, 0, tif_z);
                                                        entry1.set_src_flags(1, 0, bz);
                                                        entry1.set_scw_flags(1, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, 0, bzt);
                                                        entry2.set_trw_flags(1, 0, tif_z);
                                                        entry2.set_src_flags(0, 0, bz);
                                                        entry2.set_scw_flags(0, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (case3_4interactions_yaxis)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(0, 0, bzt);
                                                        entry1.set_trw_flags(0, 0, tif_z);
                                                        entry1.set_src_flags(1, 1, bz);
                                                        entry1.set_scw_flags(1, 1, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, 0, bzt);
                                                        entry2.set_trw_flags(0, 1, tif_z);
                                                        entry2.set_src_flags(1, 0, bz);
                                                        entry2.set_scw_flags(1, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                        PairListEntry entry3(body_idx_src);
                                                        entry3.set_tar_flags(0, 0, bzt);
                                                        entry3.set_trw_flags(1, 0, tif_z);
                                                        entry3.set_src_flags(0, 1, bz);
                                                        entry3.set_scw_flags(0, 1, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry3);

                                                        PairListEntry entry4(body_idx_src);
                                                        entry4.set_tar_flags(0, 0, bzt);
                                                        entry4.set_trw_flags(1, 1, tif_z);
                                                        entry4.set_src_flags(0, 0, bz);
                                                        entry4.set_scw_flags(0, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry4);
                                                    }
                                                    else
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(1, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(1, byt, bzt);
                                                        entry2.set_trw_flags(1, tif_y, tif_z);
                                                        entry2.set_src_flags(bx, by, bz);
                                                        entry2.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                }
                                                else if (in_rel_regx_src || in_rel_regx_tar)
                                                {
                                                    if (dist_z == 0 && in_regz_src && in_regz_tar && !in_rel_regz_src)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, 0);
                                                        entry1.set_scw_flags(sif_x, sif_y, 1);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, byt, 0);
                                                        entry2.set_trw_flags(tif_x, tif_y, 1);
                                                        entry2.set_src_flags(bx, by, 0);
                                                        entry2.set_scw_flags(sif_x, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_y == 0 && in_regy_src && in_regy_tar && !in_rel_regy_src)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, 0, bz);
                                                        entry1.set_scw_flags(sif_x, 1, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 0, bzt);
                                                        entry2.set_trw_flags(tif_x, 1, tif_z);
                                                        entry2.set_src_flags(bx, 0, bz);
                                                        entry2.set_scw_flags(sif_x, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_z == dist_1cell && in_regz_src && in_regz_tar && in_rel_regz_src != in_rel_regz_tar)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, byt, 0);
                                                        entry2.set_trw_flags(tif_x, tif_y, 0);
                                                        entry2.set_src_flags(bx, by, 0);
                                                        entry2.set_scw_flags(sif_x, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_y == dist_1cell && in_regy_src && in_regy_tar && in_rel_regy_src != in_rel_regy_tar)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 0, bzt);
                                                        entry2.set_trw_flags(tif_x, 0, tif_z);
                                                        entry2.set_src_flags(bx, 0, bz);
                                                        entry2.set_scw_flags(sif_x, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_z == dist_2cells && ((in_rel_regx_src && in_regz_tar && in_rel_regz_tar) || (in_rel_regx_tar && in_regz_src && in_rel_regz_src)))
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, 1);
                                                        entry1.set_trw_flags(tif_x, tif_y, 1);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, byt, 0);
                                                        entry2.set_trw_flags(tif_x, tif_y, 0);
                                                        entry2.set_src_flags(bx, by, 0);
                                                        entry2.set_scw_flags(sif_x, sif_y, 1);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_y == dist_2cells && ((in_rel_regx_src && in_regy_tar && in_rel_regy_tar) || (in_rel_regx_tar && in_regy_src && in_rel_regy_src)))
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, 1, bzt);
                                                        entry1.set_trw_flags(tif_x, 1, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 0, bzt);
                                                        entry2.set_trw_flags(tif_x, 0, tif_z);
                                                        entry2.set_src_flags(bx, 0, bz);
                                                        entry2.set_scw_flags(sif_x, 1, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else
                                                    {
                                                        PairListEntry entry(body_idx_src);
                                                        entry.set_tar_flags(bxt, byt, bzt);
                                                        entry.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry.set_src_flags(bx, by, bz);
                                                        entry.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry);
                                                    }
                                                }
                                            }
                                            else if (dist_y == dist_2cells)
                                            {
                                                if (in_rel_regy_src && in_rel_regy_tar)
                                                {
                                                    // untested
                                                    if (case1_4interactions_zaxis)
                                                    {
                                                        if (in_rel_regz_tar)
                                                        {
                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, 0, bzt);
                                                            entry1.set_trw_flags(tif_x, 0, tif_z);
                                                            entry1.set_src_flags(bx, 1, 0);
                                                            entry1.set_scw_flags(sif_x, 1, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, 0, bzt);
                                                            entry2.set_trw_flags(tif_x, 1, tif_z);
                                                            entry2.set_src_flags(bx, 0, 0);
                                                            entry2.set_scw_flags(sif_x, 0, 1);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(bxt, 0, 0);
                                                            entry3.set_trw_flags(tif_x, 0, 1);
                                                            entry3.set_src_flags(bx, 1, 0);
                                                            entry3.set_scw_flags(sif_x, 1, 0);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(bxt, 0, 0);
                                                            entry4.set_trw_flags(tif_x, 1, 1);
                                                            entry4.set_src_flags(bx, 0, 0);
                                                            entry4.set_scw_flags(sif_x, 0, 0);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, 0, 0);
                                                            entry1.set_trw_flags(tif_x, 0, 1);
                                                            entry1.set_src_flags(bx, 1, bz);
                                                            entry1.set_scw_flags(sif_x, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, 0, 0);
                                                            entry2.set_trw_flags(tif_x, 1, 1);
                                                            entry2.set_src_flags(bx, 0, bz);
                                                            entry2.set_scw_flags(sif_x, 0, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(bxt, 0, 0);
                                                            entry3.set_trw_flags(tif_x, 0, 0);
                                                            entry3.set_src_flags(bx, 1, 0);
                                                            entry3.set_scw_flags(sif_x, 1, 1);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(bxt, 0, 0);
                                                            entry4.set_trw_flags(tif_x, 1, 0);
                                                            entry4.set_src_flags(bx, 0, 0);
                                                            entry4.set_scw_flags(sif_x, 0, 1);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                    }
                                                    else if (case2_4interactions_zaxis)
                                                    {
                                                        if (in_rel_regz_tar)
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, 0, bzt);
                                                            entry1.set_trw_flags(tif_x, 0, tif_z);
                                                            entry1.set_src_flags(bx, 1, 0);
                                                            entry1.set_scw_flags(sif_x, 1, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, 0, bzt);
                                                            entry2.set_trw_flags(tif_x, 1, tif_z);
                                                            entry2.set_src_flags(bx, 0, 0);
                                                            entry2.set_scw_flags(sif_x, 0, 1);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, 0, 0);
                                                            entry1.set_trw_flags(tif_x, 0, 1);
                                                            entry1.set_src_flags(bx, 1, bz);
                                                            entry1.set_scw_flags(sif_x, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, 0, 0);
                                                            entry2.set_trw_flags(tif_x, 1, 1);
                                                            entry2.set_src_flags(bx, 0, bz);
                                                            entry2.set_scw_flags(sif_x, 0, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, 0, 0);
                                                        entry1.set_trw_flags(tif_x, 0, 0);
                                                        entry1.set_src_flags(bx, 1, 0);
                                                        entry1.set_scw_flags(sif_x, 1, 0);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 0, 0);
                                                        entry2.set_trw_flags(tif_x, 1, 0);
                                                        entry2.set_src_flags(bx, 0, 0);
                                                        entry2.set_scw_flags(sif_x, 0, 0);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (case3_4interactions_zaxis)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, 0, 0);
                                                        entry1.set_trw_flags(tif_x, 0, 0);
                                                        entry1.set_src_flags(bx, 1, 1);
                                                        entry1.set_scw_flags(sif_x, 1, 1);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 0, 0);
                                                        entry2.set_trw_flags(tif_x, 0, 1);
                                                        entry2.set_src_flags(bx, 1, 0);
                                                        entry2.set_scw_flags(sif_x, 1, 0);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                        PairListEntry entry3(body_idx_src);
                                                        entry3.set_tar_flags(bxt, 0, 0);
                                                        entry3.set_trw_flags(tif_x, 1, 0);
                                                        entry3.set_src_flags(bx, 0, 1);
                                                        entry3.set_scw_flags(sif_x, 0, 1);
                                                        pair_list[body_idx_tar].push_back(entry3);

                                                        PairListEntry entry4(body_idx_src);
                                                        entry4.set_tar_flags(bxt, 0, 0);
                                                        entry4.set_trw_flags(tif_x, 1, 1);
                                                        entry4.set_src_flags(bx, 0, 0);
                                                        entry4.set_scw_flags(sif_x, 0, 0);
                                                        pair_list[body_idx_tar].push_back(entry4);
                                                    }
                                                    else if (case1_4interactions_xaxis)
                                                    {
                                                        if (in_rel_regx_tar)
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, 0, bzt);
                                                            entry1.set_trw_flags(tif_x, 0, tif_z);
                                                            entry1.set_src_flags(0, 1, bz);
                                                            entry1.set_scw_flags(1, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, 0, bzt);
                                                            entry2.set_trw_flags(tif_x, 1, tif_z);
                                                            entry2.set_src_flags(0, 0, bz);
                                                            entry2.set_scw_flags(1, 0, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(0, 0, bzt);
                                                            entry3.set_trw_flags(1, 0, tif_z);
                                                            entry3.set_src_flags(0, 1, bz);
                                                            entry3.set_scw_flags(0, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(0, 0, bzt);
                                                            entry4.set_trw_flags(1, 1, tif_z);
                                                            entry4.set_src_flags(0, 0, bz);
                                                            entry4.set_scw_flags(0, 0, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                        else
                                                        {
                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, 0, bzt);
                                                            entry1.set_trw_flags(1, 0, tif_z);
                                                            entry1.set_src_flags(bx, 1, bz);
                                                            entry1.set_scw_flags(sif_x, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, 0, bzt);
                                                            entry2.set_trw_flags(1, 1, tif_z);
                                                            entry2.set_src_flags(bx, 0, bz);
                                                            entry2.set_scw_flags(sif_x, 0, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(0, 0, bzt);
                                                            entry3.set_trw_flags(0, 0, tif_z);
                                                            entry3.set_src_flags(0, 1, bz);
                                                            entry3.set_scw_flags(1, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(0, 0, bzt);
                                                            entry4.set_trw_flags(0, 1, tif_z);
                                                            entry4.set_src_flags(0, 0, bz);
                                                            entry4.set_scw_flags(1, 0, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                    }
                                                    else if (case2_4interactions_xaxis)
                                                    {
                                                        if (in_rel_regx_tar)
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, 0, bzt);
                                                            entry1.set_trw_flags(tif_x, 0, tif_z);
                                                            entry1.set_src_flags(0, 1, bz);
                                                            entry1.set_scw_flags(1, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, 0, bzt);
                                                            entry2.set_trw_flags(tif_x, 1, tif_z);
                                                            entry2.set_src_flags(0, 0, bz);
                                                            entry2.set_scw_flags(1, 0, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, 0, bzt);
                                                            entry1.set_trw_flags(1, 0, tif_z);
                                                            entry1.set_src_flags(bx, 1, bz);
                                                            entry1.set_scw_flags(sif_x, 1, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, 0, bzt);
                                                            entry2.set_trw_flags(1, 1, tif_z);
                                                            entry2.set_src_flags(bx, 0, bz);
                                                            entry2.set_scw_flags(sif_x, 0, sif_z);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(0, 0, bzt);
                                                        entry1.set_trw_flags(0, 0, tif_z);
                                                        entry1.set_src_flags(0, 1, bz);
                                                        entry1.set_scw_flags(0, 1, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, 0, bzt);
                                                        entry2.set_trw_flags(0, 1, tif_z);
                                                        entry2.set_src_flags(0, 0, bz);
                                                        entry2.set_scw_flags(0, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (case3_4interactions_xaxis)
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(0, 0, bzt);
                                                        entry1.set_trw_flags(0, 0, tif_z);
                                                        entry1.set_src_flags(1, 1, bz);
                                                        entry1.set_scw_flags(1, 1, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, 0, bzt);
                                                        entry2.set_trw_flags(1, 0, tif_z);
                                                        entry2.set_src_flags(0, 1, bz);
                                                        entry2.set_scw_flags(0, 1, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                        PairListEntry entry3(body_idx_src);
                                                        entry3.set_tar_flags(0, 0, bzt);
                                                        entry3.set_trw_flags(0, 1, tif_z);
                                                        entry3.set_src_flags(1, 0, bz);
                                                        entry3.set_scw_flags(1, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry3);

                                                        PairListEntry entry4(body_idx_src);
                                                        entry4.set_tar_flags(0, 0, bzt);
                                                        entry4.set_trw_flags(1, 1, tif_z);
                                                        entry4.set_src_flags(0, 0, bz);
                                                        entry4.set_scw_flags(0, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry4);
                                                    }
                                                    else
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, 1, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 1, bzt);
                                                        entry2.set_trw_flags(tif_x, 1, tif_z);
                                                        entry2.set_src_flags(bx, by, bz);
                                                        entry2.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                }

                                                else if (in_rel_regy_src || in_rel_regy_tar)
                                                {

                                                    if (dist_z == 0 && in_regz_src && in_regz_tar && !in_rel_regz_src)
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, 0);
                                                        entry1.set_scw_flags(sif_x, sif_y, 1);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, byt, 0);
                                                        entry2.set_trw_flags(tif_x, tif_y, 1);
                                                        entry2.set_src_flags(bx, by, 0);
                                                        entry2.set_scw_flags(sif_x, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_x == 0 && in_regx_src && in_regx_tar && !in_rel_regx_src)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(0, by, bz);
                                                        entry1.set_scw_flags(1, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, byt, bzt);
                                                        entry2.set_trw_flags(1, tif_y, tif_z);
                                                        entry2.set_src_flags(0, by, bz);
                                                        entry2.set_scw_flags(0, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_z == dist_1cell && in_regz_src && in_regz_tar && in_rel_regz_src != in_rel_regz_tar)
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, byt, 0);
                                                        entry2.set_trw_flags(tif_x, tif_y, 0);
                                                        entry2.set_src_flags(bx, by, 0);
                                                        entry2.set_scw_flags(sif_x, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_x == dist_1cell && in_regx_src && in_regx_tar && in_rel_regx_src != in_rel_regx_tar)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, byt, bzt);
                                                        entry2.set_trw_flags(0, tif_y, tif_z);
                                                        entry2.set_src_flags(0, by, bz);
                                                        entry2.set_scw_flags(0, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_z == dist_2cells && ((in_rel_regy_src && in_regz_tar && in_rel_regz_tar) || (in_rel_regy_tar && in_regz_src && in_rel_regz_src)))
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, 1);
                                                        entry1.set_trw_flags(tif_x, tif_y, 1);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, byt, 0);
                                                        entry2.set_trw_flags(tif_x, tif_y, 0);
                                                        entry2.set_src_flags(bx, by, 0);
                                                        entry2.set_scw_flags(sif_x, sif_y, 1);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_x == dist_2cells && ((in_rel_regy_src && in_regx_tar && in_rel_regx_tar) || (in_rel_regy_tar && in_regx_src && in_rel_regx_src)))
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(1, byt, bzt);
                                                        entry1.set_trw_flags(1, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, byt, bzt);
                                                        entry2.set_trw_flags(0, tif_y, tif_z);
                                                        entry2.set_src_flags(0, by, bz);
                                                        entry2.set_scw_flags(1, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);
                                                    }
                                                }
                                            }
                                            else if (dist_z == dist_2cells)
                                            {
                                                if (in_rel_regz_src && in_rel_regz_tar) // not a tested block
                                                {
                                                    if (case1_4interactions_xaxis)
                                                    {
                                                        if (in_rel_regx_tar)
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, byt, 0);
                                                            entry1.set_trw_flags(tif_x, tif_y, 0);
                                                            entry1.set_src_flags(0, by, 1);
                                                            entry1.set_scw_flags(1, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, byt, 0);
                                                            entry2.set_trw_flags(tif_x, tif_y, 1);
                                                            entry2.set_src_flags(0, by, 0);
                                                            entry2.set_scw_flags(1, sif_y, 0);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(0, byt, 0);
                                                            entry3.set_trw_flags(1, tif_y, 0);
                                                            entry3.set_src_flags(0, by, 1);
                                                            entry3.set_scw_flags(0, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(0, byt, 0);
                                                            entry4.set_trw_flags(1, tif_y, 1);
                                                            entry4.set_src_flags(0, by, 0);
                                                            entry4.set_scw_flags(0, sif_y, 0);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, byt, 0);
                                                            entry1.set_trw_flags(1, tif_y, 0);
                                                            entry1.set_src_flags(bx, by, 1);
                                                            entry1.set_scw_flags(sif_x, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, byt, 0);
                                                            entry2.set_trw_flags(1, tif_y, 1);
                                                            entry2.set_src_flags(bx, by, 0);
                                                            entry2.set_scw_flags(sif_x, sif_y, 0);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(0, byt, 0);
                                                            entry3.set_trw_flags(0, tif_y, 0);
                                                            entry3.set_src_flags(0, by, 1);
                                                            entry3.set_scw_flags(1, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(0, byt, 0);
                                                            entry4.set_trw_flags(0, tif_y, 1);
                                                            entry4.set_src_flags(0, by, 0);
                                                            entry4.set_scw_flags(1, sif_y, 0);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                    }

                                                    else if (case2_4interactions_xaxis)
                                                    {
                                                        if (in_rel_regx_tar)
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, byt, 0);
                                                            entry1.set_trw_flags(tif_x, tif_y, 0);
                                                            entry1.set_src_flags(0, by, 1);
                                                            entry1.set_scw_flags(1, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, byt, 0);
                                                            entry2.set_trw_flags(tif_x, tif_y, 1);
                                                            entry2.set_src_flags(0, by, 0);
                                                            entry2.set_scw_flags(1, sif_y, 0);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(0, byt, 0);
                                                            entry1.set_trw_flags(1, tif_y, 0);
                                                            entry1.set_src_flags(bx, by, 1);
                                                            entry1.set_scw_flags(sif_x, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(0, byt, 0);
                                                            entry2.set_trw_flags(1, tif_y, 1);
                                                            entry2.set_src_flags(bx, by, 0);
                                                            entry2.set_scw_flags(sif_x, sif_y, 0);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(0, byt, 0);
                                                        entry1.set_trw_flags(0, tif_y, 0);
                                                        entry1.set_src_flags(0, by, 1);
                                                        entry1.set_scw_flags(0, sif_y, 1);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, byt, 0);
                                                        entry2.set_trw_flags(0, tif_y, 1);
                                                        entry2.set_src_flags(0, by, 0);
                                                        entry2.set_scw_flags(0, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (case3_4interactions_xaxis) // corrected
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(0, byt, 0);
                                                        entry1.set_trw_flags(0, tif_y, 0);
                                                        entry1.set_src_flags(1, by, 1);
                                                        entry1.set_scw_flags(1, sif_y, 1);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, byt, 0);
                                                        entry2.set_trw_flags(1, tif_y, 0);
                                                        entry2.set_src_flags(0, by, 1);
                                                        entry2.set_scw_flags(0, sif_y, 1);
                                                        pair_list[body_idx_tar].push_back(entry2);

                                                        PairListEntry entry3(body_idx_src);
                                                        entry3.set_tar_flags(0, byt, 0);
                                                        entry3.set_trw_flags(0, tif_y, 1);
                                                        entry3.set_src_flags(1, by, 0);
                                                        entry3.set_scw_flags(1, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry3);

                                                        PairListEntry entry4(body_idx_src);
                                                        entry4.set_tar_flags(0, byt, 0);
                                                        entry4.set_trw_flags(1, tif_y, 1);
                                                        entry4.set_src_flags(0, by, 0);
                                                        entry4.set_scw_flags(0, sif_y, 0);
                                                        pair_list[body_idx_tar].push_back(entry4);
                                                    }
                                                    else if (case1_4interactions_yaxis)
                                                    {
                                                        if (in_rel_regy_tar)
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, byt, 0);
                                                            entry1.set_trw_flags(tif_x, tif_y, 0);
                                                            entry1.set_src_flags(bx, 0, 1);
                                                            entry1.set_scw_flags(sif_x, 1, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, byt, 0);
                                                            entry2.set_trw_flags(tif_x, tif_y, 1);
                                                            entry2.set_src_flags(bx, 0, 0);
                                                            entry2.set_scw_flags(sif_x, 1, 0);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(bxt, 0, 0);
                                                            entry3.set_trw_flags(tif_x, 1, 0);
                                                            entry3.set_src_flags(bx, 0, 1);
                                                            entry3.set_scw_flags(sif_x, 0, 1);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(bxt, 0, 0);
                                                            entry4.set_trw_flags(tif_x, 1, 1);
                                                            entry4.set_src_flags(bx, 0, 0);
                                                            entry4.set_scw_flags(sif_x, 0, 0);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                        else
                                                        {

                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, 0, 0);
                                                            entry1.set_trw_flags(tif_x, 1, 0);
                                                            entry1.set_src_flags(bx, by, 1);
                                                            entry1.set_scw_flags(sif_x, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, 0, 0);
                                                            entry2.set_trw_flags(tif_x, 1, 1);
                                                            entry2.set_src_flags(bx, by, 0);
                                                            entry2.set_scw_flags(sif_x, sif_y, 0);
                                                            pair_list[body_idx_tar].push_back(entry2);

                                                            PairListEntry entry3(body_idx_src);
                                                            entry3.set_tar_flags(bxt, 0, 0);
                                                            entry3.set_trw_flags(tif_x, 0, 0);
                                                            entry3.set_src_flags(bx, 0, 1);
                                                            entry3.set_scw_flags(sif_x, 1, 1);
                                                            pair_list[body_idx_tar].push_back(entry3);

                                                            PairListEntry entry4(body_idx_src);
                                                            entry4.set_tar_flags(bxt, 0, 0);
                                                            entry4.set_trw_flags(tif_x, 0, 1);
                                                            entry4.set_src_flags(bx, 0, 0);
                                                            entry4.set_scw_flags(sif_x, 1, 0);
                                                            pair_list[body_idx_tar].push_back(entry4);
                                                        }
                                                    }
                                                    else if (case2_4interactions_yaxis)
                                                    {
                                                        if (in_rel_regy_tar)
                                                        {
                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, byt, 0);
                                                            entry1.set_trw_flags(tif_x, tif_y, 0);
                                                            entry1.set_src_flags(bx, 0, 1);
                                                            entry1.set_scw_flags(sif_x, 1, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, byt, 0);
                                                            entry2.set_trw_flags(tif_x, tif_y, 1);
                                                            entry2.set_src_flags(bx, 0, 0);
                                                            entry2.set_scw_flags(sif_x, 1, 0);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                        else
                                                        {
                                                            PairListEntry entry1(body_idx_src);
                                                            entry1.set_tar_flags(bxt, 0, 0);
                                                            entry1.set_trw_flags(tif_x, 1, 0);
                                                            entry1.set_src_flags(bx, by, 1);
                                                            entry1.set_scw_flags(sif_x, sif_y, 1);
                                                            pair_list[body_idx_tar].push_back(entry1);

                                                            PairListEntry entry2(body_idx_src);
                                                            entry2.set_tar_flags(bxt, 0, 0);
                                                            entry2.set_trw_flags(tif_x, 1, 1);
                                                            entry2.set_src_flags(bx, by, 0);
                                                            entry2.set_scw_flags(sif_x, sif_y, 0);
                                                            pair_list[body_idx_tar].push_back(entry2);
                                                        }
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, 0, 0);
                                                        entry1.set_trw_flags(tif_x, 0, 0);
                                                        entry1.set_src_flags(bx, 0, 1);
                                                        entry1.set_scw_flags(sif_x, 0, 1);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 0, 0);
                                                        entry2.set_trw_flags(tif_x, 0, 1);
                                                        entry2.set_src_flags(bx, 0, 0);
                                                        entry2.set_scw_flags(sif_x, 0, 0);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (case3_4interactions_yaxis) // corrected
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, 0, 0);
                                                        entry1.set_trw_flags(tif_x, 0, 0);
                                                        entry1.set_src_flags(bx, 1, 1);
                                                        entry1.set_scw_flags(sif_x, 1, 1);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 0, 0);
                                                        entry2.set_trw_flags(tif_x, 1, 0);
                                                        entry2.set_src_flags(bx, 0, 1);
                                                        entry2.set_scw_flags(sif_x, 0, 1);
                                                        pair_list[body_idx_tar].push_back(entry2);

                                                        PairListEntry entry3(body_idx_src);
                                                        entry3.set_tar_flags(bxt, 0, 0);
                                                        entry3.set_trw_flags(tif_x, 0, 1);
                                                        entry3.set_src_flags(bx, 1, 0);
                                                        entry3.set_scw_flags(sif_x, 1, 0);
                                                        pair_list[body_idx_tar].push_back(entry3);

                                                        PairListEntry entry4(body_idx_src);
                                                        entry4.set_tar_flags(bxt, 0, 0);
                                                        entry4.set_trw_flags(tif_x, 1, 1);
                                                        entry4.set_src_flags(bx, 0, 0);
                                                        entry4.set_scw_flags(sif_x, 0, 0);
                                                        pair_list[body_idx_tar].push_back(entry4);
                                                    }

                                                    else
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, 1);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, byt, 1);
                                                        entry2.set_trw_flags(tif_x, tif_y, 1);
                                                        entry2.set_src_flags(bx, by, bz);
                                                        entry2.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                }
                                                else if (in_rel_regz_src || in_rel_regz_tar)
                                                {
                                                    if (dist_y == 0 && in_regy_src && in_regy_tar && !in_rel_regy_src)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, 0, bz);
                                                        entry1.set_scw_flags(sif_x, 1, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 0, bzt);
                                                        entry2.set_trw_flags(tif_x, 1, tif_z);
                                                        entry2.set_src_flags(bx, 0, bz);
                                                        entry2.set_scw_flags(sif_x, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_x == 0 && in_regx_src && in_regx_tar && !in_rel_regx_src)
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(0, by, bz);
                                                        entry1.set_scw_flags(1, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, byt, bzt);
                                                        entry2.set_trw_flags(1, tif_y, tif_z);
                                                        entry2.set_src_flags(0, by, bz);
                                                        entry2.set_scw_flags(0, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_y == dist_1cell && in_regy_src && in_regy_tar && in_rel_regy_src != in_rel_regy_tar)
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 0, bzt);
                                                        entry2.set_trw_flags(tif_x, 0, tif_z);
                                                        entry2.set_src_flags(bx, 0, bz);
                                                        entry2.set_scw_flags(sif_x, 0, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_x == dist_1cell && in_regx_src && in_regx_tar && in_rel_regx_src != in_rel_regx_tar)
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, byt, bzt);
                                                        entry2.set_trw_flags(0, tif_y, tif_z);
                                                        entry2.set_src_flags(0, by, bz);
                                                        entry2.set_scw_flags(0, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_y == dist_2cells && ((in_rel_regz_src && in_regy_tar && in_rel_regy_tar) || (in_rel_regz_tar && in_regy_src && in_rel_regy_src)))
                                                    {

                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, 1, bzt);
                                                        entry1.set_trw_flags(tif_x, 1, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(bxt, 0, bzt);
                                                        entry2.set_trw_flags(tif_x, 0, tif_z);
                                                        entry2.set_src_flags(bx, 0, bz);
                                                        entry2.set_scw_flags(sif_x, 1, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else if (dist_x == dist_2cells && ((in_rel_regz_src && in_regx_tar && in_rel_regx_tar) || (in_rel_regz_tar && in_regx_src && in_rel_regx_src)))
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(1, byt, bzt);
                                                        entry1.set_trw_flags(1, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);

                                                        PairListEntry entry2(body_idx_src);
                                                        entry2.set_tar_flags(0, byt, bzt);
                                                        entry2.set_trw_flags(0, tif_y, tif_z);
                                                        entry2.set_src_flags(0, by, bz);
                                                        entry2.set_scw_flags(1, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry2);
                                                    }
                                                    else
                                                    {
                                                        PairListEntry entry1(body_idx_src);
                                                        entry1.set_tar_flags(bxt, byt, bzt);
                                                        entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                        entry1.set_src_flags(bx, by, bz);
                                                        entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                        pair_list[body_idx_tar].push_back(entry1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    const real interaction_region_tcell = cell.radius * 5;     // two and half cell
                                    const real interaction_region_scell = adj_cell.radius * 5; // two and half cell

                                    const real interaction_region_tcelli = cell.radius * 3;     // one and half cell
                                    const real interaction_region_scelli = adj_cell.radius * 3; // one and half cell

                                    const bool bxti = !in_regx_tar || (fabs(body_tar.x[0] - adj_cell.center[0]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scelli);
                                    const bool byti = !in_regy_tar || (fabs(body_tar.x[1] - adj_cell.center[1]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scelli);
                                    const bool bzti = !in_regz_tar || (fabs(body_tar.x[2] - adj_cell.center[2]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scelli);

                                    const bool bxi = !in_regx_src || (fabs(body_src.x[0] - cell.center[0]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_tcelli);
                                    const bool byi = !in_regy_src || (fabs(body_src.x[1] - cell.center[1]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_tcelli);
                                    const bool bzi = !in_regz_src || (fabs(body_src.x[2] - cell.center[2]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_tcelli);

                                    const bool sif_x = !(dist_x == dist_3cells && in_rel_regx_src);
                                    const bool sif_y = !(dist_y == dist_3cells && in_rel_regy_src);
                                    const bool sif_z = !(dist_z == dist_3cells && in_rel_regz_src);

                                    const bool tif_x = !(dist_x == dist_3cells && in_rel_regx_tar);
                                    const bool tif_y = !(dist_y == dist_3cells && in_rel_regy_tar);
                                    const bool tif_z = !(dist_z == dist_3cells && in_rel_regz_tar);

                                    // Check if the source particle is within
                                    // the x-range of the target cell
                                    const bool srcp_in_tarc_x = fabs(body_src.x[0] - cell.center[0]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the source particle is within
                                    // the y-range of the target cell
                                    const bool srcp_in_tarc_y = fabs(body_src.x[1] - cell.center[1]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the source particle is within
                                    // the z-range of the target cell
                                    const bool srcp_in_tarc_z = fabs(body_src.x[2] - cell.center[2]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the target particle is within
                                    // the x-range of the source cell
                                    const bool tarp_in_srcc_x = fabs(body_tar.x[0] - adj_cell.center[0]) <= interaction_region_scell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the target particle is within
                                    // the y-range of the source cell
                                    const bool tarp_in_srcc_y = fabs(body_tar.x[1] - adj_cell.center[1]) <= interaction_region_scell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the target particle is within
                                    // the z-range of the source cell
                                    const bool tarp_in_srcc_z = fabs(body_tar.x[2] - adj_cell.center[2]) <= interaction_region_scell + fmm_weights_eval_.getRegAlpha();

                                    // Determine if the source particle is
                                    // entirely within the range of the target
                                    // cell
                                    const bool srcp_in_tarc = srcp_in_tarc_x && srcp_in_tarc_y && srcp_in_tarc_z;

                                    // Determine if the target particle is
                                    // entirely within the range of the source
                                    // cell
                                    const bool tarp_in_srcc = tarp_in_srcc_x && tarp_in_srcc_y && tarp_in_srcc_z;

                                    if ((srcp_in_tarc && (sif_x == false || sif_y == false || sif_z == false)) && (tarp_in_srcc && (tif_x == false || tif_y == false || tif_z == false)))
                                    {

                                        const bool is_xaxis_relevant = (dist_x == dist_2cells && (in_rel_regx_tar || in_rel_regx_src));
                                        const bool is_yaxis_relevant = (dist_y == dist_2cells && (in_rel_regy_tar || in_rel_regy_src));
                                        const bool is_zaxis_relevant = (dist_z == dist_2cells && (in_rel_regz_tar || in_rel_regz_src));

                                        const bool is_xaxis_relevant_src = (dist_x == dist_2cells && in_rel_regx_src);
                                        const bool is_yaxis_relevant_src = (dist_y == dist_2cells && in_rel_regy_src);
                                        const bool is_zaxis_relevant_src = (dist_z == dist_2cells && in_rel_regz_src);

                                        const bool is_xaxis_relevant_tar = (dist_x == dist_2cells && in_rel_regx_tar);
                                        const bool is_yaxis_relevant_tar = (dist_y == dist_2cells && in_rel_regy_tar);
                                        const bool is_zaxis_relevant_tar = (dist_z == dist_2cells && in_rel_regz_tar);

                                        if (dist_x == dist_3cells && in_rel_regx_src && in_rel_regx_tar)
                                        {
                                            if (dist_y < dist_2cells && dist_z < dist_2cells)
                                            {
                                                PairListEntry entry1(body_idx_src);
                                                entry1.set_tar_flags(bxti, byti, bzti);
                                                entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                entry1.set_src_flags(bxi, byi, bzi);
                                                entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                pair_list[body_idx_tar].push_back(entry1);
                                            }
                                            else if (dist_y < dist_2cells && is_zaxis_relevant)
                                            {
                                                if (is_zaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, bzti);
                                                    entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                    entry1.set_src_flags(bxi, byi, 0);
                                                    entry1.set_scw_flags(sif_x, sif_y, 0);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, 0);
                                                    entry1.set_trw_flags(tif_x, tif_y, 0);
                                                    entry1.set_src_flags(bxi, byi, bzi);
                                                    entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                            }
                                            else if (dist_z < dist_2cells && is_yaxis_relevant)
                                            {
                                                if (is_yaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, bzti);
                                                    entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                    entry1.set_src_flags(bxi, 0, bzi);
                                                    entry1.set_scw_flags(sif_x, 0, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, 0, bzti);
                                                    entry1.set_trw_flags(tif_x, 0, tif_z);
                                                    entry1.set_src_flags(bxi, byi, bzi);
                                                    entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                            }
                                            else if (is_yaxis_relevant && is_zaxis_relevant)
                                            {
                                                if (is_yaxis_relevant_src && is_zaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, bzti);
                                                    entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                    entry1.set_src_flags(bxi, 0, 0);
                                                    entry1.set_scw_flags(sif_x, 0, 0);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else if (is_yaxis_relevant_src && is_zaxis_relevant_tar)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, 0);
                                                    entry1.set_trw_flags(tif_x, tif_y, 0);
                                                    entry1.set_src_flags(bxi, 0, bzi);
                                                    entry1.set_scw_flags(sif_x, 0, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else if (is_yaxis_relevant_tar && is_zaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, 0, bzti);
                                                    entry1.set_trw_flags(tif_x, 0, tif_z);
                                                    entry1.set_src_flags(bxi, byi, 0);
                                                    entry1.set_scw_flags(sif_x, sif_y, 0);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, 0, 0);
                                                    entry1.set_trw_flags(tif_x, 0, 0);
                                                    entry1.set_src_flags(bxi, byi, bzi);
                                                    entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                            }
                                        }
                                        else if (dist_y == dist_3cells && in_rel_regy_src && in_rel_regy_tar)
                                        {
                                            if (dist_x < dist_2cells && dist_z < dist_2cells)
                                            {
                                                PairListEntry entry1(body_idx_src);
                                                entry1.set_tar_flags(bxti, byti, bzti);
                                                entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                entry1.set_src_flags(bxi, byi, bzi);
                                                entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                pair_list[body_idx_tar].push_back(entry1);
                                            }
                                            else if (dist_x < dist_2cells && is_zaxis_relevant)
                                            {
                                                if (is_zaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, bzti);
                                                    entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                    entry1.set_src_flags(bxi, byi, 0);
                                                    entry1.set_scw_flags(sif_x, sif_y, 0);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, 0);
                                                    entry1.set_trw_flags(tif_x, tif_y, 0);
                                                    entry1.set_src_flags(bxi, byi, bzi);
                                                    entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                            }
                                            else if (dist_z < dist_2cells && is_xaxis_relevant)
                                            {
                                                if (is_xaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, bzti);
                                                    entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                    entry1.set_src_flags(0, byi, bzi);
                                                    entry1.set_scw_flags(0, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(0, byti, bzti);
                                                    entry1.set_trw_flags(0, tif_y, tif_z);
                                                    entry1.set_src_flags(bxi, byi, bzi);
                                                    entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                            }
                                            else if (is_xaxis_relevant && is_zaxis_relevant)
                                            {
                                                if (is_xaxis_relevant_src && is_zaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, bzti);
                                                    entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                    entry1.set_src_flags(0, byi, 0);
                                                    entry1.set_scw_flags(0, sif_y, 0);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else if (is_xaxis_relevant_tar && is_zaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(0, byti, bzti);
                                                    entry1.set_trw_flags(0, tif_y, tif_z);
                                                    entry1.set_src_flags(bxi, byi, 0);
                                                    entry1.set_scw_flags(sif_x, sif_y, 0);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else if (is_xaxis_relevant_src && is_zaxis_relevant_tar)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, 0);
                                                    entry1.set_trw_flags(tif_x, tif_y, 0);
                                                    entry1.set_src_flags(0, byi, bzi);
                                                    entry1.set_scw_flags(0, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(0, byti, 0);
                                                    entry1.set_trw_flags(0, tif_y, 0);
                                                    entry1.set_src_flags(bxi, byi, bzi);
                                                    entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                            }
                                        }
                                        else if (dist_z == dist_3cells && in_rel_regz_src && in_rel_regz_tar)
                                        {

                                            if (dist_x < dist_2cells && dist_y < dist_2cells)
                                            {
                                                PairListEntry entry1(body_idx_src);
                                                entry1.set_tar_flags(bxti, byti, bzti);
                                                entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                entry1.set_src_flags(bxi, byi, bzi);
                                                entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                pair_list[body_idx_tar].push_back(entry1);
                                            }
                                            else if (is_xaxis_relevant && dist_y < dist_2cells)
                                            {
                                                if (is_xaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, bzti);
                                                    entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                    entry1.set_src_flags(0, byi, bzi);
                                                    entry1.set_scw_flags(0, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(0, byti, bzti);
                                                    entry1.set_trw_flags(0, tif_y, tif_z);
                                                    entry1.set_src_flags(bxi, byi, bzi);
                                                    entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                            }
                                            else if (is_yaxis_relevant && dist_x < dist_2cells)
                                            {
                                                if (is_yaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, bzti);
                                                    entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                    entry1.set_src_flags(bxi, 0, bzi);
                                                    entry1.set_scw_flags(sif_x, 0, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, 0, bzti);
                                                    entry1.set_trw_flags(tif_x, 0, tif_z);
                                                    entry1.set_src_flags(bxi, byi, bzi);
                                                    entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                            }
                                            else if (is_xaxis_relevant && is_yaxis_relevant)
                                            {
                                                if (is_xaxis_relevant_src && is_yaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, byti, bzti);
                                                    entry1.set_trw_flags(tif_x, tif_y, tif_z);
                                                    entry1.set_src_flags(0, 0, bzi);
                                                    entry1.set_scw_flags(0, 0, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else if (is_xaxis_relevant_src && is_yaxis_relevant_tar)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(bxti, 0, bzti);
                                                    entry1.set_trw_flags(tif_x, 0, tif_z);
                                                    entry1.set_src_flags(0, byi, bzi);
                                                    entry1.set_scw_flags(0, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else if (is_xaxis_relevant_tar && is_yaxis_relevant_src)
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(0, byti, bzti);
                                                    entry1.set_trw_flags(0, tif_y, tif_z);
                                                    entry1.set_src_flags(bxi, 0, bzi);
                                                    entry1.set_scw_flags(sif_x, 0, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
                                                }
                                                else
                                                {
                                                    PairListEntry entry1(body_idx_src);
                                                    entry1.set_tar_flags(0, 0, bzti);
                                                    entry1.set_trw_flags(0, 0, tif_z);
                                                    entry1.set_src_flags(bxi, byi, bzi);
                                                    entry1.set_scw_flags(sif_x, sif_y, sif_z);
                                                    pair_list[body_idx_tar].push_back(entry1);
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

    std::unordered_map<std::string, std::array<real, 4>> st;
    std::ofstream fout1("log_file_my_interactons.txt");

    for (size_t i = 0; i < bodies_all_.size(); i++)
    {
        const FBody &body_tar = bodies_all_[i];
        const real xt = body_tar.x[0];
        const real yt = body_tar.x[1];
        const real zt = body_tar.x[2];
        const RVec wtar_ws = w_per_atom[i];

        real pj_effective = 0.0;
        real fxj_effective = 0.0, fyj_effective = 0.0, fzj_effective = 0.0;
        int ix = 0;

        for (auto &entry : pair_list[i])
        { 

            const int body_idx_src = entry.body_idx_src;
            gmx::fmm::FBody &asrc = bodies_all_[body_idx_src];

            const bool bx_src = entry.bx_src;
            const bool by_src = entry.by_src;
            const bool bz_src = entry.bz_src;

            const bool bx_tar = entry.bx_tar;
            const bool by_tar = entry.by_tar;
            const bool bz_tar = entry.bz_tar;

            const bool is_within_sx = entry.sx_within;
            const bool is_within_sy = entry.sy_within;
            const bool is_within_sz = entry.sz_within;

            const bool is_within_tx = entry.tx_within;
            const bool is_within_ty = entry.ty_within;
            const bool is_within_tz = entry.tz_within;

            const RVec wsrc_ws = w_per_atom[body_idx_src];

            const real xs = asrc.x[0];
            const real ys = asrc.x[1];
            const real zs = asrc.x[2];

            const real dx = xt - xs;
            const real dy = yt - ys;
            const real dz = zt - zs;

            real wsrc_x = bx_src == 1 ? 1 : (is_within_sx == 1 ? wsrc_ws[0] : 1 - wsrc_ws[0]);
            real wsrc_y = by_src == 1 ? 1 : (is_within_sy == 1 ? wsrc_ws[1] : 1 - wsrc_ws[1]);
            real wsrc_z = bz_src == 1 ? 1 : (is_within_sz == 1 ? wsrc_ws[2] : 1 - wsrc_ws[2]);
            real wsrc = wsrc_x * wsrc_y * wsrc_z;

            real wtar_x = bx_tar == 1 ? 1 : (is_within_tx == 1 ? wtar_ws[0] : 1 - wtar_ws[0]);
            real wtar_y = by_tar == 1 ? 1 : (is_within_ty == 1 ? wtar_ws[1] : 1 - wtar_ws[1]);
            real wtar_z = bz_tar == 1 ? 1 : (is_within_tz == 1 ? wtar_ws[2] : 1 - wtar_ws[2]);
            const real wtar = wtar_x * wtar_y * wtar_z;

            // std::cout << asrc.x << "--" << body_tar.x << "--" << wsrc << "--" << wtar << "" << "\n";
            fout1 << asrc.x << "--" << body_tar.x << "--" << wsrc << "--" << wtar << "" << "\n";
            std::string key = std::to_string(body_idx_src) + "--" + std::to_string(i);
            // std::cout << wsrc_ws << "--" << wtar_ws << std::endl;

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

            if (st.find(key) != st.end())
            {
                st[key][0] += wtar * qinvr;
                st[key][1] += wtar * fxj;
                st[key][2] += wtar * fyj;
                st[key][3] += wtar * fzj;
            }
            else
            {
                st[key] = {wtar * qinvr, wtar * fxj, wtar * fyj, wtar * fzj};
            }

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

    fout1.close();

    std::ofstream fout("log_file_my_fps.txt");
    // Iterate through the map
    for (const auto &[key, values] : st)
    {
        // Find the position of the delimiter "--"
        size_t delimiter_pos = key.find("--");
        if (delimiter_pos == std::string::npos)
        {
            std::cerr << "Invalid key format for key: " << key << std::endl;
            continue;
        }

        // Extract the two integers from the key
        std::string part1 = key.substr(0, delimiter_pos);  // First part
        std::string part2 = key.substr(delimiter_pos + 2); // Second part
        int bidx_src = std::stoi(part1);
        int bidx_tar = std::stoi(part2);

        fout << bodies_all_[bidx_src].x << "--" << bodies_all_[bidx_tar].x << "--[" << values[0] << "," << values[1] << "," << values[2] << "," << values[3] << "]\n";
    }
    fout.close();

    return forces_and_potentials;
}

void gmx::fmm::FMMDirectInteractions::recompute_weights() { compute_weights_(); }

void gmx::fmm::FMMDirectInteractions::rebuild_and_reprocess_tree() { fmm_direct_interactions_tree_.rebuild_and_reprocess_tree(); }
