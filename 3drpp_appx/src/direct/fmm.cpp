#include "fmm.h"
#include <fstream>
#include <iomanip>

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

    pair_list_tif_within.clear();
    pair_list_tif_within.resize(bodies_all_.size());

    pair_list_sif_within.clear();
    pair_list_sif_within.resize(bodies_all_.size());

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
                        // pair_list[body_idx_tar].push_back(body_idx_src);
                        // pair_list_bxyz_src[body_idx_tar].push_back({1, 1, 1});
                        // pair_list_bxyz_tar[body_idx_tar].push_back({1, 1, 1});
                        // pair_list_interaction_type_src[body_idx_tar].push_back(1);
                        // pair_list_interaction_type_tar[body_idx_tar].push_back(1);
                    }
                    else
                    {

                        bool x_dir_flag = true;
                        bool y_dir_flag = true;
                        bool z_dir_flag = true;

                        real r1x = body_tar.x[0] - cell.center[0];
                        real r2x = body_src.x[0] - cell.center[0];

                        // both particles are in opposite direction of x within the cell
                        if (r1x * r2x < 0)
                        {
                            x_dir_flag = false;
                        }

                        // both particles are in opposite direction of y within the cell
                        real r1y = body_tar.x[1] - cell.center[1];
                        real r2y = body_src.x[1] - cell.center[1];
                        if (r1y * r2y < 0)
                        {
                            y_dir_flag = false;
                        }

                        // both particles are in opposite direction of z within the cell
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
                        // pair_list[body_idx_tar].push_back(body_idx_src);

                        // pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});
                        // pair_list_bxyz_tar[body_idx_tar].push_back({bxts, byts, bzts});

                        // pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                        // pair_list_bxyz_src[body_idx_tar].push_back({1, 1, 1});

                        // indirect weights outside the boundary of the cell
                        if (bxts == false || byts == false || bzts == false)
                        {
                            // pair_list[body_idx_tar].push_back(body_idx_src);
                            // pair_list_bxyz_src[body_idx_tar].push_back({bxts, byts, bzts});
                            // pair_list_bxyz_tar[body_idx_tar].push_back({bxts, byts, bzts});
                            // pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                            // pair_list_tif_within[body_idx_tar].push_back({0, 0, 0});
                        }
                    }
                }
            }

            const real interaction_region_tcell = cell.radius * 3; // one and half cell
            const real width_of_tcell = cell.radius * 2;
            const RVec wst = w_per_atom[body_idx_tar];

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
                            const real interaction_region_scell = adj_cell.radius * 3;
                            short num_away = 1; // Default to "one away"
                            if (std::abs(dx) == 3 || std::abs(dy) == 3 || std::abs(dz) == 3)
                            {
                                num_away = 3; // If any direction is "three away"
                            }
                            else if (std::abs(dx) == 2 || std::abs(dy) == 2 || std::abs(dz) == 2)
                            {
                                num_away = 2; // If any direction is "two away"
                            }

                            bool bxt = wst[0] != 1
                                           ? (fabs(body_tar.x[0] - adj_cell.center[0]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scell ? 1 : 0)
                                           : 1;

                            bool byt = wst[1] != 1
                                           ? (fabs(body_tar.x[1] - adj_cell.center[1]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scell ? 1 : 0)
                                           : 1;

                            bool bzt = wst[2] != 1
                                           ? (fabs(body_tar.x[2] - adj_cell.center[2]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_scell ? 1 : 0)
                                           : 1;

                            if (num_away == 1)
                            {
                                // one and half cell
                                for (const int &body_idx_src : adj_cell.bodiesIndices)
                                {
                                    const FBody &body_src = bodies_all_[body_idx_src];
                                    if (body_idx_tar != body_idx_src)
                                    {
                                        const RVec ws = w_per_atom[body_idx_src];
                                        bool bx =
                                            ws[0] < 1
                                                ? (fabs(body_src.x[0] - cell.center[0]) + fmm_weights_eval_.getRegAlpha() > interaction_region_tcell ? 0 : 1)
                                                : 1;

                                        bool by =
                                            ws[1] < 1
                                                ? (fabs(body_src.x[1] - cell.center[1]) + fmm_weights_eval_.getRegAlpha() > interaction_region_tcell ? 0 : 1)
                                                : 1;

                                        bool bz =
                                            ws[2] < 1
                                                ? (fabs(body_src.x[2] - cell.center[2]) + fmm_weights_eval_.getRegAlpha() > interaction_region_tcell ? 0 : 1)
                                                : 1;

                                        // pair_list[body_idx_tar].push_back(body_idx_src);
                                        // pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                        // pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                        // pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                        // pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});
                                    }
                                }
                            }

                            else if (num_away == 2) // num_away == 2
                            {

                                const real rtx_tc = body_tar.x[0] - cell.center[0];
                                const real rty_tc = body_tar.x[1] - cell.center[1];
                                const real rtz_tc = body_tar.x[2] - cell.center[2];

                                const real rtx_sc = body_tar.x[0] - adj_cell.center[0];
                                const real rty_sc = body_tar.x[1] - adj_cell.center[1];
                                const real rtz_sc = body_tar.x[2] - adj_cell.center[2];

                                const real dist_2cells = cell.radius * 4;

                                const real dist_x = fabs(adj_cell.center[0] - cell.center[0]);
                                const real dist_y = fabs(adj_cell.center[1] - cell.center[1]);
                                const real dist_z = fabs(adj_cell.center[2] - cell.center[2]);

                                bool tif_x = true;
                                bool tif_y = true;
                                bool tif_z = true;

                                const bool in_regx_tar = fabs(body_tar.x[0] - cell.center[0]) + fmm_weights_eval_.getRegAlpha() > cell.radius;
                                const bool in_regy_tar = fabs(body_tar.x[1] - cell.center[1]) + fmm_weights_eval_.getRegAlpha() > cell.radius;
                                const bool in_regz_tar = fabs(body_tar.x[2] - cell.center[2]) + fmm_weights_eval_.getRegAlpha() > cell.radius;

                                const bool in_rel_regx_tar = in_regx_tar && rtx_tc * rtx_sc < 0;
                                const bool in_rel_regy_tar = in_regy_tar && rty_tc * rty_sc < 0;
                                const bool in_rel_regz_tar = in_regz_tar && rtz_tc * rtz_sc < 0;

                                // if (dist_x == dist_2cells && dist_y != dist_2cells && dist_z != dist_2cells)
                                // {
                                //     if (in_rel_regx_tar)
                                //     {
                                //         tif_x = false;
                                //     }
                                // }
                                // else if (dist_y == dist_2cells && dist_x != dist_2cells && dist_z != dist_2cells)
                                // {
                                //     if (in_rel_regy_tar)
                                //     {
                                //         tif_y = false;
                                //     }
                                // }
                                // else if (dist_z == dist_2cells && dist_x != dist_2cells && dist_y != dist_2cells)
                                // {
                                //     if (in_rel_regz_tar)
                                //     {
                                //         tif_z = false;
                                //     }
                                // }
                                // else if (dist_x == dist_2cells && dist_y == dist_2cells && dist_z != dist_2cells)
                                // {
                                //     if (in_rel_regx_tar && in_rel_regy_tar)
                                //     {
                                //         tif_x = false;
                                //         tif_y = false;
                                //     }
                                // }
                                // else if (dist_y == dist_2cells && dist_z == dist_2cells && dist_x != dist_2cells)
                                // {
                                //     if (in_rel_regy_tar && in_rel_regz_tar)
                                //     {
                                //         tif_y = false;
                                //         tif_z = false;
                                //     }
                                // }
                                // else if (dist_x == dist_2cells && dist_z == dist_2cells && dist_y != dist_2cells)
                                // {
                                //     if (in_rel_regx_tar && in_rel_regz_tar)
                                //     {
                                //         tif_x = false;
                                //         tif_z = false;
                                //     }
                                // }
                                // else if (dist_x == dist_2cells && dist_y == dist_2cells && dist_z == dist_2cells)
                                // {
                                //     if (in_rel_regx_tar && in_rel_regy_tar && in_rel_regz_tar)
                                //     {
                                //         tif_x = false;
                                //         tif_y = false;
                                //         tif_z = false;
                                //     }
                                // }

                                if (dist_x == dist_2cells)
                                {
                                    if (in_rel_regx_tar)
                                    {
                                        tif_x = false;
                                    }
                                }

                                if (dist_y == dist_2cells)
                                {
                                    if (in_rel_regy_tar)
                                    {
                                        tif_y = false;
                                    }
                                }

                                if (dist_z == dist_2cells)
                                {
                                    if (in_rel_regz_tar)
                                    {
                                        tif_z = false;
                                    }
                                }

                                for (const int body_idx_src : adj_cell.bodiesIndices)
                                {
                                    const FBody &body_src = bodies_all_[body_idx_src];

                                    bool x_same_dir_tar = true;
                                    bool y_same_dir_tar = true;
                                    bool z_same_dir_tar = true;

                                    const real rsx_tc = body_src.x[0] - cell.center[0];
                                    const real rsx_sc = body_src.x[0] - adj_cell.center[0];

                                    const real rsy_tc = body_src.x[1] - cell.center[1];
                                    const real rsy_sc = body_src.x[1] - adj_cell.center[1];

                                    const real rsz_tc = body_src.x[2] - cell.center[2];
                                    const real rsz_sc = body_src.x[2] - adj_cell.center[2];

                                    // both particles are in opposite direction of x
                                    if (rtx_tc * rsx_tc < 0)
                                    {
                                        x_same_dir_tar = false;
                                    }

                                    // both particles are in opposite direction of y
                                    if (rty_tc * rsy_tc < 0)
                                    {
                                        y_same_dir_tar = false;
                                    }

                                    // both particles are in opposite direction of z

                                    if (rtz_tc * rsz_tc < 0)
                                    {
                                        z_same_dir_tar = false;
                                    }

                                    bool x_same_dir_src = true;
                                    bool y_same_dir_src = true;
                                    bool z_same_dir_src = true;

                                    if (rsx_sc * rtx_sc < 0)
                                    {
                                        x_same_dir_src = false;
                                    }

                                    if (rsy_sc * rty_sc < 0)
                                    {
                                        y_same_dir_src = false;
                                    }

                                    if (rsz_sc * rtz_sc < 0)
                                    {
                                        z_same_dir_src = false;
                                    }

                                    const RVec ws = w_per_atom[body_idx_src];

                                    bool bx = ws[0] < 1
                                                  ? (fabs(body_src.x[0] - cell.center[0]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_tcell ? 1 : 0)
                                                  : 1;

                                    bool by = ws[1] < 1
                                                  ? (fabs(body_src.x[1] - cell.center[1]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_tcell ? 1 : 0)
                                                  : 1;

                                    bool bz = ws[2] < 1
                                                  ? (fabs(body_src.x[2] - cell.center[2]) + fmm_weights_eval_.getRegAlpha() <= interaction_region_tcell ? 1 : 0)
                                                  : 1;

                                    const bool in_regx_src = fabs(body_src.x[0] - adj_cell.center[0]) + fmm_weights_eval_.getRegAlpha() > adj_cell.radius;
                                    const bool in_regy_src = fabs(body_src.x[1] - adj_cell.center[1]) + fmm_weights_eval_.getRegAlpha() > adj_cell.radius;
                                    const bool in_regz_src = fabs(body_src.x[2] - adj_cell.center[2]) + fmm_weights_eval_.getRegAlpha() > adj_cell.radius;

                                    const bool in_rel_regx_src = in_regx_src && rsx_tc * rsx_sc < 0;

                                    const bool in_rel_regy_src = in_regy_src && rsy_tc * rsy_sc < 0;

                                    const bool in_rel_regz_src = in_regz_src && rsz_tc * rsz_sc < 0;

                                    bool sif_x = true;
                                    bool sif_y = true;
                                    bool sif_z = true;

                                    // if (dist_x == dist_2cells && dist_y != dist_2cells && dist_z != dist_2cells)
                                    // {
                                    //     if (in_rel_regx_src)
                                    //     {
                                    //         sif_x = false;
                                    //     }
                                    // }
                                    // else if (dist_y == dist_2cells && dist_x != dist_2cells && dist_z != dist_2cells)
                                    // {
                                    //     if (in_rel_regy_src)
                                    //     {
                                    //         sif_y = false;
                                    //     }
                                    // }
                                    // else if (dist_z == dist_2cells && dist_x != dist_2cells && dist_y != dist_2cells)
                                    // {
                                    //     if (in_rel_regz_src)
                                    //     {
                                    //         sif_z = false;
                                    //     }
                                    // }
                                    // else if (dist_x == dist_2cells && dist_y == dist_2cells && dist_z != dist_2cells)
                                    // {
                                    //     if (in_rel_regx_src && in_rel_regy_src)
                                    //     {
                                    //         sif_x = false;
                                    //         sif_y = false;
                                    //     }
                                    // }
                                    // else if (dist_y == dist_2cells && dist_z == dist_2cells && dist_x != dist_2cells)
                                    // {
                                    //     if (in_rel_regy_src && in_rel_regz_src)
                                    //     {
                                    //         sif_y = false;
                                    //         sif_z = false;
                                    //     }
                                    // }
                                    // else if (dist_x == dist_2cells && dist_z == dist_2cells && dist_y != dist_2cells)
                                    // {
                                    //     if (in_rel_regx_src && in_rel_regz_src)
                                    //     {
                                    //         sif_x = false;
                                    //         sif_z = false;
                                    //     }
                                    // }
                                    // else if (dist_x == dist_2cells && dist_y == dist_2cells && dist_z == dist_2cells)
                                    // {
                                    //     if (in_rel_regx_src && in_rel_regy_src && in_rel_regz_src)
                                    //     {
                                    //         sif_x = false;
                                    //         sif_y = false;
                                    //         sif_z = false;
                                    //     }
                                    // }

                                    if (dist_x == dist_2cells)
                                    {
                                        if (in_rel_regx_src)
                                        {
                                            sif_x = false;
                                        }
                                    }
                                    if (dist_y == dist_2cells)
                                    {
                                        if (in_rel_regy_src)
                                        {
                                            sif_y = false;
                                        }
                                    }
                                    if (dist_z == dist_2cells)
                                    {
                                        if (in_rel_regz_src)
                                        {
                                            sif_z = false;
                                        }
                                    }

                                    // Check if the source particle is within the x-range of the target cell
                                    bool srcp_in_tarc_x = fabs(body_src.x[0] - cell.center[0]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the source particle is within the y-range of the target cell
                                    bool srcp_in_tarc_y = fabs(body_src.x[1] - cell.center[1]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the source particle is within the z-range of the target cell
                                    bool srcp_in_tarc_z = fabs(body_src.x[2] - cell.center[2]) <= interaction_region_tcell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the target particle is within the x-range of the source cell
                                    bool tarp_in_srcc_x =
                                        fabs(body_tar.x[0] - adj_cell.center[0]) <= interaction_region_scell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the target particle is within the y-range of the source cell
                                    bool tarp_in_srcc_y =
                                        fabs(body_tar.x[1] - adj_cell.center[1]) <= interaction_region_scell + fmm_weights_eval_.getRegAlpha();

                                    // Check if the target particle is within the z-range of the source cell
                                    bool tarp_in_srcc_z =
                                        fabs(body_tar.x[2] - adj_cell.center[2]) <= interaction_region_scell + fmm_weights_eval_.getRegAlpha();

                                    // Determine if the source particle is entirely within the range of the target cell
                                    bool srcp_in_tarc = srcp_in_tarc_x && srcp_in_tarc_y && srcp_in_tarc_z;

                                    // Determine if the target particle is entirely within the range of the source cell
                                    bool tarp_in_srcc = tarp_in_srcc_x && tarp_in_srcc_y && tarp_in_srcc_z;

                                    if (!is_reg_body[body_idx_tar])
                                    {
                                        if (is_reg_body[body_idx_src])
                                        {
                                            if (srcp_in_tarc && (sif_x == false || sif_y == false || sif_z == false))
                                            {
                                                // pair_list[body_idx_tar].push_back(body_idx_src);

                                                // pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                                // pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                // pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                // pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                            }
                                        }
                                    }
                                    else
                                    {

                                        if (!is_reg_body[body_idx_src])
                                        {

                                            // if (tarp_in_srcc && (tif_x == false || tif_y == false || tif_z == false))
                                            // {
                                            //     pair_list[body_idx_tar].push_back(body_idx_src);
                                            //     pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                            //     pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, tif_z});

                                            //     pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                            //     pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                            // }
                                        }
                                        else if (!srcp_in_tarc && !tarp_in_srcc)
                                        {

                                            bool xs_out = sif_x == false && srcp_in_tarc_x;
                                            bool ys_out = sif_y == false && srcp_in_tarc_y;
                                            bool zs_out = sif_z == false && srcp_in_tarc_z;

                                            bool xt_out = tif_x == false && tarp_in_srcc_x;
                                            bool yt_out = tif_y == false && tarp_in_srcc_y;
                                            bool zt_out = tif_z == false && tarp_in_srcc_z;

                                            if (tif_x == false || tif_y == false || tif_z == false)
                                            {

                                                bool sx_c = 1;
                                                bool sy_c = 1;
                                                bool sz_c = 1;

                                                bool tx_c = 1;
                                                bool ty_c = 1;
                                                bool tz_c = 1;

                                                if (srcp_in_tarc_x && sif_x == false)
                                                {
                                                    sx_c = 0;
                                                }

                                                if (srcp_in_tarc_y && sif_y == false)
                                                {
                                                    sy_c = 0;
                                                }

                                                if (srcp_in_tarc_z && sif_z == false)
                                                {
                                                    sz_c = 0;
                                                }

                                                // if (sx_c == false || sy_c == false || sz_c == false)
                                                // {
                                                //     pair_list[body_idx_tar].push_back(body_idx_src);
                                                //     pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                                //     pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, tif_z});

                                                //     pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //     pair_list_sif_within[body_idx_tar].push_back({sx_c, sy_c, sz_c});
                                                // }
                                            }
                                        }
                                        else
                                        {

                                            if (x_same_dir_tar && y_same_dir_tar && z_same_dir_tar)
                                            {
                                                if (srcp_in_tarc && (sif_x == false || sif_y == false || sif_z == false))
                                                {
                                                    pair_list[body_idx_tar].push_back(body_idx_src);

                                                    pair_list_bxyz_tar[body_idx_tar].push_back({1, 1, 1});
                                                    pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                    pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                    pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                                }

                                                if (tarp_in_srcc && (tif_x == false || tif_y == false || tif_z == false))
                                                {
                                                    pair_list[body_idx_tar].push_back(body_idx_src);
                                                    pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                                    pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, tif_z});

                                                    pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                    pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                                }
                                            }

                                            else if (!x_same_dir_tar && y_same_dir_tar && z_same_dir_tar)
                                            {
                                                // if (sif_x == false || sif_y == false || sif_z == false)
                                                // {
                                                //     if (dist_x == 0)
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 1, 1});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({0, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({1, sif_y, sif_z});

                                                //         if (in_regx_src)
                                                //         {

                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 1, 1});
                                                //             pair_list_tif_within[body_idx_tar].push_back({0, 1, 1});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({0, by, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({0, sif_y, sif_z});
                                                //         }
                                                //     }
                                                //     else
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({0, 1, 1});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                                //     }
                                                // }

                                                // if (tif_x == false || tif_y == false || tif_z == false)
                                                // {
                                                //     if (dist_x == 0)
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, bzt});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, tif_y, tif_z});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({1, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         if (in_regx_tar)
                                                //         {

                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, bzt});
                                                //             pair_list_tif_within[body_idx_tar].push_back({0, tif_y, tif_z});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({0, by, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                                //         }
                                                //     }
                                                //     else
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, bzt});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, tif_y, tif_z});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         if (in_rel_regx_src && wst[0] != 1)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, bzt});
                                                //             pair_list_tif_within[body_idx_tar].push_back({0, tif_y, tif_z});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({0, by, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({0, 1, 1});
                                                //         }
                                                //     }
                                                // }
                                            }
                                            else if (x_same_dir_tar && y_same_dir_tar && !z_same_dir_tar)
                                            {
                                                // if (sif_x == false || sif_y == false || sif_z == false)
                                                // {
                                                //     if (dist_z == 0)
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({1, 1, bzt});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, 0});
                                                //         pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, 1});

                                                //         if (in_regz_src)
                                                //         {

                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({1, 1, bzt});
                                                //             pair_list_tif_within[body_idx_tar].push_back({1, 1, 0});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, by, 0});
                                                //             pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, 0});
                                                //         }
                                                //     }
                                                //     else
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({1, 1, 0});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                                //     }
                                                // }

                                                // if (tif_x == false || tif_y == false || tif_z == false)
                                                // {
                                                //     if (dist_z == 0)
                                                //     {

                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, 0});
                                                //         pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, 1});
                                                //         pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         if (in_regz_tar)
                                                //         {

                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, 0});
                                                //             pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, 0});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, by, 0});
                                                //             pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                                //         }
                                                //     }
                                                //     else
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, 0});
                                                //         pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         if (in_rel_regz_src && wst[2] != 1)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, 0});
                                                //             pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, 0});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, by, 0});
                                                //             pair_list_sif_within[body_idx_tar].push_back({1, 1, 0});
                                                //         }
                                                //     }
                                                // }
                                            }

                                            else if (x_same_dir_tar && !y_same_dir_tar && z_same_dir_tar)
                                            {
                                                // if (sif_x == false || sif_y == false || sif_z == false)
                                                // {
                                                //     if (dist_y == 0)
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({1, byt, 1});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({sif_x, 1, sif_z});

                                                //         if (in_regy_src)
                                                //         {

                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({1, byt, 1});
                                                //             pair_list_tif_within[body_idx_tar].push_back({1, 0, 1});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({sif_x, 0, sif_z});
                                                //         }
                                                //     }
                                                //     else
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({1, 0, 1});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                                //     }
                                                // }

                                                // if (tif_x == false || tif_y == false || tif_z == false)
                                                // {
                                                //     if (dist_y == 0)
                                                //     {

                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, bzt});
                                                //         pair_list_tif_within[body_idx_tar].push_back({tif_x, 1, tif_z});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, 1, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         if (in_regy_tar)
                                                //         {

                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, bzt});
                                                //             pair_list_tif_within[body_idx_tar].push_back({tif_x, 0, tif_z});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                                //         }
                                                //     }
                                                //     else
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, bzt});
                                                //         pair_list_tif_within[body_idx_tar].push_back({tif_x, 1, tif_z});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         if (in_rel_regy_src && wst[1] != 1)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, bzt});
                                                //             pair_list_tif_within[body_idx_tar].push_back({tif_x, 0, tif_z});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({1, 0, 1});
                                                //         }
                                                //     }
                                                // }
                                            }

                                            else if (!x_same_dir_tar && y_same_dir_tar && !z_same_dir_tar)
                                            {
                                                // if (sif_x == false || sif_y == false || sif_z == false)
                                                // {
                                                //     if (dist_z == 0 && dist_x == 0)
                                                //     {
                                                //          // skeptical if block, else blocks are fine
                                                //         if (in_regz_src && in_regx_src)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, 0});
                                                //             pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});

                                                //             if (wst[2] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, 0});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 1, 0});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, by, 0});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, 1});
                                                //             }

                                                //             if (wst[0] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({0, 1, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({0, by, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, sif_y, sif_z});
                                                //             }
                                                //         }
                                                //         else
                                                //         {
                                                //             if (in_regx_src || in_regz_src)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 1, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                                //             }
                                                //             else
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 1, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({0, by, 0});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, sif_y, 1});
                                                //             }
                                                //         }
                                                //     }
                                                //     else
                                                //     {

                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({0, 1, 0});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                                //     }
                                                // }

                                                // if (tif_x == false || tif_y == false || tif_z == false)
                                                // {

                                                //     if (dist_z == 0 && dist_x == 0)
                                                //     {
                                                //         if (in_regz_tar && in_regx_tar)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                                //             pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, tif_z});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({0, by, 0});
                                                //             pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //             if (ws[2] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, 0});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, by, 0});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 1, 0});
                                                //             }

                                                //             if (ws[0] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, tif_y, tif_z});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({0, by, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({0, 1, 1});
                                                //             }
                                                //         }
                                                //         else
                                                //         {
                                                //             if (in_regx_tar || in_regz_tar)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, tif_z});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({1, by, 1});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                                //             }
                                                //             else
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, 0});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, tif_y, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({1, by, 1});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                                //             }
                                                //         }
                                                //     }
                                                //     else
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, 0});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, tif_y, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         // not sure if it is needed or not, rethink later , not such interaction yet observed in tests
                                                //         if (in_rel_regz_src && wst[2] != 1 && in_rel_regx_src && wst[0] != 1)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, 0});
                                                //             pair_list_tif_within[body_idx_tar].push_back({0, tif_y, 0});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({0, by, 0});
                                                //             pair_list_sif_within[body_idx_tar].push_back({0, 1, 0});
                                                //         }
                                                //     }
                                                // }
                                            }
                                            else if (!x_same_dir_tar && !y_same_dir_tar && z_same_dir_tar)
                                            {
                                                // std::cout << body_src.x << "--" << body_tar.x << "--" << sif_x << "--" << sif_y << "--" << sif_z << "--"
                                                //           << tif_x << "--" << tif_y << "--" << tif_z << "--" << srcp_in_tarc << "--" << tarp_in_srcc <<
                                                //           std::endl;
                                                // if (srcp_in_tarc && (sif_x == false || sif_y == false || sif_z == false))
                                                // {
                                                //     if (dist_x == 0 && dist_y == 0)
                                                //     {
                                                //         // skeptical if block, else blocks are fine
                                                //         if (in_regx_src && in_regy_src)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({0, 0, bzt});
                                                //             pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});

                                                //             if (wst[0] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({0, 1, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({0, by, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, sif_y, sif_z});
                                                //             }

                                                //             if (wst[1] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 0, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({sif_x, 1, sif_z});
                                                //             }
                                                //         }
                                                //         else
                                                //         {
                                                //             if (in_regy_src || in_regx_src)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, 1});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                                //             }
                                                //             else
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, 1});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({0, 0, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 1, sif_z});
                                                //             }
                                                //         }
                                                //     }
                                                //     else
                                                //     {

                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({0, 0, 1});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                                //     }
                                                // }

                                                // if (tarp_in_srcc && (tif_x == false || tif_y == false || tif_z == false))
                                                // {

                                                //     if (dist_y == 0 && dist_x == 0)
                                                //     {
                                                //         if (in_regy_tar && in_regx_tar)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                                //             pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, tif_z});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({0, 0, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //             if (ws[1] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({tif_x, 1, tif_z});
                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 0, 1});
                                                //             }

                                                //             if (ws[0] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({0, byt, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, tif_y, tif_z});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({0, by, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({0, 1, 1});
                                                //             }
                                                //         }
                                                //         else
                                                //         {
                                                //             if (in_regy_tar || in_regx_tar)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, tif_z});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({1, 1, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                                //             }
                                                //             else
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({0, 0, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 1, tif_z});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({1, 1, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                                //             }
                                                //         }
                                                //     }
                                                //     else
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({0, 0, bzt});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, 1, tif_z});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         // not sure if it is needed or not, rethink later , not such interaction yet observed in tests
                                                //         if (in_rel_regx_src && wst[0] != 1 && in_rel_regy_src && wst[1] != 1)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({0, 0, bzt});
                                                //             pair_list_tif_within[body_idx_tar].push_back({0, 0, tif_z});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({0, 0, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({0, 0, 1});
                                                //         }
                                                //     }
                                                // }
                                            }

                                            // other conditions of one same dir and two opposite dirs are dependent on correctness of this block
                                            else if (x_same_dir_tar && !y_same_dir_tar && !z_same_dir_tar)
                                            {

                                                // if (sif_x == false || sif_y == false || sif_z == false)
                                                // {
                                                //     if (dist_z == 0 && dist_y == 0)
                                                //     {
                                                //          // skeptical if block, else blocks are fine
                                                //         if (in_regz_src && in_regy_src)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, 0});
                                                //             pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //             pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});

                                                //             if (wst[2] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, 0});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 1, 0});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, by, 0});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, 1});
                                                //             }

                                                //             if (wst[1] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 0, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({sif_x, 1, sif_z});
                                                //             }
                                                //         }
                                                //         else
                                                //         {

                                                //             if (in_regy_src || in_regz_src)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({1, byt, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                                //             }
                                                //             else
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({1, byt, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, 0});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({sif_x, 1, 1});
                                                //             }
                                                //         }
                                                //     }
                                                //     else
                                                //     {

                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({1, 0, 0});
                                                //         pair_list_tif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({sif_x, sif_y, sif_z});
                                                //     }
                                                // }

                                                // if (tif_x == false || tif_y == false || tif_z == false)
                                                // {

                                                //     if (dist_z == 0 && dist_y == 0)
                                                //     {

                                                //         if (in_regz_tar && in_regy_tar)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                                //             pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, tif_z});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, 0});
                                                //             pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //             if (ws[2] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, 0});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, by, 0});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 1, 0});
                                                //             }

                                                //             if (ws[1] != 1)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({tif_x, 1, tif_z});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, bz});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 0, 1});
                                                //             }
                                                //         }
                                                //         else
                                                //         {
                                                //             if (in_regy_tar || in_regz_tar)
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, byt, bzt});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({tif_x, tif_y, tif_z});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, 1, 1});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                                //             }
                                                //             else
                                                //             {
                                                //                 pair_list[body_idx_tar].push_back(body_idx_src);
                                                //                 pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, 0});
                                                //                 pair_list_tif_within[body_idx_tar].push_back({tif_x, 1, 1});

                                                //                 pair_list_bxyz_src[body_idx_tar].push_back({bx, 1, 1});
                                                //                 pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});
                                                //             }
                                                //         }
                                                //     }
                                                //     else
                                                //     {
                                                //         pair_list[body_idx_tar].push_back(body_idx_src);
                                                //         pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, 0});
                                                //         pair_list_tif_within[body_idx_tar].push_back({tif_x, 1, 1});

                                                //         pair_list_bxyz_src[body_idx_tar].push_back({bx, by, bz});
                                                //         pair_list_sif_within[body_idx_tar].push_back({1, 1, 1});

                                                //         // not sure if it is needed or not, rethink later , not such interaction yet observed in tests
                                                //         if (in_rel_regz_src && wst[2] != 1 && in_rel_regy_src && wst[1] != 1)
                                                //         {
                                                //             pair_list[body_idx_tar].push_back(body_idx_src);
                                                //             pair_list_bxyz_tar[body_idx_tar].push_back({bxt, 0, 0});
                                                //             pair_list_tif_within[body_idx_tar].push_back({tif_x, 0, 0});

                                                //             pair_list_bxyz_src[body_idx_tar].push_back({bx, 0, 0});
                                                //             pair_list_sif_within[body_idx_tar].push_back({1, 0, 0});
                                                //         }
                                                //     }
                                                // }
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
        for (auto &body_src_idx : pair_list[i])
        {

            gmx::fmm::FBody &asrc = bodies_all_[body_src_idx];

            const BVec bxyz_src = pair_list_bxyz_src[i][ix];
            const BVec bxyz_tar = pair_list_bxyz_tar[i][ix];
            const BVec is_wihin_src = pair_list_sif_within[i][ix];
            const BVec is_within_tar = pair_list_tif_within[i][ix];

            const RVec wsrc_ws = w_per_atom[body_src_idx];

            const real xs = asrc.x[0];
            const real ys = asrc.x[1];
            const real zs = asrc.x[2];

            const real dx = xt - xs;
            const real dy = yt - ys;
            const real dz = zt - zs;

            const real dxs = xs - xt;
            const real dys = ys - yt;
            const real dzs = zs - zt;

            real wsrc_x = bxyz_src[0] == 1 ? 1 : (is_wihin_src[0] == 1 ? wsrc_ws[0] : 1 - wsrc_ws[0]);
            real wsrc_y = bxyz_src[1] == 1 ? 1 : (is_wihin_src[1] == 1 ? wsrc_ws[1] : 1 - wsrc_ws[1]);
            real wsrc_z = bxyz_src[2] == 1 ? 1 : (is_wihin_src[2] == 1 ? wsrc_ws[2] : 1 - wsrc_ws[2]);
            real wsrc = wsrc_x * wsrc_y * wsrc_z;

            real wtar_x = bxyz_tar[0] == 1 ? 1 : (is_within_tar[0] == 1 ? wtar_ws[0] : 1 - wtar_ws[0]);
            real wtar_y = bxyz_tar[1] == 1 ? 1 : (is_within_tar[1] == 1 ? wtar_ws[1] : 1 - wtar_ws[1]);
            real wtar_z = bxyz_tar[2] == 1 ? 1 : (is_within_tar[2] == 1 ? wtar_ws[2] : 1 - wtar_ws[2]);

            const real wtar = wtar_x * wtar_y * wtar_z;

            std::cout << asrc.x << "--" << body_tar.x << "--" << wsrc << "--" << wtar << "" << "\n";
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
