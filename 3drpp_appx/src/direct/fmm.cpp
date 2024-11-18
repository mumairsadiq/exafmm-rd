#include "fmm.h"

// assuming coordinates are simd packed (like xxx..yyy..zzz.. as per simd width)
gmx::fmm::FMMDirectInteractions::FMMDirectInteractions(
    const std::vector<RVec> coordinates, const std::vector<real> charges,
    const RVec box_center, const real box_radius,
    const size_t max_particles_per_cell, const real reg_alpha)
    : bodies_all_(coordinates, charges),
      fmm_direct_interactions_tree_(bodies_all_, box_center, box_radius,
                                    max_particles_per_cell, reg_alpha)
{
}

std::vector<std::pair<gmx::RVec, gmx::real>>
gmx::fmm::FMMDirectInteractions::execute_direct_kernel()
{

    std::vector<std::pair<gmx::RVec, real>> foreces_and_potentials(
        bodies_all_.size());

    for (size_t i = 0; i < bodies_all_.size(); i++)
    {
        foreces_and_potentials[i] = std::make_pair(RVec(0, 0, 0), 0);
    }

    FMMCells &fmm_cells = fmm_direct_interactions_tree_.get_cells();
    for (size_t i = 0; i < fmm_cells.size(); i++)
    {
        gmx::fmm::FMMCell &cell_tar = fmm_cells[i];
        auto &atomsTarAll = cell_tar.bodiesIndices;

        for (size_t j = 0; j < atomsTarAll.size(); j++)
        {
            const int atom_idx_tar = atomsTarAll[j];

            gmx::fmm::FBody &atar = bodies_all_[atom_idx_tar];
            const real xj = atar.x[0];
            const real yj = atar.x[1];
            const real zj = atar.x[2];
            real pj = 0.0;
            real fxj = 0.0, fyj = 0.0, fzj = 0.0;

            for (const auto &adjacent_cell_idx :
                 fmm_direct_interactions_tree_.get_adjacent_cells(i))
            {

                gmx::fmm::FMMCell &cell_src = fmm_cells[adjacent_cell_idx];
                auto &atomsSrcAll = cell_src.bodiesIndices;

                for (size_t k = 0; k < atomsSrcAll.size(); k++)
                {
                    int atom_idx_src = atomsSrcAll[k];
                    gmx::fmm::FBody &asrc = bodies_all_[atom_idx_src];

                    // Compute distance differences
                    const real dx = xj - asrc.x[0];
                    const real dy = yj - asrc.x[1];
                    const real dz = zj - asrc.x[2];

                    // Compute squared distance
                    real invr = dx * dx + dy * dy + dz * dz;
                    if (invr > 0) // Avoid division by zero
                    {
                        invr =
                            1.0 / std::sqrt(invr); // Compute inverse distance

                        const real qi = asrc.q; // weight source

                        real qinvr = qi * invr;
                        pj += qinvr;
                        qinvr = qinvr * invr * invr;

                        fxj += qinvr * dx;
                        fyj += qinvr * dy;
                        fzj += qinvr * dz;
                    }
                }
            }

            // Apply accumulated forces and potential to target bodies
            foreces_and_potentials[atom_idx_tar].second += pj;
            foreces_and_potentials[atom_idx_tar].first[0] -= fxj;
            foreces_and_potentials[atom_idx_tar].first[1] -= fyj;
            foreces_and_potentials[atom_idx_tar].first[2] -= fzj;
        }
    }

    for (size_t i = 0; i < fmm_cells.size(); i++)
    {
        gmx::fmm::FMMCell &cell_tar = fmm_cells[i];
        auto &atomsTarAll = cell_tar.bodiesIndices;
        for (size_t j = 0; j < atomsTarAll.size(); j++)
        {
            const int atom_idx_tar = atomsTarAll[j];

            gmx::fmm::FBody &atar = bodies_all_[atom_idx_tar];
            const real xj = atar.x[0];
            const real yj = atar.x[1];
            const real zj = atar.x[2];
            real pj = 0.0;
            real fxj = 0.0, fyj = 0.0, fzj = 0.0;

            for (size_t k = 0; k < cell_tar.bodiesIndicesRegularized.size(); k++)
            {
                int atom_idx_src = cell_tar.bodiesIndicesRegularized[k];
                gmx::fmm::FBody &asrc = bodies_all_[atom_idx_src];
                std::cout << asrc.x << "\t" << atar.x << "\t" << cell_tar.weightsBodiesRegularized[k] << "\n";

                // Compute distance differences
                const real dx = xj - asrc.x[0];
                const real dy = yj - asrc.x[1];
                const real dz = zj - asrc.x[2];

                // Compute squared distance
                real invr = dx * dx + dy * dy + dz * dz;


                invr = 1.0 / std::sqrt(invr); // Compute inverse distance

                const real qi =
                    asrc.q * cell_tar.weightsBodiesRegularized[k]; // weight source

                real qinvr = qi * invr;
                pj += qinvr;
                qinvr = qinvr * invr * invr;

                fxj += qinvr * dx;
                fyj += qinvr * dy;
                fzj += qinvr * dz;
            }

            foreces_and_potentials[atom_idx_tar].second += pj;
            foreces_and_potentials[atom_idx_tar].first[0] -= fxj;
            foreces_and_potentials[atom_idx_tar].first[1] -= fyj;
            foreces_and_potentials[atom_idx_tar].first[2] -= fzj;
        }
    }
    
    
    
    return foreces_and_potentials;
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
    fmm_direct_interactions_tree_.recompute_weights();
}

void gmx::fmm::FMMDirectInteractions::rebuild_and_reprocess_tree()
{
    fmm_direct_interactions_tree_.rebuild_and_reprocess_tree();
}
