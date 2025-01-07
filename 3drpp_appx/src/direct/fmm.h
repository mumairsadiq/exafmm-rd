#ifndef GMX_FMM_H
#define GMX_FMM_H

#include "fmmdirectinteractionstree.h"

namespace gmx
{
namespace fmm
{

class FMMDirectInteractions
{
  public:
    FMMDirectInteractions(const std::vector<RVec> coordinates, const std::vector<real> charges, const RVec box_center, const real box_radius, const size_t max_depth, const real reg_alpha);

    bool is_point_within_radius(const RVec &point1, const RVec &point2, double radius);

    // returns forces and potentials pair
    std::vector<std::pair<RVec, real>> execute_direct_kernel();

    // returns forces and potentials pair
    std::vector<std::pair<RVec, real>> execute_direct_kernel_simd();

    void recompute_weights();

    void rebuild_and_reprocess_tree();

  private:
    FBodies bodies_all_;
    FMMWeightEvaluator fmm_weights_eval_;
    FMMDirectInteractionsTree fmm_direct_interactions_tree_;

    std::vector<std::vector<int>> pair_list;

    // Boolean flag to determine the source weights for each particle pair:
    // - If bxy_src[i][d] == 1, use a weight of 1.
    // - Otherwise, the decision is based on sif[i][d]:
    //   - If sif[i][d] == 1, use weight w.
    //   - If sif[i][d] == 0, use weight (1 - w).
    std::vector<std::vector<BVec>> pair_list_bxyz_src;

    // Boolean flag indicating whether the weight for a particle pair 
    // should be taken within the source cell or outside it.
    // only valid if bxy_src[i][d] != 1
    std::vector<std::vector<BVec>> pair_list_sif_within; 

    // Boolean flag indicating whether the target weights for a particle pair
    // - If bxy_tar[i][d] == 1, use a weight of 1.
    // - Otherwise, the decision is based on tif[i][d]:
    //   - If tif[i][d] == 1, use weight w.
    //   - If tif[i][d] == 0, use weight (1 - w).
    std::vector<std::vector<BVec>> pair_list_tif_within;

    // Boolean flag indicating whether the weight for a particle pair
    // should be taken within the target cell or outside it.
    // only valid if bxy_tar[i][d] != 1
    std::vector<std::vector<BVec>> pair_list_bxyz_tar;
    
    // weight values for each atom within its original cell
    std::vector<RVec> w_per_atom;

    void compute_weights_();
};

} // namespace fmm
} // namespace gmx

#endif