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

    // not relevant any more, need to be removed (never used in updated code)
    std::vector<std::vector<real>> pair_list_w_tar;

    // flags indicating whether to use a weight (w=1) or compute weights based on the atom's position in its original cell
    std::vector<std::vector<BVec>> pair_list_bxyz_src; // Source-particle weighting
    std::vector<std::vector<BVec>> pair_list_bxyz_tar; // Target-particle weighting

    // Flags (tif_i) indicates whether the target has weight outside in the direction of adjacent cells beyond its original cell
    std::vector<std::vector<BVec>> pair_list_tif_within; // Target interactions confined within the original cell or not

    // Flags indicating whether the source has weight outside in the direction of adjacent cells beyond its original cell
    std::vector<std::vector<BVec>> pair_list_sif_within; // Source interactions confined within the original cell or not

    // weight values for each atom within its original cell
    std::vector<RVec> w_per_atom;

    void compute_weights_();
};

} // namespace fmm
} // namespace gmx

#endif