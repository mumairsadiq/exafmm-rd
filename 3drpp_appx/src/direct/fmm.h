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
    FMMDirectInteractions(const std::vector<RVec> coordinates,
                          const std::vector<real> charges,
                          const RVec box_center, const real box_radius,
                          const size_t max_particles_per_cell,
                          const real reg_alpha);

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

    std::vector<std::vector<int>> atoms_interactions_list;
    std::vector<std::vector<real>> atoms_interactions_weights;

    void compute_weights_();
};

} // namespace fmm
} // namespace gmx

#endif