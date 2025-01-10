#ifndef GMX_FMM_H
#define GMX_FMM_H

#include "fmmdirectinteractionstree.h"

namespace gmx
{
namespace fmm
{

// Structure to hold flags for source and target weights for each particle pair
struct PairListEntry
{

    // Source body ID
    int body_idx_src;

    // Boolean flag to determine the source weights for each particle pair:
    // - If bx_src, by_src, or bz_src is true, use a weight of 1.
    // - Otherwise, the decision is based on sx_within, sy_within, or sz_within:
    //   - If sx_within, sy_within, or sz_within is true, use weight w.
    //   - If not, use weight (1 - w).
    bool bx_src, by_src, bz_src;

    // Boolean flag indicating whether the weight for a particle pair
    // should be taken within the source cell or outside it.
    // Only valid if bx_src, by_src, or bz_src is false.
    bool sx_within, sy_within, sz_within;

    // Boolean flag to determine the target weights for each particle pair:
    // - If bx_tar, by_tar, or bz_tar is true, use a weight of 1.
    // - Otherwise, the decision is based on tx_within, ty_within, or tz_within:
    //   - If tx_within, ty_within, or tz_within is true, use weight w.
    //   - If not, use weight (1 - w).
    bool bx_tar, by_tar, bz_tar;

    // Boolean flag indicating whether the weight for a particle pair
    // should be taken within the target cell or outside it.
    // Only valid if bx_tar, by_tar, or bz_tar is false.
    bool tx_within, ty_within, tz_within;

    // Default constructor (sets all flags to true and initializes body_idx_src)
    PairListEntry(int bd_src_id)
        : body_idx_src(bd_src_id), bx_src(true), by_src(true), bz_src(true), bx_tar(true), by_tar(true), bz_tar(true), sx_within(true), sy_within(true), sz_within(true), tx_within(true), ty_within(true),
          tz_within(true)
    {
    }

    // Setter methods for source flags
    void set_src_flags(bool x, bool y, bool z)
    {
        bx_src = x;
        by_src = y;
        bz_src = z;
    }

    // Setter methods for target flags
    void set_tar_flags(bool x, bool y, bool z)
    {
        bx_tar = x;
        by_tar = y;
        bz_tar = z;
    }

    // Setter methods for source within-cell flags
    void set_scw_flags(bool x, bool y, bool z)
    {
        sx_within = x;
        sy_within = y;
        sz_within = z;
    }

    // Setter methods for target within-cell flags
    void set_trw_flags(bool x, bool y, bool z)
    {
        tx_within = x;
        ty_within = y;
        tz_within = z;
    }
};

// Class to manage FMMDirectInteractions
class FMMDirectInteractions
{
  public:
    FMMDirectInteractions(const std::vector<RVec> coordinates, const std::vector<real> charges, const RVec box_center, const real box_radius, const size_t max_depth, const real reg_alpha);

    bool is_point_within_radius(const RVec &point1, const RVec &point2, double radius);

    // Returns forces and potentials pair
    std::vector<std::pair<RVec, real>> execute_direct_kernel();

    // Returns forces and potentials pair
    std::vector<std::pair<RVec, real>> execute_direct_kernel_simd();

    void recompute_weights();

    void rebuild_and_reprocess_tree();

  private:
    FBodies bodies_all_;
    FMMWeightEvaluator fmm_weights_eval_;
    FMMDirectInteractionsTree fmm_direct_interactions_tree_;

    // List of particle pairs with their associated weight flags
    std::vector<std::vector<PairListEntry>> pair_list;

    // Weight values for each atom within its original cell
    std::vector<RVec> w_per_atom;

    void compute_weights_();
};

} // namespace fmm
} // namespace gmx

#endif