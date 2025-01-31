#ifndef GMX_FMM_H
#define GMX_FMM_H

#include "fmmdirectinteractionstree.h"
#include <functional>

namespace gmx
{
namespace fmm
{

// Structure to hold source flags
struct PairListEntrySrcFlags
{
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

    // Default constructor (sets all flags to true)
    PairListEntrySrcFlags() : bx_src(true), by_src(true), bz_src(true), sx_within(true), sy_within(true), sz_within(true) {}

    // Constructor with parameters
    PairListEntrySrcFlags(bool x, bool y, bool z, bool sx, bool sy, bool sz) : bx_src(x), by_src(y), bz_src(z), sx_within(sx), sy_within(sy), sz_within(sz) {}

    // Setter methods for source flags
    void set_src_flags(bool x, bool y, bool z)
    {
        bx_src = x;
        by_src = y;
        bz_src = z;
    }

    // Setter methods for source within-cell flags
    void set_scw_flags(bool x, bool y, bool z)
    {
        sx_within = x;
        sy_within = y;
        sz_within = z;
    }

    // Define equality operator for unordered_map key comparison
    bool operator==(const PairListEntrySrcFlags &other) const
    {
        return bx_src == other.bx_src && by_src == other.by_src && bz_src == other.bz_src && sx_within == other.sx_within && sy_within == other.sy_within &&
               sz_within == other.sz_within;
    }
};

// Custom hash function for PairListEntrySrcFlags
struct PairListEntrySrcFlagsHash
{
    std::size_t operator()(const PairListEntrySrcFlags &flags) const
    {
        std::size_t h1 = std::hash<bool>{}(flags.bx_src);
        std::size_t h2 = std::hash<bool>{}(flags.by_src);
        std::size_t h3 = std::hash<bool>{}(flags.bz_src);
        std::size_t h4 = std::hash<bool>{}(flags.sx_within);
        std::size_t h5 = std::hash<bool>{}(flags.sy_within);
        std::size_t h6 = std::hash<bool>{}(flags.sz_within);

        // Mix the hashes together using bitwise operations
        return (h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3) ^ (h5 << 4) ^ (h6 << 5)) * 2654435761u;
    }
};

// Structure to hold target flags
struct PairListEntryTargetFlags
{
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

    // Default constructor (sets all flags to true)
    PairListEntryTargetFlags() : bx_tar(true), by_tar(true), bz_tar(true), tx_within(true), ty_within(true), tz_within(true) {}

    // Constructor with parameters
    PairListEntryTargetFlags(bool x, bool y, bool z, bool tx, bool ty, bool tz) : bx_tar(x), by_tar(y), bz_tar(z), tx_within(tx), ty_within(ty), tz_within(tz) {}

    // Setter methods for target flags
    void set_tar_flags(bool x, bool y, bool z)
    {
        bx_tar = x;
        by_tar = y;
        bz_tar = z;
    }

    // Setter methods for target within-cell flags
    void set_trw_flags(bool x, bool y, bool z)
    {
        tx_within = x;
        ty_within = y;
        tz_within = z;
    }
};

// Define the nested hash map structure
using PairListMap = std::unordered_map<PairListEntrySrcFlags, PairListEntryTargetFlags, PairListEntrySrcFlagsHash>;

// Vector of unordered_maps indexed by `body_idx_tar`

class FMMDirectInteractions
{
  public:
    FMMDirectInteractions(const std::vector<RVec> coordinates, const std::vector<real> charges, const RVec box_center, const real box_radius, const size_t max_depth,
                          const real reg_alpha);

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

    std::vector<std::unordered_map<int, PairListMap>> pair_list;
    // weight values for each atom within its original cell
    std::vector<RVec> w_per_atom;

    void compute_weights_();
};

} // namespace fmm
} // namespace gmx

#endif