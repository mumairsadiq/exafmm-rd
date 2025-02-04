#ifndef GMX_FMM_H
#define GMX_FMM_H

#include "fmmdirectinteractionstree.h"
#include <array>
#include <functional>

constexpr int MAX_ENTRIES = 8;

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

    // Define equality operator for unordered_map key comparison
    bool operator==(const PairListEntrySrcFlags &other) const
    {
        return bx_src == other.bx_src && by_src == other.by_src && bz_src == other.bz_src && sx_within == other.sx_within && sy_within == other.sy_within &&
               sz_within == other.sz_within;
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
};

constexpr int MAX_ENTRIES_IN_FIXED_MAP = 8;

struct FixedPairListMap
{
    struct Entry
    {
        int key = -1;
        PairListEntrySrcFlags entry_src;
        PairListEntryTargetFlags entry_tar;
    };

    std::array<Entry, MAX_ENTRIES_IN_FIXED_MAP> data = {};
    int size = 0;

    static constexpr int encode(const PairListEntrySrcFlags &entry)
    {
        return (entry.bx_src << 5) | (entry.by_src << 4) | (entry.bz_src << 3) | (entry.sx_within << 2) | (entry.sy_within << 1) | entry.sz_within;
    }

    PairListEntryTargetFlags *find(const PairListEntrySrcFlags &entry_src)
    {
        int key = encode(entry_src);
        if (size > 0 && data[0].key == key)
            return &data[0].entry_tar;
        if (size > 1 && data[1].key == key)
            return &data[1].entry_tar;
        if (size > 2 && data[2].key == key)
            return &data[2].entry_tar;
        if (size > 3 && data[3].key == key)
            return &data[3].entry_tar;
        if (size > 4 && data[4].key == key)
            return &data[4].entry_tar;
        if (size > 5 && data[5].key == key)
            return &data[5].entry_tar;
        if (size > 6 && data[6].key == key)
            return &data[6].entry_tar;
        if (size > 7 && data[7].key == key)
            return &data[7].entry_tar;
        return nullptr;
    }

    void insert(const PairListEntrySrcFlags &entry_src, const PairListEntryTargetFlags &entry_tar)
    {
        int key = encode(entry_src);
        PairListEntryTargetFlags *existing = find(entry_src);
        if (existing)
        {
            *existing = entry_tar;
            return;
        }
        if (size < MAX_ENTRIES_IN_FIXED_MAP)
        {
            data[size++] = {key, entry_src, entry_tar};
        }
    }

    struct Iterator
    {
        Entry *ptr;
        bool operator!=(const Iterator &other) const { return ptr != other.ptr; }
        void operator++() { ++ptr; }
        std::pair<const PairListEntrySrcFlags &, PairListEntryTargetFlags &> operator*() { return {ptr->entry_src, ptr->entry_tar}; }
    };

    struct ConstIterator
    {
        const Entry *ptr;
        bool operator!=(const ConstIterator &other) const { return ptr != other.ptr; }
        void operator++() { ++ptr; }
        std::pair<const PairListEntrySrcFlags &, const PairListEntryTargetFlags &> operator*() const { return {ptr->entry_src, ptr->entry_tar}; }
    };

    Iterator begin() { return {data.data()}; }
    Iterator end() { return {data.data() + size}; }

    ConstIterator begin() const { return {data.data()}; }
    ConstIterator end() const { return {data.data() + size}; }
};

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

    PairListEntry() = default;

    PairListEntry(int bd_src_id, bool bx_s, bool by_s, bool bz_s, bool sx_w, bool sy_w, bool sz_w, bool bx_t, bool by_t, bool bz_t, bool tx_w, bool ty_w, bool tz_w)
        : body_idx_src(bd_src_id), bx_src(bx_s), by_src(by_s), bz_src(bz_s), sx_within(sx_w), sy_within(sy_w), sz_within(sz_w), bx_tar(bx_t), by_tar(by_t), bz_tar(bz_t),
          tx_within(tx_w), ty_within(ty_w), tz_within(tz_w)
    {
    }
};

// Vector of unordered_maps indexed by `body_idx_tar`

class FMMDirectInteractions
{
  public:
    FMMDirectInteractions(const std::vector<RVec> coordinates, const std::vector<real> charges, const RVec box_center, const real box_radius, const size_t max_depth,
                          const real reg_alpha);

    bool is_point_within_radius(const RVec &point1, const RVec &point2, double radius);

    // returns forces and potentials pair
    void execute_direct_kernel(real *forces_and_potentials);

    void recompute_weights();

    void rebuild_and_reprocess_tree();

  private:
    FBodies bodies_all_;
    FMMWeightEvaluator fmm_weights_eval_;
    FMMDirectInteractionsTree fmm_direct_interactions_tree_;

    std::vector<std::vector<PairListEntry>> pair_list;

    // weight values for each atom within its original cell
    std::vector<RVec> w_per_atom;

    void compute_weights_();
};

} // namespace fmm
} // namespace gmx

#endif