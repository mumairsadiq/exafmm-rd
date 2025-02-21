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
struct WeightFlags
{
    // Boolean flag to determine the weights for each particle pair:
    // - If bx, by, or bz is true, use a weight of 1.
    // - Otherwise, the decision is based on x_in, y_in, or z_in:
    //   - If x_in, y_in, or z_in is true, use weight w.
    //   - If not, use weight (1 - w).
    bool bx, by, bz;

    // Boolean flag indicating whether the weight for a particle pair
    // should be taken within the source cell or outside it.
    // Only valid if bx, by, or bz is false.
    bool x_in, y_in, z_in;

    // Default constructor (sets all flags to true)
    WeightFlags() : bx(true), by(true), bz(true), x_in(true), y_in(true), z_in(true) {}

    // Constructor with parameters
    WeightFlags(bool bx, bool by, bool bz, bool x_in, bool y_in, bool z_in) : bx(bx), by(by), bz(bz), x_in(x_in), y_in(y_in), z_in(z_in) {}

    // Define equality operator for unordered_map key comparison
    bool operator==(const WeightFlags &other) const { return bx == other.bx && by == other.by && bz == other.bz && x_in == other.x_in && y_in == other.y_in && z_in == other.z_in; }
};

constexpr int MAX_ENTRIES_IN_FIXED_MAP = 4;

struct FixedPairListMap
{
    struct Entry
    {
        int key = -1;
        WeightFlags wsrc_flgs;
        WeightFlags wtar_flgs;
    };

    std::array<Entry, MAX_ENTRIES_IN_FIXED_MAP> data = {};
    size_t size = 0;

    static constexpr int encode(const WeightFlags &entry) { return (entry.bx << 5) | (entry.by << 4) | (entry.bz << 3) | (entry.x_in << 2) | (entry.y_in << 1) | entry.z_in; }

    const size_t get_size() const { return size; }

    // return target flags for given source flags
    WeightFlags *find(const WeightFlags &wsrc_flgs)
    {
        int key = encode(wsrc_flgs);
        if (size > 0 && data[0].key == key)
            return &data[0].wtar_flgs;
        if (size > 1 && data[1].key == key)
            return &data[1].wtar_flgs;
        if (size > 2 && data[2].key == key)
            return &data[2].wtar_flgs;
        if (size > 3 && data[3].key == key)
            return &data[3].wtar_flgs;
        if (size > 4 && data[4].key == key)
            return &data[4].wtar_flgs;
        return nullptr;
    }

    void insert(const WeightFlags &entry_src, const WeightFlags &entry_tar)
    {
        int key = encode(entry_src);
        WeightFlags *existing = find(entry_src);
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
        std::pair<const WeightFlags &, WeightFlags &> operator*() { return {ptr->wsrc_flgs, ptr->wtar_flgs}; }
    };

    struct ConstIterator
    {
        const Entry *ptr;
        bool operator!=(const ConstIterator &other) const { return ptr != other.ptr; }
        void operator++() { ++ptr; }
        std::pair<const WeightFlags &, const WeightFlags &> operator*() const { return {ptr->wsrc_flgs, ptr->wtar_flgs}; }
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

class FMMDirectInteractions
{
  public:
    FMMDirectInteractions(const std::vector<RVec> coordinates, const std::vector<real> charges, const RVec box_center, const real box_radius, const size_t max_depth,
                          const real reg_alpha);

    bool is_point_within_radius(const RVec &point1, const RVec &point2, double radius);

    // returns forces and potentials pair
    void execute_direct_kernel(real *forces_and_potentials);

    void compute_group_interactions_(const std::vector<FBody> &sbodies, const std::vector<size_t> &group_counts, const std::vector<size_t> &group_prefix_sum,
                                     real *forces_potentials);

    void recompute_weights();

    void rebuild_and_reprocess_tree();

    u_int32_t get_num_groups();

  private:
    FBodies bodies_all_;
    FMMWeightEvaluator fmm_weights_eval_;
    FMMDirectInteractionsTree fmm_direct_interactions_tree_;

    std::vector<std::vector<PairListEntry>> pair_list;

    void compute_weights_();

    std::vector<int> group_bodies;

    std::unordered_map<std::string, uint32_t> group_map;
    uint32_t next_group_id = 0;
    u_int32_t get_group_id(int ocell_idx, int a_cells_idxs[], size_t valid_size);
};

} // namespace fmm
} // namespace gmx

#endif