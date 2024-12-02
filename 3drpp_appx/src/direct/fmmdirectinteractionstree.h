

#ifndef GMX_FMM_DIRECT_TREE_H
#define GMX_FMM_DIRECT_TREE_H

#include "fmmtree.h"
#include "fmmweight_evaluator.h"

namespace gmx
{
namespace fmm
{

class FMMDirectInteractionsTree : public FMMTree
{
  public:
    FMMDirectInteractionsTree(const FBodies &bodies, const RVec box_center,
                              const real box_radius,
                              const size_t cell_limit_param,
                              const bool is_tree_uniform = true);

    // New methods specific to direct interactions
    const FPIndices &get_adjacent_cells(size_t i) const;

    const std::unordered_set<int> &get_adjacent_cells_set(size_t i) const;

    void rebuild_and_reprocess_tree();

  private:
    void process_tree_();

    void find_all_adjacent_cells_();

    int check_adjacent_parent_(int ca_idx, int cb_idx);

    // Data members

    std::vector<FPIndices> adjacent_cells_info_;
    std::vector<std::unordered_set<int>> adjacent_cells_info_set_;
    const int num_neighbours;
};

} // namespace fmm
} // namespace gmx

#endif