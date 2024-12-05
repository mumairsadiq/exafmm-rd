

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
    FMMDirectInteractionsTree(const FBodies &bodies, const RVec box_center, const real box_radius, const size_t max_depth_);

    int get_neighbour_idx(RVec neighbor_center);

    void rebuild_and_reprocess_tree();

  private:
    std::unordered_map<long, int> cells_map;

    void process_tree_();
};

} // namespace fmm
} // namespace gmx

#endif