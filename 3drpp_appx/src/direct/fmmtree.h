

#ifndef GMX_FMM_TREE_H
#define GMX_FMM_TREE_H

#include "fmmtypes.h"

namespace gmx
{
namespace fmm
{
class IFMMTree
{
  public:
    virtual FMMCells &get_cells() = 0;
    virtual void find_min_max_level() = 0;
    virtual ~IFMMTree() = default;
};

class FMMTree : public IFMMTree
{
  public:
    FMMTree(const FBodies &bodies, const RVec box_center, const real box_radius,
            const size_t max_cell_param, const bool is_tree_uniform = true);

    // Overriding methods from IFMMTree
    FMMCells &get_cells() override;
    FMMCells getLeaves();

    gmx::RVec get_box_center() const;

    gmx::real get_box_radius() const;

    void find_min_max_level() override;

  private:
    RVec get_child_cell_x_(RVec x_par, real r_par, int octant);

  protected:
    // Data members
    const FBodies &bodies_;
    FMMCells fmm_cells_;

    int min_depth_;
    int max_depth_;

    const RVec box_center_;
    const real box_radius_;

    // Determines the cell subdivision parameter: it represents
    // - If using a uniform tree:  then maximum depth of the tree.
    // - If using an adaptive tree: then maximum number of particles per cell
    const size_t cell_limit_param_;

    // methods
    void build_tree_non_uniform();
    void build_tree_uniform();
};

} // namespace fmm

} // namespace gmx

#endif