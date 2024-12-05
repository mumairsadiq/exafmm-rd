

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

    virtual ~IFMMTree() = default;

    virtual gmx::RVec get_box_center() const = 0;

    virtual gmx::real get_box_radius() const = 0;
};

class FMMTree : public IFMMTree
{
  public:
    FMMTree(const FBodies &bodies, const RVec box_center, const real box_radius, const size_t max_depth);

    // Overriding methods from IFMMTree
    FMMCells &get_cells() override;

    gmx::RVec get_box_center() const override;

    gmx::real get_box_radius() const override;

  private:
    RVec get_child_cell_x_(RVec x_par, real r_par, int octant);

  protected:
    // Data members
    const FBodies &bodies_;

    FMMCells fmm_cells_;

    const RVec box_center_;

    const real box_radius_;

    const size_t max_depth_;

    // methods

    void build_tree_uniform();

    FMMCells get_leaves();
};

} // namespace fmm

} // namespace gmx

#endif