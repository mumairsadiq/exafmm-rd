

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
            const size_t max_particles_per_cell);

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

    const size_t max_particles_per_cell_;

    // methods
    void build_tree();
};

} // namespace fmm

} // namespace gmx

#endif