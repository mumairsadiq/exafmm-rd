#ifndef GMX_FMM_WEIGHTS_EVALUATOR_H
#define GMX_FMM_WEIGHTS_EVALUATOR_H

namespace gmx
{
namespace fmm
{
/**
 * @class FMMWeightEvaluator
 * This class calculates weights for the Fast Multipole Method (FMM).
 * It computes the regularization weights based on the distance between
 * particles and their interaction region.
 */
class FMMWeightEvaluator
{

  private:
    inline real compute_reg_w(real x) { return 0.25 * (2 + 3 * x - x * x * x); }

    inline real compute_w_single(real dx, real R)
    {
        real r = std::abs(dx) - R + reg_alpha_;
        real x;

        if (r <= 0)
            x = 1;
        else if (r > 2 * reg_alpha_)
            x = -1;
        else
            x = 1 - r / reg_alpha_;

        return compute_reg_w(x);
    }

    inline real compute_w(RVec dx, real R)
    {
        real w = 1;
        for (int d = 0; d < 3; d++)
        {
            w *= compute_w_single(dx[d], R);
        }
        return w;
    }

  public:
    /**
     * Constructor
     * @param reg_alpha A regularization parameter that determines the
     * smoothness of the weight function.
     */
    FMMWeightEvaluator(const RVec box_center, const real box_radius,
                       const real reg_alpha)
        : box_center_(box_center), box_radius_(box_radius),
          reg_alpha_(reg_alpha)
    {
    }

    /**
     * Compute weights separately for each dimension in 3D space.
     * @param dx The distance vector (x, y, z).
     * @param R The radius of the interaction region.
     * @return A vector of weights, one for each dimension (x, y, z).
     */
    inline RVec compute_w_xyz(RVec dx, real R)
    {
        RVec ws; // Result vector to store weights for each dimension.
        for (int d = 0; d < 3; d++)
        {
            ws[d] =
                compute_w_single(dx[d], R); // Compute weight for dimension `d`.
        }
        return ws; // Return the vector of weights.
    }

    inline real compute_weight(RVec body_x, RVec center, real radius,
                               bool is_periodic = false)
    {

        const RVec dx = body_x - center;
        RVec ws = compute_w_xyz(dx, radius);
        if (is_periodic == false)
        {
            RVec dx_simcenter_inter = body_x - box_center_;
            dx_simcenter_inter[0] = fabs(dx_simcenter_inter[0]);
            dx_simcenter_inter[1] = fabs(dx_simcenter_inter[1]);
            dx_simcenter_inter[2] = fabs(dx_simcenter_inter[2]);
            const RVec dx_simcenter(dx_simcenter_inter[0] + reg_alpha_,
                                    dx_simcenter_inter[1] + reg_alpha_,
                                    dx_simcenter_inter[2] + reg_alpha_);

            for (int d = 0; d <= 2; d++)
            {
                if (dx_simcenter[d] >= box_radius_)
                {
                    ws[d] = 1;
                }
            }
        }
        return ws[0] * ws[1] * ws[2]; // Return the weight
    }

    inline const real &getRegAlpha() const { return this->reg_alpha_; }

  private:
    const RVec box_center_;
    const real box_radius_;
    const real reg_alpha_;
};
} // namespace fmm
} // namespace gmx

#endif
