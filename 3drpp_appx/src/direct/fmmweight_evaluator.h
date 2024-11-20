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
  public:
    /**
     * Constructor
     * @param reg_alpha A regularization parameter that determines the
     * smoothness of the weight function.
     */
    FMMWeightEvaluator(const real reg_alpha) : reg_alpha_(reg_alpha) {}

    /**
     * Compute the regularization weight function for a given value.
     * This function defines how the weight decays as the distance increases.
     * @param x Input value for the weight function.
     * @return Computed weight based on the regularization formula.
     */
    inline real compute_reg_w(real x) { return 0.25 * (2 + 3 * x - x * x * x); }

    /**
     * Compute the weight for a single dimension (1D).
     * @param dx The distance in one dimension.
     * @param R The radius of the interaction region.
     * @return Weight for the given dimension based on the distance.
     */
    inline real compute_w_single(real dx, real R)
    {
        // Calculate the relative distance from the boundary of the interaction
        // region.
        real r = std::abs(dx) - R + reg_alpha_;
        real x;

        // Determine the value of `x` based on the distance `r`.
        if (r <= 0) // Inside the interaction region.
            x = 1;
        else if (r > 2 * reg_alpha_) // Far outside the interaction region.
            x = -1;
        else // Smooth transition in the regularization zone.
            x = 1 - r / reg_alpha_;

        // Compute the regularization weight using `x`.
        return compute_reg_w(x);
    }

    /**
     * Compute the overall weight for a 3D vector.
     * This combines weights from all three dimensions.
     * @param dx The distance vector (x, y, z).
     * @param R The radius of the interaction region.
     * @return The overall weight based on the 3D distance.
     */
    inline real compute_w(RVec dx, real R)
    {
        real w = 1;                 // Start with a weight of 1.
        for (int d = 0; d < 3; d++) // Multiply the weights from each dimension.
        {
            w *= compute_w_single(dx[d], R);
        }
        return w; // Return the final weight.
    }

    /**
     * Compute weights separately for each dimension in 3D space.
     * @param dx The distance vector (x, y, z).
     * @param R The radius of the interaction region.
     * @return A vector of weights, one for each dimension (x, y, z).
     */
    inline RVec compute_w_xyz(RVec dx, real R)
    {
        RVec res; // Result vector to store weights for each dimension.
        for (int d = 0; d < 3; d++)
        {
            res[d] =
                compute_w_single(dx[d], R); // Compute weight for dimension `d`.
        }
        return res; // Return the vector of weights.
    }

    /**
     * Get the regularization parameter (reg_alpha).
     * @return The value of reg_alpha used in the calculations.
     */
    inline const real &getRegAlpha() const { return this->reg_alpha_; }

  private:
    const real
        reg_alpha_; // Regularization parameter for smooth weight transitions.
};
} // namespace fmm
} // namespace gmx

#endif
