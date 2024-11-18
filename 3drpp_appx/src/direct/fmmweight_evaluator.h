
#ifndef GMX_FMM_WEIGHTS_EVALUATOR_H
#define GMX_FMM_WEIGHTS_EVALUATOR_H

namespace gmx
{
namespace fmm
{
class FMMWeightEvaluator
{
  public:
    FMMWeightEvaluator(const real reg_alpha) : reg_alpha_(reg_alpha) {}

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

    inline RVec compute_w_xyz(RVec dx, real R)
    {
        RVec res;
        for (int d = 0; d < 3; d++)
        {
            res[d] = compute_w_single(dx[d], R);
        }
        return res;
    }
    inline const real &getRegAlpha() const { return this->reg_alpha_; }

  private:
    const real reg_alpha_;
};
} // namespace fmm
} // namespace gmx

#endif