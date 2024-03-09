#pragma once
#include "type.h"
#include "body.h"
#include "argument.h"
#include "tree.h"
#include <complex>

namespace rtfmm
{
class EwaldSolver
{
private:
    using complexr = std::complex<real>;
    struct Wave 
    {
        vec3r K;
        complexr val;
    };
public:
    EwaldSolver(const Bodies3& bs_, const Argument& args_);
    Bodies3 solve();
private:
    void real_part(int this_cell_idx, int that_cell_idx);

    void fourier_part();

    void self_correction();

    void child(int this_cell_idx, int that_cell_idx);

    void real_p2p(int this_cell_idx, int that_cell_idx, vec3r offset);

    void DFT(std::vector<Wave>& ws, std::vector<Body3>& bs);

    void IDFT(std::vector<Wave>& ws, std::vector<Body3>& bs);

    Argument args;
    Bodies3 bs;

    int ksize;
    real scale;
    real alpha;
    real cutoff;
    std::vector<Wave> waves;
    std::vector<Cell3> cells;
};
}