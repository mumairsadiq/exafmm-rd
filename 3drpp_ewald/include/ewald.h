#pragma once
#include "type.h"
#include "body.h"
#include "argument.h"
#include "tree.h"
#include <complex>

namespace rtfmm
{
struct IndexAndOffset
{
    IndexAndOffset(int index_, vec3r offset_) : index(index_), offset(offset_){}
    int index;
    vec3r offset;
};

class RealMap
{
public:
    void append(int tar, int src, const vec3r& offset)
    {
        real_map[tar].push_back(IndexAndOffset(src, offset));
    }
    std::vector<IndexAndOffset> get_list(int tar)
    {
        return real_map[tar];
    }
    std::map<int, std::vector<IndexAndOffset>> real_map;
};

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
    EwaldSolver(const Bodies3& tar_bs_, const Argument& args_);
    EwaldSolver(const Bodies3& tar_bs_, const Bodies3& src_bs_, const Argument& args_);
    Bodies3 solve();
private:
    void real_part(int this_cell_idx, int that_cell_idx, RealMap& real_map);

    void fourier_part();

    void self_correction();

    void child(int this_cell_idx, int that_cell_idx, RealMap& real_map);

    void real_p2p(RealMap& real_map);

    void real_p2p_kernel(int this_cell_idx, int that_cell_idx, vec3r offset);

    void DFT(std::vector<Wave>& ws, std::vector<Body3>& bs);

    void IDFT(std::vector<Wave>& ws, std::vector<Body3>& bs);

    Argument args;
    Bodies3 tar_bs;
    Bodies3 src_bs;

    int ksize;
    real scale;
    real alpha;
    real cutoff;
    std::vector<Wave> waves;
    std::vector<Cell3> tar_cells;
    std::vector<Cell3> src_cells;
};
}