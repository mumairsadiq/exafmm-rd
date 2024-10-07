#include "epkifmm.h"
#include <iostream>

rtfmm::EpkiFMM::EpkiFMM(const Bodies3& bs_, const Argument& args_) : bs(bs_), args(args_)
{
    std::cout<<"EPKIFMM()"<<std::endl;
    assert_exit(bs.size() == args.n, "LaplaceFMM init body size error");
}

rtfmm::Bodies3 rtfmm::EpkiFMM::solve()
{
    std::cout<<"EPKIFMM::solve()"<<std::endl;

    if(args.dipole_correction)
        dipole_correction(bs, args.cycle);
        
    if(args.divide_4pi)
        scale_bodies(bs);

    return bs;
}