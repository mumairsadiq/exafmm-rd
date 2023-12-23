#pragma once
#include "type.h"
#include "body.h"
#include "tree.h"
#include "argument.h"
#include "kernel.h"

namespace rtfmm
{

class LaplaceFMM
{
public:
    LaplaceFMM(const Bodies3& bs_, const Argument& args_);

    Bodies3 solve();

private:
    Argument args;
    Bodies3 bs;
    Cells3 cells;
    LaplaceKernel kernel;
};

}