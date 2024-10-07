#pragma once

#include "type.h"
#include "body.h"
#include "argument.h"

namespace rtfmm
{

class EpkiFMM
{
public:
    EpkiFMM(const Bodies3& bs_, const Argument& args_);
    Bodies3 solve();
private:
    Bodies3 bs;
    Argument args;
};

}