#pragma once

#include "type.h"
#include "body.h"
#include "argument.h"
#include "fmm.h"

namespace rtfmm
{

class EpkiFMM : public LaplaceFMM
{
public:
    EpkiFMM(const Bodies3& bs_, const Argument& args_);
    Bodies3 solve();

private:
    void calculate_root_check();


    void calculate_root_check2();


    void calculate_root_check3();


    void calculate_root_check4();

    void calculate_root_check5();

    void calculate_root_check6();
};

}