#pragma once

#include "type.h"
#include "body.h"
#include "tree.h"

namespace rtfmm
{
    
void make_source_cell(Cell3& cell_src, const Bodies3& bs, int image = 0, real reg = 0.0)
{
    if(image == 0 && reg == 0)
    {
        cell_src.bodies = bs;
    }
}

}