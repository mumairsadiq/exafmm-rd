#pragma once

#include "type.h"
#include "body.h"
#include "tree.h"

namespace rtfmm
{
    
void make_src_cell(Cell3& cell_src, const Bodies3& bs, vec3r root_x, real root_r)
{
    cell_src.x = root_x;
    cell_src.r = root_r;
    cell_src.bodies = bs;
}

void make_reg_src_cell(Cells3& cell_srcs, const Bodies3& bs, vec3r root_x, real root_r, real cycle, real rega)
{
    for(int z = -1; z <= 1; z++)
    {
        for(int y = -1; y <= 1; y++)
        {
            for(int x = -1; x <= 1; x++)
            {
                vec3r offset = vec3r(x,y,z) * cycle;
                if(offset != vec3r(0,0,0))
                {
                    Cell3 c;
                    c.x = root_x + offset;
                    c.r = root_r;
                    cell_srcs.push_back(c);
                }
            }
        }
    }
}

}