#pragma once

#include "type.h"
#include "body.h"
#include "tree.h"
#include "fmm.h"
#include "reg.h"

namespace rtfmm
{

void make_tar_cell(Cell3& cell_tar, const Bodies3& bs, vec3r root_x, real root_r)
{
    cell_tar.x = root_x;
    cell_tar.r = root_r;
    cell_tar.bodies = bs;
}

void make_reg_tar_cell(Cell3& cell_tar, const Bodies3& bs, vec3r root_x, real root_r, real cycle, real rega)
{
    cell_tar.x = root_x;
    cell_tar.r = root_r;
    
    for(int i = 0; i < bs.size(); i++)
    {
        const Body3& b = bs[i];
        real w = reg::get_w_xyz(b.x - cell_tar.x, cell_tar.r, rega).mul();
        cell_tar.bodies.push_back(b);
        cell_tar.weights.push_back(w);
        if(w != 1)
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
                            Body3 bb = b;
                            bb.x += offset;
                            real w = reg::get_w_xyz(bb.x - cell_tar.x, cell_tar.r, rega).mul();
                            if(w > 0)
                            {
                                cell_tar.bodies.push_back(bb);
                                cell_tar.weights.push_back(w);
                            }
                        }
                    }
                }
            }
        }
    }
}
    
void make_src_cell(Cell3& cell_src, const Bodies3& bs, vec3r root_x, real root_r)
{
    cell_src.x = root_x;
    cell_src.r = root_r;
    cell_src.bodies = bs;
}

void make_reg_src_cell(Cells3& cell_srcs, const Bodies3& bs, vec3r root_x, real root_r, real cycle, real rega, int image)
{
    if(image == 0)
    {
        cell_srcs.push_back(Cell3());
        make_reg_tar_cell(cell_srcs[0], bs, root_x, root_r, cycle, rega); // when image is 0, src and tar are same, see #1
    }
    else
    {
        /*Cell3 root;
        root.x = root_x;
        root.r = root_r;
        root.bodies = bs;
        cell_srcs.push_back(root);
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

                        // TODO : add reg body
                    }
                }
            }
        }*/
    }
}

void fusion_bodies(Cell3& cell_tar, Bodies3& res)
{
    int bnum = cell_tar.bodies.size();
    for(int i = 0; i < bnum; i++)
    {
        Body3& b = cell_tar.bodies[i];
        res[b.idx].p += b.p;
        res[b.idx].f += b.f;
    }
}

}