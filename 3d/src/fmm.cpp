#include "fmm.h"

rtfmm::LaplaceFMM::LaplaceFMM()
{
    printf("nonpara-LaplaceFMM\n");
}

rtfmm::LaplaceFMM::LaplaceFMM(
    Bodies3& bs,
    int P_,
    vec3r x_,
    real r_
) : bodies(bs), P(P_), x(x_), r(r_)
{
    printf("LaplaceFMM\n");

    num_body = bodies.size();

    rtfmm::Cell3 cell;
    cell.idx = 0;
    cell.depth = 0;
    cell.r = r;
    cell.x = {0,0,0};
    cell.crange = {0,0};
    cell.brange = {0,num_body};

    cells.push_back(cell);

    rtfmm::print_bodies(bs, cells[0].brange.number, cells[0].brange.offset);
}

void rtfmm::LaplaceFMM::solve()
{
    printf("solve\n");
}
