#include "type.h"
#include "body.h"
#include "tree.h"
#include "fmm.h"
#include "argument.h"
#include "ewald.h"
#include "direct.h"
#include <omp.h>
#include "timer.h"
#include "epkifmm.h"

int main(int argc, char* argv[])
{
    rtfmm::title(argv[0]);

    rtfmm::Argument args(argc, argv);
    args.show(args.res_filepath);

    return 0;
}