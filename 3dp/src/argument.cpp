#include "argument.h"

int rtfmm::verbose;

rtfmm::Argument::Argument(int argc, char* argv[])
{
    
    cmdline::parser cmd;
    cmd.add<int>("P", 'P', "point number per edge of surface point box", false, 4, cmdline::range(2, 32));
    cmd.add<int>("nbody", 'n', "number of bodies", false, 1000, cmdline::range(1, (int)std::pow(2,16)));
    cmd.add<int>("ncrit", 'm', "minimum number of bodies per leaf box", false, 16);
    cmd.add<int>("timing", 't', "if measure execution time of computational steps", false, 1);
    cmd.add<int>("images", 'i', "periodic images depth", false, 0, cmdline::range(0, 20));
    cmd.add<int>("verbose", 'v', "is verbose", false, 1);
    cmd.add<real>("cycle", 'c', "cycle of images(or box_r)", false, 2 * M_PI);
    cmd.parse_check(argc, argv);
    P = cmd.get<int>("P");
    n = cmd.get<int>("nbody");
    ncrit = cmd.get<int>("ncrit");
    timing = cmd.get<int>("timing");
    images = cmd.get<int>("images");
    verbose = cmd.get<int>("verbose");
    cycle = cmd.get<real>("cycle");

    x = vec3r(0,0,0);
    r = cycle / 2;
}

void rtfmm::Argument::show()
{
    printf("[P=%d, n=%d, r=%.4f, x=(%.3f,%.3f,%.3f), ncrit=%d, timing=%d, images=%d, cycle=%.4f, verbose=%d]\n", P, n, r, x[0], x[1], x[2], ncrit, timing, images, cycle, verbose);
}