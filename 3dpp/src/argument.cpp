#include "argument.h"

#define CONTAINS(str, sub_str) (str.find(sub_str) != std::string::npos)

int rtfmm::verbose;

rtfmm::Argument::Argument(int argc, char* argv[])
{
    
    cmdline::parser cmd;
    cmd.add<int>("P", 'P', "point number per edge of surface point box", false, 4, cmdline::range(2, 32));
    cmd.add<int>("nbody", 'n', "number of bodies", false, 1000);
    cmd.add<int>("ncrit", 'm', "minimum number of bodies per leaf box", false, 16);
    cmd.add<int>("timing", 't', "if measure execution time of computational steps", false, 1);
    cmd.add<int>("images", 'i', "periodic images depth", false, 0, cmdline::range(0, 20));
    cmd.add<int>("verbose", 'v', "is verbose", false, 1);
    cmd.add<real>("cycle", 'c', "cycle of images(or box_r)", false, 2 * M_PI);
    cmd.add<int>("num_compare", 0, "compare number of target body(also direct calculation number)", false, -1);
    cmd.add<int>("ewald_ksize", 0, "ksize of ewald DFT", false, 11);
    cmd.add<int>("th_num", 0, "number of omp threads", false, 4);
    cmd.add<int>("seed", 0, "random seed", false, 5);
    cmd.add<std::string>("algorithm", 'a', "algorithms", false, "fde");
    cmd.add<int>("check_tree", 0, "if check tree", false, 1);
    cmd.add<int>("use_fft", 0, "if use fft in M2L", false, 1);
    cmd.add<int>("use_precompute", 0, "if use precomputed kernel", false, 1);
    cmd.parse_check(argc, argv);
    P = cmd.get<int>("P");
    n = cmd.get<int>("nbody");
    ncrit = cmd.get<int>("ncrit");
    timing = cmd.get<int>("timing");
    images = cmd.get<int>("images");
    verbose = cmd.get<int>("verbose");
    cycle = cmd.get<real>("cycle");
    use_fft = cmd.get<int>("use_fft");
    use_precompute = cmd.get<int>("use_precompute");

    num_compare = cmd.get<int>("num_compare");
    if(num_compare == -1) num_compare = n;
    ewald_ksize = cmd.get<int>("ewald_ksize");
    th_num = cmd.get<int>("th_num");
    seed = cmd.get<int>("seed");

    std::string algo = cmd.get<std::string>("algorithm");
    enable_fmm = CONTAINS(algo, "f") ? 1 : 0;
    enable_direct = CONTAINS(algo, "d") ? 1 : 0;
    enable_ewald = CONTAINS(algo, "e") ? 1 : 0;
    check_tree = cmd.get<int>("check_tree");

    x = vec3r(0,0,0);
    r = cycle / 2;
}

void rtfmm::Argument::show()
{
    //printf("[P=%d, n=%d, ncrit=%d, images=%d, cycle=%.4f, ewald_ksize=%d]\n", P, n, ncrit, images, cycle, ewald_ksize);
    //printf("[(f,d,e)=(%d,%d,%d), num_compare=%d , th_num=%d, seed=%d, check_tree=%d]\n\n", enable_fmm, enable_direct, enable_ewald, num_compare, th_num, seed, check_tree);
    //printf("[timing=%d, verbose=%d, use_fft=%d, use_precompute=%d]\n", timing, verbose, use_fft, use_precompute);
    printf("input parameters=[\n");
    printf("  %-15s = %d\n","P",P);
    printf("  %-15s = %d\n","n",n);
    printf("  %-15s = %d\n","ncrit",ncrit);
    printf("  %-15s = %d\n","images",images);
    printf("  %-15s = %.4f\n","cycle",cycle);
    printf("  %-15s = %d\n","ewald_ksize",ewald_ksize);
    printf("  %-15s = (%d,%d,%d)\n","(f,d,e)",enable_fmm,enable_direct,enable_ewald);
    printf("  %-15s = %d\n","num_compare",num_compare);
    printf("  %-15s = %d\n","th_num",th_num);
    printf("  %-15s = %d\n","seed",seed);
    printf("  %-15s = %d\n","check_tree",check_tree);
    printf("  %-15s = %d\n","timing",timing);
    printf("  %-15s = %d\n","verbose",verbose);
    printf("  %-15s = %d\n","use_fft",use_fft);
    printf("  %-15s = %d\n","use_precompute",use_precompute);
    printf("]\n\n");

}