#include "argument.h"
#include <fstream>

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
    cmd.add<int>("verbose", 'v', "is verbose", false, 0);
    cmd.add<real>("cycle", 'c', "cycle of images(or box_r)", false, 2 * M_PI);
    cmd.add<int>("num_compare", 0, "compare number of target body(also direct calculation number)", false, -1);
    cmd.add<int>("ewald_ksize", 0, "ksize of ewald DFT", false, 11);
    cmd.add<int>("th_num", 0, "number of omp threads", false, 8);
    cmd.add<int>("seed", 0, "random seed", false, 5);
    cmd.add<std::string>("algorithm", 'a', "algorithms", false, "fde");
    cmd.add<int>("check_tree", 0, "if check tree", false, 0);
    cmd.add<int>("use_fft", 0, "if use fft in M2L", false, 1);
    cmd.add<int>("use_precompute", 0, "if use precomputed kernel", false, 1);
    cmd.add<int>("use_simd", 0, "if use simd in P2P", false, 1);
    cmd.add<int>("dipole_correction", 0, "if use dipole correction for FMM/direct", false, 1);
    cmd.add<int>("zero_netcharge", 0, "if zero net charge", false, 1);
    cmd.add<int>("print_body_number", 0, "print body number", false, 3);
    cmd.add<int>("divide_4pi", 0, "if divide 4pi", false, 0);
    cmd.add<int>("setting_t", 0, "set some parameters to match the result with exafmm-t", false, 0);
    cmd.add<real>("rega", 'r', "regularization alpha", false, 0);
    cmd.add<int>("override", 0, "override gtest setting", false, 0);
    cmd.add<real>("x0", 0, "x0", false, 0.0);
    cmd.add<real>("y0", 0, "y0", false, 0.0);
    cmd.add<real>("z0", 0, "z0", false, 0.0);
    cmd.add<int>("body0_idx", 0, "index of body moved(for reg plot check)", false, -1);
    cmd.add<int>("check_body_idx", 0, "index of body check(for reg plot check)", false, 0);
    cmd.add<std::string>("reg_image0_type", 0, "type of fmm when reg=true,image=false", false, "c");
    cmd.add<std::string>("res_filepath", 0, "path of result file", true, "");

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
    dipole_correction = cmd.get<int>("dipole_correction");
    zero_netcharge = cmd.get<int>("zero_netcharge");
    print_body_number = cmd.get<int>("print_body_number");
    divide_4pi = cmd.get<int>("divide_4pi");
    use_simd = cmd.get<int>("use_simd");
    setting_t = cmd.get<int>("setting_t");
    rega = cmd.get<real>("rega");

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
    override_gtest_setting = cmd.get<int>("override");
    check_body_idx = cmd.get<int>("check_body_idx");

    x0 = cmd.get<real>("x0");
    y0 = cmd.get<real>("y0");
    z0 = cmd.get<real>("z0");
    body0_idx = cmd.get<int>("body0_idx");

    reg_image0_type = cmd.get<std::string>("reg_image0_type");

    res_filepath = cmd.get<std::string>("res_filepath");

    x = vec3r(0,0,0);
    r = cycle / 2;

    if(setting_t)
    {
        dipole_correction = 0;
        zero_netcharge = 0;
        divide_4pi = 1;
    }
}

void rtfmm::Argument::show(std::string res_filepath)
{
    std::string args = "";
    args += rtfmm::format("input parameters=[\n");
    args += rtfmm::format("\t%-20s = %d\n","override",override_gtest_setting);
    args += rtfmm::format("\t%-20s = %d\n","P",P);
    args += rtfmm::format("\t%-20s = %d\n","n",n);
    args += rtfmm::format("\t%-20s = %d\n","ncrit",ncrit);
    args += rtfmm::format("\t%-20s = %d\n","images",images);
    args += rtfmm::format("\t%-20s = %.4f\n","cycle",cycle);
    args += rtfmm::format("\t%-20s = %.4f\n","rega",rega);
    args += rtfmm::format("\t%-20s = %d\n","ewald_ksize",ewald_ksize);
    args += rtfmm::format("\t%-20s = (%d,%d,%d)\n","(f,d,e)",enable_fmm,enable_direct,enable_ewald);
    args += rtfmm::format("\t%-20s = %d\n","num_compare",num_compare);
    args += rtfmm::format("\t%-20s = %d\n","th_num",th_num);
    args += rtfmm::format("\t%-20s = %d\n","seed",seed);
    args += rtfmm::format("\t%-20s = %d\n","check_tree",check_tree);
    args += rtfmm::format("\t%-20s = %d\n","timing",timing);
    args += rtfmm::format("\t%-20s = %d\n","verbose",verbose);
    args += rtfmm::format("\t%-20s = %d\n","use_fft",use_fft);
    args += rtfmm::format("\t%-20s = %d\n","use_precompute",use_precompute);
    args += rtfmm::format("\t%-20s = %d\n","use_simd",use_simd);
    args += rtfmm::format("\t%-20s = %d\n","dipole_correction",dipole_correction);
    args += rtfmm::format("\t%-20s = %d\n","zero_netcharge",zero_netcharge);
    args += rtfmm::format("\t%-20s = %d\n","print_body_number",print_body_number);
    args += rtfmm::format("\t%-20s = %d\n","divide_4pi",divide_4pi);
    args += rtfmm::format("\t%-20s = %d\n","setting_t",setting_t);
    args += rtfmm::format("\t%-20s = %s\n","reg_image0_type",reg_image0_type.c_str());
    args += rtfmm::format("\t%-20s = [%.6f,%.6f,%.6f] for body[%d], check body[%d]\n","[x0,y0,z0]",x0,y0,z0,body0_idx,check_body_idx);
    args += rtfmm::format("\t%-20s = %s\n","res_filepath",res_filepath.c_str());
    args += rtfmm::format("]");

    std::cout << args << std::endl << std::endl;

    std::ofstream outputfile(res_filepath);
    outputfile << args << std::endl << std::endl;
    outputfile.close();
}