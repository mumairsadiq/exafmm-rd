#include "gtest/gtest.h"
#include <stdexcept>
#include <iostream>
#include "type.h"
#include "body.h"
#include "tree.h"
#include "fmm.h"
#include "argument.h"
#include "ewald.h"
#include <omp.h>


TEST(FmmTest, n100_p4) 
{
    rtfmm::Argument args;

    args.n = 1000;
    args.P = 4;
    args.images = 0;
    args.rega = 0;

    args.show();

    omp_set_dynamic(0);
    omp_set_num_threads(args.th_num);
    RTLOG("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x, args.seed, args.zero_netcharge);

    /* solve by FMM */
    rtfmm::Bodies3 res_fmm;
    if(args.enable_fmm)
    {
        rtfmm::LaplaceFMM fmm(bs, args);
        TIME_BEGIN(FMM);
        res_fmm = fmm.solve();
        if(args.timing) {TIME_END_stdout(FMM);}
    }

    /* solve by ewald */
    rtfmm::Bodies3 res_ewald;
    if(args.enable_ewald)
    {
        rtfmm::EwaldSolver ewald(bs, args);
        TIME_BEGIN(EWALD);
        res_ewald = ewald.solve();
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_ewald);
        if(args.timing) {TIME_END_stdout(EWALD);}
    }

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    if(args.enable_direct)
    {
        rtfmm::LaplaceKernel kernel;
        rtfmm::Cell3 cell_src;
        cell_src.brange = {0, args.n};
        rtfmm::Cell3 cell_tar;
        cell_tar.brange = {0, args.num_compare};
        TIME_BEGIN(DIRECT);
        kernel.direct(res_direct, res_direct, args.images, args.cycle);
        if(args.dipole_correction)
            rtfmm::dipole_correction(res_direct, args.cycle);
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_direct);
        if(args.timing) {TIME_END_stdout(DIRECT);}
    }

    /* compare */
    if(rtfmm::verbose)
    {
        if(args.enable_fmm) rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "fmm");
        if(args.enable_direct) rtfmm::print_bodies(res_direct, args.print_body_number, 0, "direct");
        if(args.enable_ewald) rtfmm::print_bodies(res_ewald, args.print_body_number, 0, "ewald");
    }
    if(args.enable_fmm && args.enable_direct)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_fmm, res_direct, "FMM", "Direct", args.num_compare);
        res.show();
        EXPECT_LE(res.l2f, 1.1e-4);
        EXPECT_LE(res.l2e, 2.9e-4);
    }
    if(args.enable_fmm && args.enable_ewald)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_fmm, res_ewald, "FMM", "Ewald", args.num_compare);
        res.show();
    }
    if(args.enable_direct && args.enable_ewald)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_direct, res_ewald, "Direct", "Ewald", args.num_compare);
        res.show();
    }
}

TEST(FmmTest, n24000_p6) 
{
    rtfmm::Argument args;

    args.n = 24000;
    args.num_compare = 24000;
    args.P = 6;
    args.images = 0;
    args.rega = 0;
    args.ncrit = 128;
    args.dipole_correction = 0;
    args.zero_netcharge = 0;
    args.divide_4pi = 1;
    args.setting_t = 1;

    args.show();

    omp_set_dynamic(0);
    omp_set_num_threads(args.th_num);
    RTLOG("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x, args.seed, args.zero_netcharge);

    /* solve by FMM */
    rtfmm::Bodies3 res_fmm;
    if(args.enable_fmm)
    {
        rtfmm::LaplaceFMM fmm(bs, args);
        TIME_BEGIN(FMM);
        res_fmm = fmm.solve();
        if(args.timing) {TIME_END_stdout(FMM);}
    }

    /* solve by ewald */
    rtfmm::Bodies3 res_ewald;
    if(args.enable_ewald)
    {
        rtfmm::EwaldSolver ewald(bs, args);
        TIME_BEGIN(EWALD);
        res_ewald = ewald.solve();
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_ewald);
        if(args.timing) {TIME_END_stdout(EWALD);}
    }

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    if(args.enable_direct)
    {
        rtfmm::LaplaceKernel kernel;
        rtfmm::Cell3 cell_src;
        cell_src.brange = {0, args.n};
        rtfmm::Cell3 cell_tar;
        cell_tar.brange = {0, args.num_compare};
        TIME_BEGIN(DIRECT);
        kernel.direct(res_direct, res_direct, args.images, args.cycle);
        if(args.dipole_correction)
            rtfmm::dipole_correction(res_direct, args.cycle);
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_direct);
        if(args.timing) {TIME_END_stdout(DIRECT);}
    }

    /* compare */
    if(rtfmm::verbose)
    {
        if(args.enable_fmm) rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "fmm");
        if(args.enable_direct) rtfmm::print_bodies(res_direct, args.print_body_number, 0, "direct");
        if(args.enable_ewald) rtfmm::print_bodies(res_ewald, args.print_body_number, 0, "ewald");
    }
    if(args.enable_fmm && args.enable_direct)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_fmm, res_direct, "FMM", "Direct", args.num_compare);
        res.show();
        EXPECT_LE(res.l2f, 5e-6);
        EXPECT_LE(res.l2e, 5e-7);
    }
    if(args.enable_fmm && args.enable_ewald)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_fmm, res_ewald, "FMM", "Ewald", args.num_compare);
        res.show();
    }
    if(args.enable_direct && args.enable_ewald)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_direct, res_ewald, "Direct", "Ewald", args.num_compare);
        res.show();
    }
}

TEST(FmmTest, n24000_p10_image5) 
{
    rtfmm::Argument args;

    args.enable_direct = 0;
    args.n = 24000;
    args.num_compare = 24000;
    args.P = 10;
    args.images = 5;
    args.rega = 0;
    args.ncrit = 128;
    args.dipole_correction = 1;
    args.zero_netcharge = 1;
    args.divide_4pi = 0;
    args.setting_t = 0;

    args.show();

    omp_set_dynamic(0);
    omp_set_num_threads(args.th_num);
    RTLOG("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x, args.seed, args.zero_netcharge);

    /* solve by FMM */
    rtfmm::Bodies3 res_fmm;
    if(args.enable_fmm)
    {
        rtfmm::LaplaceFMM fmm(bs, args);
        TIME_BEGIN(FMM);
        res_fmm = fmm.solve();
        if(args.timing) {TIME_END_stdout(FMM);}
    }

    /* solve by ewald */
    rtfmm::Bodies3 res_ewald;
    if(args.enable_ewald)
    {
        rtfmm::EwaldSolver ewald(bs, args);
        TIME_BEGIN(EWALD);
        res_ewald = ewald.solve();
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_ewald);
        if(args.timing) {TIME_END_stdout(EWALD);}
    }

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    if(args.enable_direct)
    {
        rtfmm::LaplaceKernel kernel;
        rtfmm::Cell3 cell_src;
        cell_src.brange = {0, args.n};
        rtfmm::Cell3 cell_tar;
        cell_tar.brange = {0, args.num_compare};
        TIME_BEGIN(DIRECT);
        kernel.direct(res_direct, res_direct, args.images, args.cycle);
        if(args.dipole_correction)
            rtfmm::dipole_correction(res_direct, args.cycle);
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_direct);
        if(args.timing) {TIME_END_stdout(DIRECT);}
    }

    /* compare */
    if(rtfmm::verbose)
    {
        if(args.enable_fmm) rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "fmm");
        if(args.enable_direct) rtfmm::print_bodies(res_direct, args.print_body_number, 0, "direct");
        if(args.enable_ewald) rtfmm::print_bodies(res_ewald, args.print_body_number, 0, "ewald");
    }
    if(args.enable_fmm && args.enable_direct)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_fmm, res_direct, "FMM", "Direct", args.num_compare);
        res.show();
    }
    if(args.enable_fmm && args.enable_ewald)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_fmm, res_ewald, "FMM", "Ewald", args.num_compare);
        EXPECT_LE(res.l2f, 3e-7);
        EXPECT_LE(res.l2e, 2e-5);
        res.show();
    }
    if(args.enable_direct && args.enable_ewald)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_direct, res_ewald, "Direct", "Ewald", args.num_compare);
        res.show();
    }
}

TEST(FmmTest, n24000_p6_reg0001) 
{
    rtfmm::Argument args;

    args.n = 24000;
    args.num_compare = 24000;
    args.P = 6;
    args.images = 0;
    args.rega = 0.001;
    args.ncrit = 128;
    args.dipole_correction = 0;
    args.zero_netcharge = 0;
    args.divide_4pi = 1;
    args.setting_t = 1;

    args.show();

    omp_set_dynamic(0);
    omp_set_num_threads(args.th_num);
    RTLOG("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x, args.seed, args.zero_netcharge);

    /* solve by FMM */
    rtfmm::Bodies3 res_fmm;
    if(args.enable_fmm)
    {
        rtfmm::LaplaceFMM fmm(bs, args);
        TIME_BEGIN(FMM);
        res_fmm = fmm.solve();
        if(args.timing) {TIME_END_stdout(FMM);}
    }

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    if(args.enable_direct)
    {
        rtfmm::LaplaceKernel kernel;
        rtfmm::Cell3 cell_src;
        cell_src.brange = {0, args.n};
        rtfmm::Cell3 cell_tar;
        cell_tar.brange = {0, args.num_compare};
        TIME_BEGIN(DIRECT);
        kernel.direct(res_direct, res_direct, args.images, args.cycle);
        if(args.dipole_correction)
            rtfmm::dipole_correction(res_direct, args.cycle);
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_direct);
        if(args.timing) {TIME_END_stdout(DIRECT);}
    }

    /* compare */
    if(rtfmm::verbose)
    {
        if(args.enable_fmm) rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "fmm");
        if(args.enable_direct) rtfmm::print_bodies(res_direct, args.print_body_number, 0, "direct");
    }
    if(args.enable_fmm && args.enable_direct)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_fmm, res_direct, "FMM", "Direct", args.num_compare);
        res.show();
        EXPECT_LE(res.l2f, 5e-6);
        EXPECT_LE(res.l2e, 5e-7);
    }
}

TEST(FmmTest, n24000_p6_reg001) 
{
    rtfmm::Argument args;

    args.n = 24000;
    args.num_compare = 24000;
    args.P = 6;
    args.images = 0;
    args.rega = 0.01;
    args.ncrit = 128;
    args.dipole_correction = 0;
    args.zero_netcharge = 0;
    args.divide_4pi = 1;
    args.setting_t = 1;

    args.show();

    omp_set_dynamic(0);
    omp_set_num_threads(args.th_num);
    RTLOG("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x, args.seed, args.zero_netcharge);

    /* solve by FMM */
    rtfmm::Bodies3 res_fmm;
    if(args.enable_fmm)
    {
        rtfmm::LaplaceFMM fmm(bs, args);
        TIME_BEGIN(FMM);
        res_fmm = fmm.solve();
        if(args.timing) {TIME_END_stdout(FMM);}
    }

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    if(args.enable_direct)
    {
        rtfmm::LaplaceKernel kernel;
        rtfmm::Cell3 cell_src;
        cell_src.brange = {0, args.n};
        rtfmm::Cell3 cell_tar;
        cell_tar.brange = {0, args.num_compare};
        TIME_BEGIN(DIRECT);
        kernel.direct(res_direct, res_direct, args.images, args.cycle);
        if(args.dipole_correction)
            rtfmm::dipole_correction(res_direct, args.cycle);
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_direct);
        if(args.timing) {TIME_END_stdout(DIRECT);}
    }

    /* compare */
    if(rtfmm::verbose)
    {
        if(args.enable_fmm) rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "fmm");
        if(args.enable_direct) rtfmm::print_bodies(res_direct, args.print_body_number, 0, "direct");
    }
    if(args.enable_fmm && args.enable_direct)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_fmm, res_direct, "FMM", "Direct", args.num_compare);
        res.show();
        EXPECT_LE(res.l2f, 5e-6);
        EXPECT_LE(res.l2e, 5e-7);
    }
}

TEST(FmmTest, n24000_p6_reg01) 
{
    rtfmm::Argument args;

    args.n = 24000;
    args.num_compare = 24000;
    args.P = 6;
    args.images = 0;
    args.rega = 0.1;
    args.ncrit = 128;
    args.dipole_correction = 0;
    args.zero_netcharge = 0;
    args.divide_4pi = 1;
    args.setting_t = 1;

    args.show();

    omp_set_dynamic(0);
    omp_set_num_threads(args.th_num);
    RTLOG("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x, args.seed, args.zero_netcharge);

    /* solve by FMM */
    rtfmm::Bodies3 res_fmm;
    if(args.enable_fmm)
    {
        rtfmm::LaplaceFMM fmm(bs, args);
        TIME_BEGIN(FMM);
        res_fmm = fmm.solve();
        if(args.timing) {TIME_END_stdout(FMM);}
    }

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    if(args.enable_direct)
    {
        rtfmm::LaplaceKernel kernel;
        rtfmm::Cell3 cell_src;
        cell_src.brange = {0, args.n};
        rtfmm::Cell3 cell_tar;
        cell_tar.brange = {0, args.num_compare};
        TIME_BEGIN(DIRECT);
        kernel.direct(res_direct, res_direct, args.images, args.cycle);
        if(args.dipole_correction)
            rtfmm::dipole_correction(res_direct, args.cycle);
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_direct);
        if(args.timing) {TIME_END_stdout(DIRECT);}
    }

    /* compare */
    if(rtfmm::verbose)
    {
        if(args.enable_fmm) rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "fmm");
        if(args.enable_direct) rtfmm::print_bodies(res_direct, args.print_body_number, 0, "direct");
    }
    if(args.enable_fmm && args.enable_direct)
    {
        rtfmm::BodyCompareResult res = rtfmm::compare(res_fmm, res_direct, "FMM", "Direct", args.num_compare);
        res.show();
        EXPECT_LE(res.l2f, 5e-6);
        EXPECT_LE(res.l2e, 5e-7);
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
