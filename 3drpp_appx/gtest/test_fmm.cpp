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


TEST(FmmTest, default_param) 
{
    rtfmm::title("rtfmm_3dpp_test_reg 1");
    EXPECT_EQ("empty", "empty") << "string error!";    
}

TEST(FmmTest, empty_string) 
{
    rtfmm::title("rtfmm_3dpp_test_reg 2");
    EXPECT_EQ("empty",  "empty");
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    rtfmm::Argument args(argc, argv);
    args.show();
    return RUN_ALL_TESTS();
}
