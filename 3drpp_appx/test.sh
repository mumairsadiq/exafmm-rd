#!/bin/bash

if [ "$1" == "unit" ]; then
    echo "run unit test with parameter : ${@:2}"
    cd build/gtest
    ./test_fmm "${@:2}"
elif [ "$1" == "reg" ]; then
    echo "run regularization check with parameter : ${@:2}"
    cd build/gtest
    for z0 in $(seq -0.01 0.001 0.01)
    do
        ./test_fmm --gtest_filter=*basic --override 1 --body0_idx 0 --x0 0.1 --y0 0.1 --z0 $z0 ${@:2}
    done
    
fi