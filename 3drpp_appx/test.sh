#!/bin/bash

if [ "$1" == "unit" ]; then
    echo "run unit test with parameter : ${@:2}"
    cd build/gtest
    ./test_fmm "${@:2}"
elif [ "$1" == "reg" ]; then
    echo "run regularization check with parameter : ${@:2}"
    cd build/gtest
    ./test_fmm --gtest_filter=*basic "${@:2}"
fi