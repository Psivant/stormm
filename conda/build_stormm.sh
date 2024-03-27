#!/bin/bash

set -ex

CMAKE_FLAGS="${CMAKE_ARGS} -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_BUILD_TYPE=Release -DSTORMM_ENABLE_RDKIT=YES -DSTORMM_BUILD_TESTS=NO"

if [[ "$cuda_compiler_version" == "None" ]]; then
  CMAKE_FLAGS="${CMAKE_ARGS} -DSTORMM_ENABLE_CUDA=NO"
else
  CMAKE_FLAGS="${CMAKE_ARGS} -DSTORMM_ENABLE_CUDA=YES"
fi

mkdir -p build
cd build
cmake ${CMAKE_FLAGS} ..
make -j$CPU_COUNT
make -j$CPU_COUNT install
