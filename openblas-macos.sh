#!/bin/bash

VER=0.3.3

if [ ! -f OpenBLAS-${VER}.tar.gz ]; then
    curl -fkLo OpenBLAS-$VER.tar.gz https://github.com/xianyi/OpenBLAS/archive/v$VER.tar.gz || exit 1
fi

TOP=$(pwd)

LIBDIR=/usr/local/lib

OPTS="DYNAMIC_ARCH=1 DYNAMIC_OLDER=0 NO_CBLAS=1 NO_LAPACKE=1 NO_SHARED=1"

THREAD0="USE_THREAD=0"
THREAD1="USE_THREAD=1 NUM_THREADS=4 GEMM_MULTITHREADING_THRESHOLD=50 NO_AFFINITY=1"

rm -rf OpenBLAS-$VER
tar zxf OpenBLAS-${VER}.tar.gz
cd OpenBLAS-$VER

make $OPTS $THREAD1 libs netlib || exit 1

sudo cp libopenblas*-r${VER}.a ${LIBDIR}/libopenblas.a
cd $TOP
