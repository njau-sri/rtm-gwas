#!/bin/bash

VER=0.3.3

if [ ! -f OpenBLAS-${VER}.tar.gz ]; then
    wget -O OpenBLAS-${VER}.tar.gz https://github.com/xianyi/OpenBLAS/archive/v${VER}.tar.gz || exit 1
fi

if [[ $# -ne 1 ]]; then
    echo "usage: $0 mingw32|mingw64"
fi

TOP=$(pwd)

OPTS="USE_THREAD=0 DYNAMIC_ARCH=1 DYNAMIC_OLDER=0 NO_CBLAS=1 NO_LAPACKE=1 NO_SHARED=1"

if [[ $1 == mingw32 ]]; then
    OPTS="$OPTS NUM_THREADS=8 USE_NETLIB_GEMV=1"
    CROSS=i686-w64-mingw32-
    LIBDIR=/usr/i686-w64-mingw32/sys-root/mingw/lib
elif [[ $1 == mingw64 ]]; then
    OPTS="$OPTS NUM_THREADS=24"
    CROSS=x86_64-w64-mingw32-
    LIBDIR=/usr/x86_64-w64-mingw32/sys-root/mingw/lib
else
    echo "!!! ONLY FOR i686-w64-mingw32 AND x86_64-w64-mingw32 !!!"
    exit 1
fi

rm -rf OpenBLAS-$VER
tar zxf OpenBLAS-${VER}.tar.gz
cd OpenBLAS-$VER

make $OPTS HOSTCC=gcc CC=${CROSS}gcc FC=${CROSS}gfortran libs netlib || exit 1

cp -f libopenblas*-r${VER}.a ${LIBDIR}/libopenblas.a
cd $TOP
