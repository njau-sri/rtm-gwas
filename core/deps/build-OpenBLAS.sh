#!/bin/bash

VER=0.2.18

if [ $# -ne 1 ]; then
    echo "usage: $0 lnx64|osx64|win32|win64"
    exit 1
fi

if [ -d OpenBLAS-$1 ]; then
    exit 0
fi

if [ ! -f OpenBLAS-$VER.tar.gz ]; then
    curl -fkL https://github.com/xianyi/OpenBLAS/archive/v$VER.tar.gz -o OpenBLAS-$VER.tar.gz || exit 1
fi

TOP=`pwd`

rm -rf OpenBLAS-$VER
tar zxf OpenBLAS-$VER.tar.gz
cd OpenBLAS-$VER

BUILD_OPTS="DYNAMIC_ARCH=1 USE_THREAD=1 NUM_THREADS=4 NO_CBLAS=1 NO_LAPACKE=1"

if [ $1 == "lnx64" ]; then
    make $BUILD_OPTS BINARY=64 USE_OPENMP=1 libs netlib || exit 1
elif [ $1 == "osx64" ]; then
    make $BUILD_OPTS BBINARY=64 libs netlib || exit 1
elif [ $1 == "win32" ]; then
    make $BUILD_OPTS BBINARY=32  USE_OPENMP=1 CC=i686-w64-mingw32-gcc FC=i686-w64-mingw32-gfortran HOSTCC=gcc libs netlib || exit 1
elif [ $1 == "win64" ]; then
    make $BUILD_OPTS BBINARY=64 USE_OPENMP=1 CC=x86_64-w64-mingw32-gcc FC=x86_64-w64-mingw32-gfortran HOSTCC=gcc libs netlib || exit 1
else
    echo "*** invalid target $1"
    exit 1
fi

mkdir $TOP/OpenBLAS-$1
cp libopenblas.a $TOP/OpenBLAS-$1
