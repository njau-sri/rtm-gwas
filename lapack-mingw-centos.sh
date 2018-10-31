#!/bin/bash

VER=3.8.0

if [[ ! -f lapack-${VER}.tar.gz ]]; then
    wget http://www.netlib.org/lapack/lapack-${VER}.tar.gz || exit 1
fi

if [[ $# -ne 1 ]]; then
    echo "usage: $0 mingw32|mingw64"
fi

TOP=$(pwd)

CROSS=
LIBDIR=
if [[ $1 == mingw32 ]]; then
    CROSS=i686-w64-mingw32-
    LIBDIR=/usr/i686-w64-mingw32/sys-root/mingw/lib
elif [[ $1 == mingw64 ]]; then
    CROSS=x86_64-w64-mingw32-
    LIBDIR=/usr/x86_64-w64-mingw32/sys-root/mingw/lib
else
    echo "!!! ONLY FOR i686-w64-mingw32 AND x86_64-w64-mingw32 !!!"
    exit 1
fi

rm -rf lapack-$VER
tar zxf lapack-${VER}.tar.gz
cd lapack-$VER

cp make.inc.example make.inc
sed -i "s/= gfortran/= ${CROSS}gfortran/" make.inc
sed -i "s/= ar/= ${CROSS}ar/" make.inc
sed -i "s/= ranlib/= ${CROSS}ranlib/" make.inc

make blaslib lapacklib || exit 1

cp librefblas.a liblapack.a ${LIBDIR}/
cd $TOP
