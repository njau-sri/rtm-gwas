#!/bin/bash

VER=3.8.0

if [[ ! -f lapack-${VER}.tar.gz ]]; then
    curl -O http://www.netlib.org/lapack/lapack-${VER}.tar.gz || exit 1
fi

TOP=$(pwd)
LIBDIR=/usr/local/lib

rm -rf lapack-$VER
tar zxf lapack-${VER}.tar.gz
cd lapack-$VER

cp make.inc.example make.inc

make blaslib lapacklib || exit 1

sudo cp -f librefblas.a liblapack.a ${LIBDIR}/
cd $TOP
