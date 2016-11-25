#!/bin/bash
#
# Fedora
#   yum install gcc gcc-c++ gcc-gfortran glibc-static libstdc++-static libgfortran-static
#   yum install mingw32-gcc mingw32-gcc-c++ mingw32-gcc-gfortran
#   yum install mingw64-gcc mingw64-gcc-c++ mingw64-gcc-gfortran
#   yum install mingw32-winpthreads-static mingw32-libgomp
#   yum install mingw64-winpthreads-static mingw64-libgomp
#   yum install bzip2-devel xz-devel pcre-devel libcurl-devel
#
# OS X
#   xcode-select --install
#   INSTALL [gfortran]  http://hpc.sourceforge.net/
#   INSTALL [liblzma]   http://tukaani.org/xz/
#   INSTALL [pcre]      http://www.pcre.org/      --enable-utf
#

VER=3.3.0

if [ $# -ne 1 ]; then
    echo "usage: $0 lnx64|osx64|win32|win64"
    exit 1
fi

if [ -d Rmath-$1 ]; then
    exit 0
fi

if [ ! -f R-$VER.tar.gz ]; then
    curl -fkLO http://cran.rstudio.com/src/base/R-3/R-$VER.tar.gz || exit 1
fi

TOP=`pwd`

rm -rf R-$VER
tar zxf R-$VER.tar.gz
cd R-$VER

if [ $1 == "lnx64" ]; then
    ./configure --with-readline=no --with-x=no
    cd src/nmath/standalone
    make CFLAGS=-O2 static || exit 1
elif [ $1 == "osx64" ]; then
    ./configure --with-readline=no --with-x=no
    cd src/nmath/standalone
    make CFLAGS=-O2 static || exit 1
elif [ $1 == "win32" ]; then
    cd src/gnuwin32
    make MkRules || exit 1
    cd ../include
    make WIN=32 BINPREF=i686-w64-mingw32- EOPTS= -f Makefile.win config.h Rconfig.h Rmath.h || exit 1
    cd ../nmath/standalone
    make WIN=32 BINPREF=i686-w64-mingw32- EOPTS= -f Makefile.win shared implib || exit 1
elif [ $1 == "win64" ]; then
    cd src/gnuwin32
    make MkRules || exit 1
    cd ../include
    make WIN=64 BINPREF64=x86_64-w64-mingw32- EOPTS= -f Makefile.win config.h Rconfig.h Rmath.h || exit 1
    cd ../nmath/standalone
    make WIN=64 BINPREF64=x86_64-w64-mingw32- EOPTS= -f Makefile.win shared implib || exit 1
else
    echo "*** invalid target $1"
    exit 1
fi

mkdir $TOP/Rmath-$1
cp libRmath.a ../../include/Rmath.h $TOP/Rmath-$1

if [ $1 == "win32" ]; then
    cp libRmath.dll.a Rmath.dll $TOP/Rmath-$1
fi
