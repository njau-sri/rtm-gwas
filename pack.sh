#!/bin/bash -x

VER=2019.5.dev

if [[ $# -lt 1 ]]; then
    echo usage: $0 win32/win64/lnx64/macos
    exit 0
fi

PKG=rtm-gwas-$VER-$1

CC=gcc
CXX=g++

CFLAGS="-O2 -std=c11 -DNDEBUG -Wall"
CXXFLAGS="-O2 -std=c++11 -fopenmp -DNDEBUG -Wall"

LDFLAGS="-fopenmp"

MKL=
LAPACK="-llapack -lblas -lgfortran -lquadmath"

if [[ $1 == win32 ]]; then

    CC=i686-w64-mingw32-gcc
    CXX=i686-w64-mingw32-g++

    LDFLAGS="$LDFLAGS -s -static"

elif [[ $1 == win64 ]]; then

    CC=x86_64-w64-mingw32-gcc
    CXX=x86_64-w64-mingw32-g++

    LDFLAGS="$LDFLAGS -s -static"

elif [[ $1 == lnx64 ]]; then

    LDFLAGS="$LDFLAGS -s -static"

    MKL1=/opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a
    MKL2=/opt/intel/mkl/lib/intel64/libmkl_sequential.a
    MKL3=/opt/intel/mkl/lib/intel64/libmkl_core.a
    MKL="-Wl,--start-group $MKL1 $MKL2 $MKL3 -Wl,--end-group -lpthread -lm -ldl"

elif [[ $1 == macos ]]; then

    export PATH="/usr/local/opt/llvm/bin:$PATH"
    export PATH="/usr/local/opt/qt/bin:$PATH"

    CC=clang
    CXX=clang++

    CXXFLAGS="$CXXFLAGS -I/usr/local/include"

    LDFLAGS="$LDFLAGS -L/usr/local/lib"

    LAPACK="-framework Accelerate"

    MKL1=/opt/intel/mkl/lib/libmkl_intel_lp64.a
    MKL2=/opt/intel/mkl/lib/libmkl_sequential.a
    MKL3=/opt/intel/mkl/lib/libmkl_core.a
    MKL="$MKL1 $MKL2 $MKL3 -lpthread -lm -ldl"

else

    echo usage: $0 win32/win64/lnx64/mac64
    exit 1

fi

rm -rf $PKG
mkdir $PKG
cd $PKG

for src in ../rtm-gwas-snpldb/src/*.cpp; do
    $CXX -c $CXXFLAGS $src || exit 1
done
$CXX $LDFLAGS -o rtm-gwas-snpldb *.o || exit 1
rm *.o

for src in ../rtm-gwas-gsc/src/*.c; do
    $CC -c $CFLAGS $src || exit 1
done
for src in ../rtm-gwas-gsc/src/*.cpp; do
    $CXX -c $CXXFLAGS $src || exit 1
done
$CXX $LDFLAGS -o rtm-gwas-gsc *.o $LAPACK || exit 1
rm *.o

for src in ../rtm-gwas-assoc/src/*.c; do
    $CC -c $CFLAGS $src || exit 1
done
for src in ../rtm-gwas-assoc/src/*.cpp; do
    $CXX -c $CXXFLAGS $src || exit 1
done
if [[ $1 == lnx64 ]] || [[ $1 == macos ]]; then
    $CXX $LDFLAGS -o rtm-gwas-assoc *.o $MKL || exit 1
else
    $CXX $LDFLAGS -o rtm-gwas-assoc *.o $LAPACK || exit 1
fi
rm *.o

if [[ $1 = win32 ]]; then

    mingw32-qmake-qt4 "CONFIG += static" ../rtm-gwas-gui/src || exit 1
    make release || exit 1
    mv release/rtm-gwas-gui.exe .
    mingw-strip -s rtm-gwas-gui.exe
    make distclean
    rm -rf debug release object_script*

elif [[ $1 = win64 ]]; then

    mingw64-qmake-qt4 "CONFIG += static" ../rtm-gwas-gui/src || exit 1
    make release || exit 1
    mv release/rtm-gwas-gui.exe .
    mingw-strip -s rtm-gwas-gui.exe
    make distclean
    rm -rf debug release object_script*

elif [[ $1 = lnx64 ]]; then

    qmake-qt5 ../rtm-gwas-gui/src || exit 1
    make || exit 1
    strip -s rtm-gwas-gui
    make clean
    rm Makefile

elif [[ $1 == "macos" ]]; then

    qmake ../rtm-gwas-gui/src || exit 1
    make || exit 1
    macdeployqt rtm-gwas-gui.app
    make clean
    rm Makefile

fi
