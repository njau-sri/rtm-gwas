# Compiling Environment

## 1. Targetting Linux

CentOS 7 Host

    yum install zip unzip
    yum install gcc gcc-c++ glibc-static libstdc++-static libgfortran-static
    yum install boost-static qt5-qtbase-devel
    yum install blas-static lapack-static

## 2. Targetting Windows

Fedora 30 Host

### Windows 32-bit

    dnf install make mingw32-gcc mingw32-gcc-c++ mingw32-gcc-gfortran
    dnf install mingw32-libgomp mingw32-winpthreads-static
    dnf install mingw32-boost-static
    dnf install mingw32-qt-static

```
LAPACK=lapack-3.9.0

rm -rf $LAPACK
tar zxf $LAPACK.tar.gz
cd $LAPACK

cp make.inc.example make.inc
sed -i "s/= gfortran/= i686-w64-mingw32-gfortran/" make.inc
sed -i "s/= ar/= i686-w64-mingw32-ar/" make.inc
sed -i "s/= ranlib/= i686-w64-mingw32-ranlib/" make.inc
sed -i "s/librefblas/libblas/" make.inc

make blaslib lapacklib
cp libblas.a liblapack.a /usr/i686-w64-mingw32/sys-root/mingw/lib
```

### Windows 64-bit

    dnf install make mingw64-gcc mingw64-gcc-c++ mingw64-gcc-gfortran
    dnf install mingw64-libgomp mingw64-winpthreads-static
    dnf install mingw64-boost-static
    dnf install mingw64-qt-static

```
LAPACK=lapack-3.9.0

rm -rf $LAPACK
tar zxf $LAPACK.tar.gz
cd $LAPACK

cp make.inc.example make.inc
sed -i "s/= gfortran/= x86_64-w64-mingw32-gfortran/" make.inc
sed -i "s/= ar/= x86_64-w64-mingw32-ar/" make.inc
sed -i "s/= ranlib/= x86_64-w64-mingw32-ranlib/" make.inc
sed -i "s/librefblas/libblas/" make.inc

make blaslib lapacklib
cp libblas.a liblapack.a /usr/x86_64-w64-mingw32/sys-root/mingw/lib
```

## 3. Targetting macOS

macOS Mojave Host

### command line tools

    xcode-select --install

sdk header

    /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg

### homebrew

https://brew.sh/

    brew install llvm libomp boost qt

## 4. MKL

https://software.intel.com/en-us/intel-mkl
