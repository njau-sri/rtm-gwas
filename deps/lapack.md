# Math Kernel Library (MKL)

https://software.intel.com/en-us/mkl

# NETLIB-LAPACK and OpenBLAS

### CentOS 7

For Linux

    yum install blas-static lapack-static
    yum install openblas-static

For Windows

[lapack-mingw-centos.sh](lapack-mingw-centos.sh)

[openblas-mingw-centos.sh](openblas-mingw-centos.sh)

### macOS

NETLIB-LAPACK

[lapack-macos.sh](lapack-macos.sh)

Static linking

    LIB=/usr/local/lib
    LIBGFORTRAN=${LIB}/libgfortran.a
    LIBQUADMATH=${LIB}/libquadmath.a
    LIBGCC=${LIB}/gcc/x86_64-apple-darwin17.5.0/8.1.0/libgcc.a

    -llapack -lrefblas $LIBGFORTRAN $LIBQUADMATH $LIBGCC

OpenBLAS

[openblas-macos.sh](openblas-macos.sh)

Static linking

    LIB=/usr/local/lib
    LIBGFORTRAN=${LIB}/libgfortran.a
    LIBQUADMATH=${LIB}/libquadmath.a
    LIBGCC=${LIB}/gcc/x86_64-apple-darwin17.5.0/8.1.0/libgcc.a

    -lopenblas $LIBGFORTRAN $LIBQUADMATH $LIBGCC
