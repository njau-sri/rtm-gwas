# CentOS 7

### Linux NETLIB-LAPACK, OpenBLAS

    yum install blas-static lapack-static
    yum install openblas-static

### Windows NETLIB-LAPACK

[lapack-mingw-centos.sh](lapack-mingw-centos.sh)

### Windows OpenBLAS

[openblas-mingw-centos.sh](openblas-mingw-centos.sh)

# macOS High Sierra

### NETLIB-LAPACK

[lapack-macos.sh](lapack-macos.sh)

Static linking

    LIB=/usr/local/lib
    LIBGFORTRAN=${LIB}/libgfortran.a
    LIBQUADMATH=${LIB}/libquadmath.a
    LIBGCC=${LIB}/gcc/x86_64-apple-darwin17.5.0/8.1.0/libgcc.a

    -llapack -lrefblas $LIBGFORTRAN $LIBQUADMATH $LIBGCC


### OpenBLAS

[openblas-macos.sh](openblas-macos.sh)

Static linking

    LIB=/usr/local/lib
    LIBGFORTRAN=${LIB}/libgfortran.a
    LIBQUADMATH=${LIB}/libquadmath.a
    LIBGCC=${LIB}/gcc/x86_64-apple-darwin17.5.0/8.1.0/libgcc.a

    -lopenblas $LIBGFORTRAN $LIBQUADMATH $LIBGCC
