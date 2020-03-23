# Compiling Environment

## Linux

gcc, g++, boost, qt5, blas, lapack

## Windows

MSYS2 MinGW

http://www.msys2.org/

Windows 32-bit

    pacman -S make mingw-w64-i686-gcc mingw-w64-i686-gcc-fortran mingw-w64-i686-lapack mingw-w64-i686-boost

Windows 64-bit

    pacman -S make mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-lapack mingw-w64-x86_64-boost

## macOS

### command line tools

    xcode-select --install

sdk header

    /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg

### homebrew

https://brew.sh/

    brew install llvm libomp boost qt
