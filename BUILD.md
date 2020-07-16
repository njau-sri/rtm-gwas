# Compiling Environment

## Linux

gcc, g++, boost, blas, lapack, eigen3, qt5

## Windows

MSYS2 MinGW

http://www.msys2.org/

Windows 32-bit

    pacman -S mingw-w64-i686-make mingw-w64-i686-gcc
    pacman -S mingw-w64-i686-boost
    pacman -S mingw-w64-i686-lapack mingw-w64-i686-eigen3
    pacman -S mingw-w64-i686-qt5

Windows 64-bit

    pacman -S mingw-w64-x86_64-make mingw-w64-x86_64-gcc
    pacman -S mingw-w64-x86_64-boost
    pacman -S mingw-w64-x86_64-lapack mingw-w64-x86_64-eigen3
    pacman -S mingw-w64-x86_64-qt5

## macOS Mojave

### command line tools

    xcode-select --install

sdk header

    /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg

### homebrew

https://brew.sh/

    brew install llvm@9 libomp boost eigen qt
