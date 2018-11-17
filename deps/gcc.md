# 1 CentOS 7

[CentOS-7-x86_64-Minimal-1804.iso](http://mirrors.ustc.edu.cn/centos/7/isos/x86_64/CentOS-7-x86_64-Minimal-1804.iso)

### 1.1 Developer Tools

    yum update
    yum install epel-release
    yum install wget

    yum install gcc gcc-c++ gcc-gfortran
    yum install glibc-static libstdc++-static libgfortran-static
    yum install boost-static

    yum install mingw32-gcc mingw32-gcc-c++ mingw32-gcc-gfortran
    yum install mingw32-libgomp mingw32-winpthreads-static
    yum install mingw32-boost-static

    yum install mingw64-gcc mingw64-gcc-c++ mingw64-gcc-gfortran
    yum install mingw64-libgomp mingw64-winpthreads-static
    yum install mingw64-boost-static

### 1.2 Toolchain

    # Linux
    gcc
    g++
    gfortran
    strip
    
    # MINGW32
    i686-w64-mingw32-gcc
    i686-w64-mingw32-g++
    i686-w64-mingw32-gfortran
    x86_64-w64-mingw32-strip
    
    # MINGW64
    x86_64-w64-mingw32-gcc
    x86_64-w64-mingw32-g++
    x86_64-w64-mingw32-gfortran
    x86_64-w64-mingw32-strip

# 2 macOS High Sierra

### 2.1 Developer Tools

Install command line developer tools

    xcode-select --install

Install [gfortran-8.1-bin.tar.gz](http://hpc.sourceforge.net/)

    sudo tar -zxf gfortran-8.1-bin.tar.gz -C /

Install Homebrew

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    
    brew install boost

### 2.2 Toolchain

    gcc
    g++
    clang
    clang++
    gfortran
