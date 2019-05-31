# Compiling Environment

## 1 Targetting Linux and Windows

CentOS 7 Host

### Install CentOS 7

http://mirrors.ustc.edu.cn/centos/7/isos/x86_64/CentOS-7-x86_64-Minimal-1810.iso

```shell
yum install epel-release
yum update
yum install wget git zip unzip
```

Menu Player -> Manage -> Install VMware Tools

```shell
mkdir cdrom
mount /dev/cdrom cdrom/
cp cdrom/VMwareTools-*.tar.gz .
umount cdrom/
rmdir cdrom

yum install open-vm-tools
tar -zxf VMwareTools-*.tar.gz
cd vmware-tools-distrib
./vmware-install.pl -d

cd
echo 'vmhgfs-fuse .host: /mnt/hgfs' >> .bash_profile
```

### Targetting Linux 64-bit

**gcc**

```shell
yum install gcc gcc-c++ gcc-gfortran
yum install glibc-static libstdc++-static libgfortran-static
```

**boost**

```shell
yum install boost-static
```

**lapack**

```shell
yum install blas-static lapack-static
```

**qt5**

```shell
yum install qt5-qtbase-devel
```

**mkl**

```shell
MKLROOT=/opt/intel/mkl
MKLLIB=${MKLROOT}/lib/intel64

LIB1=${MKLLIB}/libmkl_intel_lp64.a
LIB2=${MKLLIB}/libmkl_sequential.a
LIB3=${MKLLIB}/libmkl_core.a

LDLIBS += -Wl,--start-group $LIB1 $LIB2 $LIB3 -Wl,--end-group -lpthread -lm -ldl
```

### Targetting Windows 32-bit

**gcc**

```shell
yum install mingw32-gcc mingw32-gcc-c++ mingw32-gcc-gfortran
yum install mingw32-libgomp mingw32-winpthreads-static
```

**boost**

```shell
yum install mingw32-boost-static
```

**lapack**

```shell
LAPACK=lapack-3.8.0

rm -rf $LAPACK
tar zxf $LAPACK.tar.gz
cd $LAPACK

cp make.inc.example make.inc
sed -i "s/= gfortran/= i686-w64-mingw32-gfortran/" make.inc
sed -i "s/= ar/= i686-w64-mingw32-ar/" make.inc
sed -i "s/= ranlib/= i686-w64-mingw32-ranlib/" make.inc

make blaslib lapacklib
cp librefblas.a liblapack.a /usr/i686-w64-mingw32/sys-root/mingw/lib
```

**qt4**

```bash
yum install mingw32-qt-static
```

### Targetting Windows 64-bit

**gcc**

```shell
yum install mingw64-gcc mingw64-gcc-c++ mingw64-gcc-gfortran
yum install mingw64-libgomp mingw64-winpthreads-static
```

**boost**

```shell
yum install mingw64-boost-static
```

**lapack**

```sh
LAPACK=lapack-3.8.0

rm -rf $LAPACK
tar zxf $LAPACK.tar.gz
cd $LAPACK

cp make.inc.example make.inc
sed -i "s/= gfortran/= x86_64-w64-mingw32-gfortran/" make.inc
sed -i "s/= ar/= x86_64-w64-mingw32-ar/" make.inc
sed -i "s/= ranlib/= x86_64-w64-mingw32-ranlib/" make.inc

make blaslib lapacklib
cp librefblas.a liblapack.a /usr/x86_64-w64-mingw32/sys-root/mingw/lib
```

**qt4**

```shell
yum install mingw64-qt-static
```

## 2 Targetting macOS

macOS Mojave Host

### vmware

unlocker: https://github.com/DrDonk/unlocker

add `smc.version = 0` line to *.vmx file

### command line developer tools

```shell
xcode-select --install
```

### homebrew

https://brew.sh/

https://raw.githubusercontent.com/Homebrew/install/master/install

replace

> BREW_REPO = "https://github.com/Homebrew/brew"

with

> BREW_REPO = "http://mirrors.ustc.edu.cn/brew.git/"

```shell
export HOMEBREW_BOTTLE_DOMAIN=https://mirrors.ustc.edu.cn/homebrew-bottles
```

### llvm

```shell
brew install llvm@7 libomp
rm /usr/local/opt/libomp/lib/*.dylib
```

sdk header

```shell
open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg
```

clang

```shell
export PATH="/usr/local/opt/llvm@7/bin:$PATH"
export LDFLAGS="-L/usr/local/opt/llvm@7/lib"
export CPPFLAGS="-I/usr/local/opt/llvm@7/include"
```

openmp

```shell
-fopenmp
```

lapack

```shell
-framework Accelerate
```

### boost

```shell
brew install boost
```

### mkl

```shell
MKLROOT=/opt/intel/mkl
MKLLIB=${MKLROOT}/lib

LIB1=${MKLLIB}/libmkl_intel_lp64.a
LIB2=${MKLLIB}/libmkl_sequential.a
LIB3=${MKLLIB}/libmkl_core.a

LDLIBS += $LIB1 $LIB2 $LIB3 -lpthread -lm -ldl
```

### qt5

```shell
brew install qt
```