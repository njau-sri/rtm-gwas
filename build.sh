#!/bin/bash

VER=v1.1

if [ $# -ne 1 ]; then
    echo "usage: $0 lnx64|osx64|win32|win64"
    exit 1
fi

TOP=`pwd`

rm -rf build
mkdir -p build/rtm-gwas

cd $TOP/core/deps
if [ ! -d Rmath-$1 ]; then
    chmod +x build-Rmath.sh
    ./build-Rmath.sh $1 || exit 1
fi
if [ ! -d OpenBLAS-$1 ]; then
    chmod +x build-OpenBLAS.sh
    ./build-OpenBLAS.sh $1 || exit 1
fi

cd $TOP/core
make -f Makefile.lnx64 clean
make -f Makefile.$1 || exit 1
cp rtm-gwas* $TOP/build/rtm-gwas

cd $TOP/gui
make distclean
if [ $1 == "lnx64" ]; then
    qmake-qt4 -spec macx-g++
    make || exit 1
	strip -s rtm-gwas-gui
    cp rtm-gwas-gui $TOP/build/rtm-gwas
elif [ $1 == "osx64" ]; then
    qmake-4.8 -spec macx-g++
    sed -i.bak 's/ -g -gdwarf-2 / -O2 /' Makefile
    make || exit 1
    macdeployqt-4.8 rtm-gwas-gui.app
    mv $TOP/build/rtm-gwas/rtm-gwas rtm-gwas-gui.app/Contents/MacOS
    rmdir $TOP/build/rtm-gwas
    mv rtm-gwas-gui.app $TOP/build/rtm-gwas.app
elif [ $1 == "win32" ]; then
    i686-w64-mingw32-qmake-qt4
    make || exit 1
	i686-w64-mingw32-strip -s release/rtm-gwas-gui.exe
    cp release/rtm-gwas-gui.exe $TOP/build/rtm-gwas
elif [ $1 == "win64" ]; then
    x86_64-w64-mingw32-qmake-qt4
    make || exit 1
	x86_64-w64-mingw32-strip -s release/rtm-gwas-gui.exe
    cp release/rtm-gwas-gui.exe $TOP/build/rtm-gwas
fi

cd $TOP/build

if [ $1 == "osx64" ]; then
sed -i.bak '4s/<dict>/<dict>\
	<key>CFBundleDisplayName<\/key>\
	<string>RTM-GWAS<\/string>/' rtm-gwas.app/Contents/Info.plist
fi

if [ $1 == "win32" ] || [ $1 == "win64" ]; then
    zip -r rtm-gwas-$VER-$1.zip rtm-gwas
    mv rtm-gwas-$VER-$1.zip $TOP
else
    if [ $1 == "osx64" ]; then
        tar zcf rtm-gwas-$VER-$1.tar.gz rtm-gwas.app
    else
        tar zcf rtm-gwas-$VER-$1.tar.gz rtm-gwas
    fi
    mv rtm-gwas-$VER-$1.tar.gz $TOP
fi
