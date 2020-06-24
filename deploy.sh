#!/bin/bash

# Host          Target
# --------------------------
# CentOS 7      lnx64
# Fedora 30     win32, win64
# macOS 10.14   macos

VER=$(cat VERSION)
PKG=rtm-gwas-$VER-$1

if [[ $# -lt 1 ]]; then
    echo usage: $0 lnx64/macos/win32/win64
    exit 0
fi

if [[ $1 != lnx64 && $1 != macos && $1 != win32 && $1 != win64 ]]; then
    echo ***error: invalid argument: $1
    exit 1
fi

if [[ $1 = macos ]]; then
    export PATH="/usr/local/opt/llvm/bin:$PATH"
    export PATH="/usr/local/opt/qt/bin:$PATH"
fi

make -C rtm-gwas-snpldb -f Makefile.deploy $1 || exit 1
make -C rtm-gwas-gsc -f Makefile.deploy $1 || exit 1
make -C rtm-gwas-assoc -f Makefile.deploy $1 || exit 1

cd rtm-gwas-gui
if [[ $1 = lnx64 ]]; then
    qmake-qt5 || exit 1
    make || exit 1
    strip -s rtm-gwas-gui
elif [[ $1 = "macos" ]]; then
    qmake || exit 1
    make || exit 1
    macdeployqt rtm-gwas-gui.app
elif [[ $1 = win32 ]]; then
    mingw32-qmake-qt4 "CONFIG += static" || exit 1
    make release || exit 1
    mingw-strip -s release/rtm-gwas-gui.exe
elif [[ $1 = win64 ]]; then
    mingw64-qmake-qt4 "CONFIG += static" || exit 1
    make release || exit 1
    mingw-strip -s release/rtm-gwas-gui.exe
fi
cd ..

rm -rf $PKG
mkdir $PKG

if [[ $1 = lnx64 || $1 = macos ]]; then
    mv -v rtm-gwas-snpldb/rtm-gwas-snpldb $PKG
    mv -v rtm-gwas-gsc/rtm-gwas-gsc $PKG
    mv -v rtm-gwas-assoc/rtm-gwas-assoc $PKG
    if [[ $1 = lnx64 ]]; then
        mv -v rtm-gwas-gui/rtm-gwas-gui $PKG
        cp -v rtm-gwas-gui/qt-runtime-linux.txt $PKG
    else
        mv -v rtm-gwas-gui/rtm-gwas-gui.app $PKG
    fi
    tar zcf $PKG.tar.gz $PKG
else
    mv -v rtm-gwas-snpldb/rtm-gwas-snpldb.exe $PKG
    mv -v rtm-gwas-gsc/rtm-gwas-gsc.exe $PKG
    mv -v rtm-gwas-assoc/rtm-gwas-assoc.exe $PKG
    mv -v rtm-gwas-gui/release/rtm-gwas-gui.exe $PKG
    zip -q -r $PKG.zip $PKG
fi
