#!/bin/bash

if [[ $# -ne 1 ]]; then
    echo "ERROR: invalid argument"
    echo "./install.sh win32|win64|glnx64|macos"
    exit 1
fi

rm -rf $1

if [[ -f Makefile ]]; then
    make distclean
fi

if [[ $1 = "win32" ]]; then

    mingw32-qmake-qt4 "CONFIG += static" src || exit 1
    make release || exit 1

    mkdir $1
    cp release/rtm-gwas-gui.exe $1/
    mingw-strip -s $1/rtm-gwas-gui.exe

elif [[ $1 = "win64" ]]; then

    mingw64-qmake-qt4 "CONFIG += static" src || exit 1
    make release || exit 1

    mkdir $1
    cp release/rtm-gwas-gui.exe $1/
    mingw-strip -s $1/rtm-gwas-gui.exe

elif [[ $1 = "glnx64" ]]; then

    qmake-qt5 src || exit 1
    make || exit 1

    # sudo apt install libqt5widgets5
    # sudo yum install qt5-qtbase-gui
    mkdir $1
    cp rtm-gwas-gui $1/
    strip -s $1/rtm-gwas-gui

elif [[ $1 == "macos" ]]; then

    export PATH="/usr/local/opt/qt/bin:$PATH"
    export LDFLAGS="-L/usr/local/opt/qt/lib"
    export CPPFLAGS="-I/usr/local/opt/qt/include"

    qmake src || exit 1
    make || exit 1

    mkdir $1
    macdeployqt rtm-gwas-gui.app
    mv rtm-gwas-gui.app $1/

else

    echo "ERROR: invalid argument: $1"
    echo "./install.sh win32|win64|glnx64|macos"
    exit 1

fi
