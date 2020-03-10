#!/bin/bash

if [[ $# -ne 1 ]]; then
    echo "ERROR: invalid argument"
    echo "./install.sh win32|win64|glnx64|macos"
    exit 1
fi

rm -rf $1
mkdir $1

if [[ $1 = "win32" ]]; then

    # fedora 30, dnf install make mingw32-qt-static

    make distclean
    mingw32-qmake-qt4 "CONFIG += static" src || exit 1
    make release || exit 1
    cp release/rtm-gwas-gui.exe $1/
    mingw-strip -s $1/rtm-gwas-gui.exe
    make distclean

elif [[ $1 = "win64" ]]; then

    # fedora 30, dnf install make mingw64-qt-static

    make distclean
    mingw64-qmake-qt4 "CONFIG += static" src || exit 1
    make release || exit 1
    cp release/rtm-gwas-gui.exe $1/
    mingw-strip -s $1/rtm-gwas-gui.exe
    make distclean

elif [[ $1 = "glnx64" ]]; then

    # centos 7, yum install qt5-qtbase-devel

    make distclean
    qmake-qt5 src || exit 1
    make || exit 1
    # sudo apt install libqt5widgets5
    # sudo yum install qt5-qtbase-gui
    cp rtm-gwas-gui $1/
    strip -s $1/rtm-gwas-gui
    make distclean

elif [[ $1 == "macos" ]]; then

    # macOS Mojave, brew install qt

    export PATH="/usr/local/opt/qt/bin:$PATH"

    make distclean
    qmake src || exit 1
    make || exit 1
    macdeployqt rtm-gwas-gui.app
    mv rtm-gwas-gui.app $1/
    make distclean

else

    echo "ERROR: invalid argument: $1"
    echo "./install.sh win32|win64|glnx64|macos"
    exit 1

fi
