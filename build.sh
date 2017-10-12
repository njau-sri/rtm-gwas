#!/bin/bash

rm -rf rtm-gwas-$1
mkdir rtm-gwas-$1

make distclean

if [ $1 == "lnx64" ]; then

    qmake-qt4
    make

    strip -s rtm-gwas
    mv rtm-gwas rtm-gwas-$1/

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-qmake-qt4
    make

    i686-w64-mingw32-strip -s release/rtm-gwas.exe
    mv release/rtm-gwas.exe rtm-gwas-$1/

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-qmake-qt4
    make

    x86_64-w64-mingw32-strip -s release/rtm-gwas.exe
    mv release/rtm-gwas.exe rtm-gwas-$1/

fi
