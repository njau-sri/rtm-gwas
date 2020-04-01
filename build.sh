#!/bin/bash

MAKE=make
QMAKE=qmake-qt5

if [[ $# -ne 0 ]]; then
    $MAKE -C rtm-gwas-snpldb clean
    $MAKE -C rtm-gwas-gsc clean
    $MAKE -C rtm-gwas-assoc clean
    $MAKE -C rtm-gwas-gui distclean
    rm -rf rtm-gwas-gui/Makefile*
    rm -rf rtm-gwas-gui/object_*
    rm -rf rtm-gwas-gui/debug
    rm -rf rtm-gwas-gui/release
    exit 0
fi

mkdir bin

$MAKE -C rtm-gwas-snpldb || exit 1
mv -v rtm-gwas-snpldb/rtm-gwas-snpldb bin

$MAKE -C rtm-gwas-gsc || exit 1
mv -v rtm-gwas-gsc/rtm-gwas-gsc bin

$MAKE -C rtm-gwas-assoc || exit 1
mv -v rtm-gwas-assoc/rtm-gwas-assoc bin

if [[ $(uname -o) == *Linux* ]]; then
    if [[ $(cat /etc/*-release) == *Ubuntu* ]]; then
        QMAKE=/usr/lib/qt5/bin/qmake
    fi
elif [[ $(uname -o) == *Msys* ]]; then
    QMAKE=qmake
fi

cd rtm-gwas-gui
$QMAKE src || exit 1
if [[ -f Makefile.Release ]]; then
    $MAKE release || exit 1
    mv -v release/rtm-gwas-gui ../bin
else
    $MAKE -C rtm-gwas-gui || exit 1
    mv -v rtm-gwas-gui ../bin
fi
