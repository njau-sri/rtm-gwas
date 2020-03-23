#!/bin/bash

MY_MAKE=make
MY_QMAKE=qmake-qt5

if [[ $# -ne 0 ]]; then
    $MY_MAKE -C rtm-gwas-snpldb clean
    $MY_MAKE -C rtm-gwas-gsc clean
    $MY_MAKE -C rtm-gwas-assoc clean
    $MY_MAKE -C rtm-gwas-gui distclean
    exit 0
fi

mkdir bin

$MY_MAKE -C rtm-gwas-snpldb || exit 1
mv -v rtm-gwas-snpldb/rtm-gwas-snpldb bin

$MY_MAKE -C rtm-gwas-gsc || exit 1
mv -v rtm-gwas-gsc/rtm-gwas-gsc bin

$MY_MAKE -C rtm-gwas-assoc || exit 1
mv -v rtm-gwas-assoc/rtm-gwas-assoc bin

if [[ $(uname -o) == *Linux* ]]; then
    if [[ $(cat /etc/*-release) == *Ubuntu* ]]; then
        MY_QMAKE=/usr/lib/qt5/bin/qmake
    fi
elif [[ $(uname -o) == *Msys* ]]; then
    MY_QMAKE=qmake
fi

cd rtm-gwas-gui
$MY_QMAKE src || exit 1
cd ..

$MY_MAKE -C rtm-gwas-gui || exit 1
mv -v rtm-gwas-gui/rtm-gwas-gui bin
