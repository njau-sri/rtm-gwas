#!/bin/bash

export PATH="/usr/local/opt/qt/bin:$PATH"
export LDFLAGS="-L/usr/local/opt/qt/lib"
export CPPFLAGS="-I/usr/local/opt/qt/include"

rm -rf rtm-gwas-mac
mkdir rtm-gwas-mac

make distclean

qmake
make

macdeployqt rtm-gwas.app
mv rtm-gwas.app rtm-gwas-mac
