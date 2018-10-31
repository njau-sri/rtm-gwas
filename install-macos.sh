#!/bin/bash

export PATH="/usr/local/opt/qt/bin:$PATH"
export LDFLAGS="-L/usr/local/opt/qt/lib"
export CPPFLAGS="-I/usr/local/opt/qt/include"

rm -rf macos
mkdir macos

make distclean

qmake
make

macdeployqt rtm-gwas.app
cp -r rtm-gwas.app macos/
