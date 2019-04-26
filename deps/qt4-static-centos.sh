#!/bin/bash

tar zxf qt-everywhere-opensource-src-4.8.7.tar.gz

mkdir qt4-static && cd "$_"

../qt-everywhere-opensource-src-4.8.7/configure -prefix /opt/qt4-static -release -opensource -static \
-no-qt3support -no-xmlpatterns -no-multimedia -no-audio-backend -no-phonon -no-phonon-backend \
-no-svg -no-webkit -no-javascript-jit -no-script -no-scripttools -no-declarative -no-declarative-debug \
-qt-zlib -no-gif -no-libtiff -no-libpng -no-libmng -no-libjpeg -no-openssl \
-nomake tools -nomake examples -nomake demos -nomake docs -nomake translations \
-no-pch
