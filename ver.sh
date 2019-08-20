#!/bin/bash

sed -i -e 's/'"$1"'/'"$2"'/' \
install.sh \
rtm-gwas-snpldb/src/snpldb.cpp \
rtm-gwas-gsc/src/gsc.cpp \
rtm-gwas-assoc/src/assoc.cpp \
rtm-gwas-gui/src/mainwindow.cpp
