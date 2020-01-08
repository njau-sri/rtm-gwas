#!/bin/bash

a=(
install.sh
rtm-gwas-snpldb/src/snpldb.cpp
rtm-gwas-gsc/src/gsc.cpp
rtm-gwas-assoc/src/assoc.cpp
rtm-gwas-gui/src/mainwindow.cpp
)

for e in ${a[@]}
do
    sed -i '0,/'"$1"'/s//'"$2"'/' $e
done

if [[ $(uname) == *"MINGW"* || $(uname) == *"MSYS"* ]]; then
    for e in ${a[@]}
    do
        sed -i 's/$/\r/' $e
    done
fi
