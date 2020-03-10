#!/bin/bash

if [[ $# < 2 ]]; then
    echo "usage: ./ver.sh curr_version new_version <unix>"
    exit 1
fi

a=(
install.sh
rtm-gwas-snpldb/src/version.h
rtm-gwas-gsc/src/version.h
rtm-gwas-assoc/src/version.h
rtm-gwas-gui/src/version.h
)

for e in ${a[@]}
do
    sed -i '0,/'"$1"'/s//'"$2"'/' $e
    if [[ $# < 3 ]]; then
        sed -i 's/$/\r/' $e
    fi
done
