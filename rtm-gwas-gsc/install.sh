#!/bin/bash

rm -rf $1
mkdir $1

make clean

if [[ $1 == "glnx64" ]]; then

    make || exit 1
    cp rtm-gwas-gsc $1/

elif [[ $1 == "win32" ]]; then

    make win32 || exit 1
    cp rtm-gwas-gsc.exe $1/

elif [[ $1 == "win64" ]]; then

    make win64 || exit 1
    cp rtm-gwas-gsc.exe $1/

elif [[ $1 == "macos" ]]; then

    export PATH="/usr/local/opt/llvm@7/bin:$PATH"
    make macos || exit 1
    cp rtm-gwas-gsc $1/

else

    exit 1

fi
