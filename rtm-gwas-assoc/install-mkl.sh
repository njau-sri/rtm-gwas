#!/bin/bash

rm -rf $1
mkdir $1

if [[ $1 == "glnx64" ]]; then

    make -f Makefile.mkl || exit 1
    cp rtm-gwas-assoc $1/

elif [[ $1 == "macos" ]]; then

    export PATH="/usr/local/opt/llvm@7/bin:$PATH"

    make -f Makefile.mkl macos || exit 1
    cp rtm-gwas-assoc $1/

else

    exit 1

fi
