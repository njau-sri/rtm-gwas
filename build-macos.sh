#!/bin/bash

# brew install llvm@9 libomp boost qt
export PATH="/usr/local/opt/llvm@9/bin:$PATH"
export PATH="/usr/local/opt/qt/bin:$PATH"

if [[ $# -ne 0 ]]; then
    make -C rtm-gwas-gconv -f Makefile.macos clean
    make -C rtm-gwas-snpldb -f Makefile.macos clean
    make -C rtm-gwas-gsc -f Makefile.macos clean
    make -C rtm-gwas-assoc -f Makefile.macos clean
    make -C rtm-gwas-ld -f Makefile.macos clean
    make -C rtm-gwas-gui distclean
    exit 0
fi

mkdir bin

make -C rtm-gwas-gconv -f Makefile.macos || exit 1
mv -v rtm-gwas-gconv/rtm-gwas-gconv bin

make -C rtm-gwas-snpldb -f Makefile.macos || exit 1
mv -v rtm-gwas-snpldb/rtm-gwas-snpldb bin

make -C rtm-gwas-gsc -f Makefile.macos || exit 1
mv -v rtm-gwas-gsc/rtm-gwas-gsc bin

make -C rtm-gwas-assoc -f Makefile.macos || exit 1
mv -v rtm-gwas-assoc/rtm-gwas-assoc bin

make -C rtm-gwas-ld -f Makefile.macos || exit 1
mv -v rtm-gwas-ld/rtm-gwas-ld bin

cd rtm-gwas-gui
qmake || exit 1
cd ..

make -C rtm-gwas-gui || exit 1
macdeployqt rtm-gwas-gui/rtm-gwas-gui.app
mv -v rtm-gwas-gui/rtm-gwas-gui.app bin
