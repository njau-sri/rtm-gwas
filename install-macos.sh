#!/bin/bash


BIN=rtm-gwas-v1.5.1-macos

TOP=$(pwd)


rm -rf rtm-gwas-snpldb
git clone https://github.com/njau-sri/rtm-gwas-snpldb.git

rm -rf rtm-gwas-gsc
git clone https://github.com/njau-sri/rtm-gwas-gsc.git

rm -rf rtm-gwas-assoc
git clone https://github.com/njau-sri/rtm-gwas-assoc.git

rm -rf rtm-gwas-gui
git clone https://github.com/njau-sri/rtm-gwas-gui.git


cd ${TOP}/rtm-gwas-snpldb
chmod +x install-*.sh
./install-macos.sh

cd ${TOP}/rtm-gwas-gsc
chmod +x install-*.sh
./install-macos.sh

cd ${TOP}/rtm-gwas-assoc
chmod +x install-*.sh
./install-macos-mkl.sh

cd ${TOP}/rtm-gwas-gui
chmod +x install-*.sh
./install-macos.sh


APP=rtm-gwas-gui.app


cd ${TOP}
rm -rf $BIN
mkdir $BIN

mv ${TOP}/rtm-gwas-gui/macos/${APP} ${BIN}/
mv ${TOP}/rtm-gwas-snpldb/macos/rtm-gwas-snpldb ${BIN}/${APP}/Contents/MacOS/
mv ${TOP}/rtm-gwas-gsc/macos/rtm-gwas-gsc ${BIN}/${APP}/Contents/MacOS/
mv ${TOP}/rtm-gwas-assoc/macos/rtm-gwas-assoc ${BIN}/${APP}/Contents/MacOS/

tar zcf ${BIN}.tar.gz ${BIN}
