#!/bin/bash


BIN=rtm-gwas-1.5-macos

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
chmod +x install-macos.sh
./install-macos.sh

cd ${TOP}/rtm-gwas-gsc
chmod +x install-macos.sh
./install-macos.sh

cd ${TOP}/rtm-gwas-assoc
chmod +x install-macos.sh
./install-macos.sh

cd ${TOP}/rtm-gwas-gui
chmod +x install-macos.sh
./install-macos.sh


cd ${TOP}
rm -rf $BIN
mkdir $BIN

mv ${TOP}/rtm-gwas-snpldb/macos/rtm-gwas-snpldb* ${BIN}/
mv ${TOP}/rtm-gwas-gsc/macos/rtm-gwas-gsc* ${BIN}/
mv ${TOP}/rtm-gwas-assoc/macos/rtm-gwas-assoc* ${BIN}/
mv ${TOP}/rtm-gwas-gui/macos/rtm-gwas-gui* ${BIN}/

tar zcf ${BIN}.tar.gz ${BIN}
