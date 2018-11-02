#!/bin/bash


BIN=rtm-gwas-v1.5-$1

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
chmod +x install-centos.sh
./install-centos.sh $1

cd ${TOP}/rtm-gwas-gsc
chmod +x install-centos.sh
./install-centos.sh $1

cd ${TOP}/rtm-gwas-assoc
chmod +x install-centos.sh
./install-centos.sh $1

cd ${TOP}/rtm-gwas-gui
chmod +x install-centos.sh
./install-centos.sh $1


cd ${TOP}
rm -rf $BIN
mkdir $BIN

mv ${TOP}/rtm-gwas-snpldb/$1/rtm-gwas-snpldb* ${BIN}/
mv ${TOP}/rtm-gwas-gsc/$1/rtm-gwas-gsc* ${BIN}/
mv ${TOP}/rtm-gwas-assoc/$1/rtm-gwas-assoc* ${BIN}/
mv ${TOP}/rtm-gwas-gui/$1/rtm-gwas-gui* ${BIN}/

tar zcf ${BIN}.tar.gz ${BIN}
