#!/bin/bash


BIN=rtm-gwas-2019.0.dev-$1

TOP=$(pwd)


rm -rf rtm-gwas-snpldb
git clone https://github.com/njau-sri/rtm-gwas-snpldb.git

rm -rf rtm-gwas-gsc
git clone https://github.com/njau-sri/rtm-gwas-gsc.git

rm -rf rtm-gwas-assoc
git clone https://github.com/njau-sri/rtm-gwas-assoc.git

rm -rf rtm-gwas-gui
git clone https://github.com/njau-sri/rtm-gwas-gui.git


cd $TOP/rtm-gwas-snpldb
chmod +x install-cross.sh
./install-cross.sh $1 || exit 1

cd $TOP/rtm-gwas-gsc
chmod +x install-cross.sh
./install-cross.sh $1 || exit 1

cd $TOP/rtm-gwas-assoc
if [ $1 == "glnx64" ]; then
    chmod +x install-glnx64-mkl.sh
    ./install-glnx64-mkl.sh $1 || exit 1
else
    chmod +x install-cross-lapack.sh
    ./install-cross-lapack.sh $1 || exit 1
fi

cd $TOP/rtm-gwas-gui
chmod +x install-cross.sh
./install-cross.sh $1 || exit 1


cd $TOP
rm -rf $BIN
mkdir $BIN

mv $TOP/rtm-gwas-snpldb/$1/rtm-gwas-snpldb* $BIN/
mv $TOP/rtm-gwas-gsc/$1/rtm-gwas-gsc* $BIN/
mv $TOP/rtm-gwas-assoc/$1/rtm-gwas-assoc* $BIN/
mv $TOP/rtm-gwas-gui/$1/rtm-gwas-gui* $BIN/

if [ $1 == "glnx64" ]; then
    tar zcf $BIN.tar.gz $BIN
else
    zip -qr $BIN.zip $BIN
fi
