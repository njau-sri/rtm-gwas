#!/bin/bash


VER=2019.2.dev


export RTM_GWAS_VERSION=$VER

PKG=rtm-gwas-$VER-$1

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
chmod +x install.sh
./install.sh $1 || exit 1

cd $TOP/rtm-gwas-gsc
chmod +x install.sh
./install.sh $1 || exit 1

cd $TOP/rtm-gwas-assoc
if [[ $1 == win* ]]; then
    chmod +x install.sh
    ./install.sh $1 || exit 1
else
    chmod +x install-mkl.sh
    ./install-mkl.sh $1 || exit 1
fi

cd $TOP/rtm-gwas-gui
chmod +x install.sh
./install.sh $1 || exit 1


cd $TOP
rm -rf $PKG
mkdir $PKG

mv $TOP/rtm-gwas-snpldb/$1/rtm-gwas-snpldb* $PKG/
mv $TOP/rtm-gwas-gsc/$1/rtm-gwas-gsc* $PKG/
mv $TOP/rtm-gwas-assoc/$1/rtm-gwas-assoc* $PKG/
mv $TOP/rtm-gwas-gui/$1/rtm-gwas-gui* $PKG/

if [[ $1 == win* ]]; then
    zip -qr $PKG.zip $PKG
else
    tar zcf $PKG.tar.gz $PKG
fi
