#include "parameter.h"

QString Parameter::root;
QString Parameter::work;
QString Parameter::open;

QString Parameter::vcf;
QString Parameter::pheno;
QString Parameter::covar;
QString Parameter::block;
QString Parameter::grm;

int Parameter::txtsize = 10;
int Parameter::logsize = 1000;

int Parameter::openmp = 0;

bool Parameter::delete_onexit = false;
