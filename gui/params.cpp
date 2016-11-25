#include "params.h"

QString Params::exe;

QString Params::work_dir;
QString Params::open_dir;

QString Params::geno;
QString Params::pheno;
QString Params::covar;
QString Params::kin;
QString Params::block;

int Params::gen_type = 0;
int Params::txt_size = 10;
int Params::log_size = 1000;

bool Params::delete_file_onexit = false;
