#ifndef PARAMS_H
#define PARAMS_H

#include <QString>

struct Params
{
    static QString exe;
    static QString work_dir;
    static QString open_dir;

    static QString geno;
    static QString pheno;
    static QString covar;
    static QString kin;
    static QString block;

    static int gen_type;
    static int txt_size;
    static int log_size;

    static bool delete_file_onexit;
};

#endif // PARAMS_H
