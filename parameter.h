#ifndef PARAMETER_H
#define PARAMETER_H

#include <QString>

struct Parameter
{
    static QString root;
    static QString work;
    static QString open;

    static QString geno;
    static QString pheno;
    static QString covar;
    static QString block;

    static int txtsize;
    static int logsize;

    static bool delete_onexit;
};

#endif // PARAMETER_H
