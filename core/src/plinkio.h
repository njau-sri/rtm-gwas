#ifndef PLINKIO_H
#define PLINKIO_H

#include "main.h"

bool is_compatible_plink(const Genotype &gt);

int read_plink(const string &prefix, Genotype &gt);

int write_plink(const Genotype &gt, const string &prefix);

#endif // PLINKIO_H
