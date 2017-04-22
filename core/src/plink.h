#ifndef PLINK_H
#define PLINK_H

#include "main.h"

bool is_compatible_plink(const Genotype &gt);

int read_plink(const string &prefix, Genotype &gt);

int write_plink(const Genotype &gt, const string &prefix);

#endif // PLINK_H
