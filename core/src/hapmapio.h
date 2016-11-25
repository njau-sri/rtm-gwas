#ifndef HAPMAPIO_H
#define HAPMAPIO_H

#include "main.h"

bool is_compatible_hapmap(const Genotype &gt);

int read_hapmap(const string &filename, Genotype &gt);

int write_hapmap(const Genotype &gt, const string &filename);

#endif // HAPMAPIO_H
