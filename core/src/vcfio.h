#ifndef VCFIO_H
#define VCFIO_H

#include "main.h"

int read_vcf(const string &filename, Genotype &gt);

int write_vcf(const Genotype &gt, const string &filename, bool force_diploid = false);

#endif // VCFIO_H
