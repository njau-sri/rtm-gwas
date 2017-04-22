#ifndef VCF_H
#define VCF_H

#include "main.h"

int read_vcf(const string &filename, Genotype &gt);

int write_vcf(const Genotype &gt, const string &filename, bool diploid = true);

#endif // VCF_H
