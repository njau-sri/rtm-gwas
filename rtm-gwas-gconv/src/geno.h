#ifndef GENO_H
#define GENO_H

#include "vcf.h"

// General Genotype File Format
//
//   white-space (' ' or '\t') delimited
//
//   Locus  Chromosome  Position   ind1      ind2
//   Mk1    1           100        630/630   909/909
//   Mk2    1           200        557/711   711/711
//   Mk3    1           300        445/445   445/668
//   Mk4    2           1001       307/307   340/340
//   Mk5    2           1002       264/273   264/264
//
//   - allele code:       number or character or string
//   - allele separator:  '/' or ':' or WHITESPACE
//   - missing genotype:  '?' or '.' or '-', or 'N'

int read_geno(const std::string &filename, Genotype &gt);

int write_geno(const Genotype &gt, const std::string &filename);

#endif // GENO_H
