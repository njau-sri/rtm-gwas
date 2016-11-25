#ifndef ASSOCLM_H
#define ASSOCLM_H

#include "main.h"

void assoc_LM(bool haploid, const vector< vector<allele_t> > &geno,  const vector<size_t> &obs,
              const vector<double> &pheno,
              const vector< vector<double> > &addcov,
              const vector< vector<double> > &intcov,
              vector<double> &result);

#endif // ASSOCLM_H
