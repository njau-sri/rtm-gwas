#ifndef ASSOCLMM_H
#define ASSOCLMM_H

#include "main.h"

void assoc_LMM(bool haploid, const vector< vector<allele_t> > &geno, const vector<size_t> &obs,
               const vector<double> &pheno,
               const vector< vector<double> > &addcov,
               const vector< vector<double> > &intcov,
               const vector< vector<double> > &kin,
               vector<double> &result);

#endif // ASSOCLMM_H
