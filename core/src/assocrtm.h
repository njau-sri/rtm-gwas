#ifndef ASSOCRTM_H
#define ASSOCRTM_H

#include "main.h"

void assoc_RTM(int multtest, int maxstep, double alpha, double maxrsq,
               bool haploid, const vector< vector<allele_t> > &geno, const vector<size_t> &obs,
               const vector<double> &pheno,
               const vector< vector<double> > &addcov,
               const vector< vector<double> > &intcov,
               vector<size_t> &result, vector<double> &ps);

#endif // ASSOCRTM_H
