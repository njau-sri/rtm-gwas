#ifndef RTMGWASASSOC_H
#define RTMGWASASSOC_H

#include "types.h"
#include "vcf.h"

int assoc_stepwise(const Genotype &gt, const std::vector<isize_t> &gi,
                   const std::vector<double> &y,
                   const std::vector< std::vector<double> > &ac,
                   const std::vector< std::vector<double> > &ic,
                   std::vector<double> &ps, std::vector<isize_t> &in);

int rtm_gwas_assoc(int argc, char *argv[]);

#endif // RTMGWASASSOC_H
