#ifndef HAPLOTYPE_H
#define HAPLOTYPE_H


#include "vcf.h"


// sidx   global index of snps to be grouped


std::size_t infer_haplotype(double maf, const Genotype &gt, const std::vector<std::size_t> &sidx, Genotype &ht);


// igt   individual genotype
// pgt   parent genotype
// fam   ril family index


std::size_t infer_haplotype_ril(const Genotype &igt, const Genotype &pgt,
                                const std::vector<std::size_t> &sidx,
                                const std::vector< std::vector<std::size_t> > &fam,
                                Genotype &iht, Genotype &pht);


#endif // HAPLOTYPE_H
