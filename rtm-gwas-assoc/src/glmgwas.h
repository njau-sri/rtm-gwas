#ifndef GLMGWAS_H
#define GLMGWAS_H

#include "types.h"
#include "vcf.h"
#include "pheno.h"

void intersect(Genotype &gt, Phenotype &pt, Covariate &ct, std::vector<isize_t> &gi);

// design matrix for _ENV_ and _BLK_ factor in RCBD
void design(const Phenotype &pt, std::vector< std::vector<double> > &ac, std::vector< std::vector<double> > &ic);

// out   model: p r2 / g: p r2 / ge: p r2
int assoc_glm(const Genotype &gt, const std::vector<isize_t> &gi,
              const std::vector<double> &y,
              const std::vector< std::vector<double> > &ac,
              const std::vector< std::vector<double> > &ic,
              std::vector< std::vector<double> > &out);

int assoc_glm_omp(const Genotype &gt, const std::vector<isize_t> &gi,
                  const std::vector<double> &y,
                  const std::vector< std::vector<double> > &ac,
                  const std::vector< std::vector<double> > &ic,
                  std::vector< std::vector<double> > &out);

int glm_gwas(int argc, char *argv[]);

#endif // GLMGWAS_H
