#ifndef HAPBLOCK_H
#define HAPBLOCK_H

#include "types.h"

struct HapBlockGabriel
{
    struct Parameter
    {
        double tol = 1e-5;
        double frac = 0.95;
        int maxlen = 100000;
        int llim = 70;
        int ulim = 98;
        int recomb = 90;
        int batch = 10000;
        int thread = 0;
        int maxit = 1000;
    };
    Parameter par;

    std::vector< std::pair<int, int> > block;
};

// Gabriel et al., Science, 2002, 296:2225-2229. DOI: 10.1126/science.1069424
// Barrett et al., Bioinformatics, 2005, 21:263-265. DOI: 10.1093/bioinformatics/bth457
// pos  SNP position
// geno SNP genotype, 0:AA, 1:Aa, 2:aa, -1:missing

int find_hapblock_gabriel(const std::vector<int> &pos, const std::vector< std::vector<int8_t> > &geno, HapBlockGabriel *out);

#endif // HAPBLOCK_H
