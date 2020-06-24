#ifndef HAPLOTYPE_H
#define HAPLOTYPE_H

#include "vcf.h"
#include "types.h"

struct Haplotype
{
    struct Parameter
    {
        double maf = 0.01;
        double identity = 0.0;
    };
    Parameter par;

    std::string chr;
    std::vector<int> dat;
    std::vector<int> pdat;
    std::vector<std::string> hap;
    int start = 0;
    int stop = 0;
    int size = 0; // number of snps
    int filter = 0;
};

int form_haplotype(const Genotype &gt, const std::vector<isize_t> &sidx, Haplotype *out);

// RIL, p1 p2 i1 i2 ...
int form_haplotype_ril(const Genotype &gt, const std::vector<isize_t> &sidx, Haplotype *out);

// NAM, pc p1 p2 p3 f1-i1 f1-i2 ... f2-i1 f2-i2 ... f3-i1 f3-i2 ...
// fam, number of lines in each family
int form_haplotype_nam(const Genotype &gt, const std::vector<isize_t> &sidx, const std::vector<isize_t> &fam,
                        Haplotype *out);

#endif // HAPLOTYPE_H
