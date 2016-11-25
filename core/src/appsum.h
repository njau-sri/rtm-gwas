#ifndef APPSUM_H
#define APPSUM_H

#include "main.h"

class AppSum
{
public:
    int run(int argc, char *argv[]);

private:
    int perform();

    void load_genotype();

    void calc_save_stats_locus() const;

    void calc_save_stats_indiv() const;

    void calc_save_stats_allele() const;

private:
    struct Params
    {
        string vcf;
        string ped;
        string hmp;
        string geno;
        string out;
        bool allele = false;
        bool indiv = false;
    };

    Params m_par;
    Genotype m_gt;
};

#endif // APPSUM_H
