#ifndef APPDATA_H
#define APPDATA_H

#include "main.h"

class AppData
{
public:
    int run(int argc, char *argv[]);

private:
    int perform();

    void load_genotype();

    void apply_filter();

    void save_genotype() const;

    void filter_ind();

    void filter_mis_ind();

    void filter_loc();

    void filter_allele_number();

    void filter_mis_loc();

    void filter_maf();

    void filter_het();

private:
    struct Params
    {
        string vcf;
        string ped;
        string hmp;
        string geno;
        string format;
        string out;
        string loc;
        string ind;

        double maf = 0.0;
        double het = 1.0;
        double mloc = 1.0;
        double mind = 1.0;

        int min_alleles = 0;
        int max_alleles = 255;

        bool diploid = false;
    };

    Params m_par;
    Genotype m_gt;
    vector<bool> m_keep_ind;
    vector<bool> m_keep_loc;
};

#endif // APPDATA_H
