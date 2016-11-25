#ifndef APPASSOC_H
#define APPASSOC_H

#include "main.h"

class AppAssoc
{
public:
    int run(int argc, char *argv[]);

private:
    int perform();

    void load_genotype();

    void load_phenotype();

    void load_covariate();

    void load_kinship();

    void merge_data();

    void parse_factor();

    void LM();

    void LMM();

    void RTM();

    string fit(int phe, const vector<size_t> &loc) const;

private:
    struct Params
    {
        string vcf;
        string ped;
        string hmp;
        string geno;
        string pheno;
        string covar;
        string kin;
        string out;
        string method;
        string multtest;
        double maxrsq = 0.99;
        double alpha = 0.05;
        double preselect = 0.05;
        int maxqtl = 0;
        int sstype = 1;
        bool no_gxe = false;
    };

    Params m_par;
    Genotype m_gt;
    Phenotype m_pt;
    Phenotype m_ct;
    SquareData m_kin;
    vector<size_t> m_obs;
    vector< vector<double> > m_addcov;
    vector< vector<double> > m_intcov;
};

#endif // APPASSOC_H
