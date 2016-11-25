#ifndef APPGSC_H
#define APPGSC_H

#include "main.h"

class AppGSC
{
public:
    int run(int argc, char *argv[]);

private:
    int perform();

    void load_gentoype();

    void calc_gsc_matrix(vector< vector<double> > &x) const;

    void eigen_decomposition(vector<double> &eval, vector < vector<double> > &evec) const;

private:
    struct Params
    {
        string vcf;
        string ped;
        string hmp;
        string geno;
        string out;
        int top = 10;
    };

    Params m_par;
    Genotype m_gt;
    SquareData m_gsc;
    Phenotype m_evec;
};

#endif // APPGSC_H
