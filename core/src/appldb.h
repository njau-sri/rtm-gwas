#ifndef APPLDB_H
#define APPLDB_H

#include "main.h"

class AppLDB
{
public:
    int run(int argc, char *argv[]);

private:
    int perform();

    void save_block_info() const;

    void recode_save_allele(Genotype &gt) const;

    void load_genotype();

    void load_block();

    void make_snpldb(const vector<size_t> &snps, Genotype &gt) const;

    vector< vector<size_t> > index_locus(const vector<string> &chrid) const;

    void find_block(const vector<size_t> &loc, vector<int> &start, vector<int> &stop) const;

    int calc_cild(size_t l1, size_t l2, int &low, int &high) const;

    double test_block(int x, int y, int len, const vector< vector<int> > &cild) const;

private:
    struct Params
    {
        int cut_low_var[5] = {0, 0, 80, 50, 50};
        int max_dist_var[5] = {0, 0, 20000, 30000, 1000000};
        string vcf;
        string ped;
        string hmp;
        string geno;
        string block;
        string out;
        double mhf = 0.01;
        double frac = 0.95;
        int low = 70;
        int high = 98;
        int rec = 90;
        int maxlen = 200000;
    };

    Params m_par;
    Genotype m_gt;
    vector<string> m_chrid;
    vector<int> m_start;
    vector<int> m_stop;
    vector<int> m_length;
    vector<int> m_size;
};

#endif // APPLDB_H
