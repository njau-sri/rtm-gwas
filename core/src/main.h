#ifndef MAIN_H
#define MAIN_H

#include <cstddef>
#include <string>
#include <vector>

using std::ptrdiff_t;
using std::size_t;

using std::string;
using std::vector;

using allele_t = unsigned char;

struct Genotype
{
    vector<string> ind;
    vector<string> loc;
    vector<string> chr;
    vector<int> pos;
    vector<double> dist;
    vector< vector<allele_t> > dat;
    vector< vector<string> > allele;
    int ploidy = 0;
};

struct Phenotype
{
    vector<string> ind;
    vector<string> phe;
    vector<string> env;
    vector<string> blk;
    vector< vector<double> > dat;
};

struct SquareData
{
    vector<string> ind;
    vector< vector<double> > dat;
};

#endif // MAIN_H
