#ifndef PHENO_H
#define PHENO_H


#include <string>
#include <vector>


struct Phenotype
{
    std::vector<std::string> ind;
    std::vector<std::string> phe;
    std::vector<std::string> env;
    std::vector<std::string> blk;
    std::vector< std::vector<double> > dat;
};


struct Covariate
{
    std::vector<std::string> ind;
    std::vector<std::string> phe;
    std::vector< std::vector<double> > dat;
};


struct SquareData
{
    std::vector<std::string> ind;
    std::vector< std::vector<double> > dat;
};


int read_pheno(const std::string &filename, Phenotype &pt);

int write_pheno(const Phenotype &pt, const std::string &filename);

int read_covar(const std::string &filename, Covariate &ct);

int write_covar(const Covariate &ct, const std::string &filename);

int read_square(const std::string &filename, SquareData &sd);

int write_square(const SquareData &sd, const std::string &filename);


#endif // PHENO_H
