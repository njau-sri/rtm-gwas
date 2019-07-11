#ifndef BLOCK_GABRIEL_H
#define BLOCK_GABRIEL_H


#include <vector>
#include <utility>


struct BlockGabriel
{
    double frac = 0.95;
    int maxlen = 100000;
    int llim = 70;
    int ulim = 98;
    int recomb = 90;
    int batch = 10000;
};


// pos  snp position
// dat  snp genotype, 0,1,2,3 -> AA,Aa,aa,NN


int find_block(const BlockGabriel &par, const std::vector<int> &pos, const std::vector< std::vector<char> > &dat,
               std::vector< std::pair<int,int> > &bpos);


int find_block_omp(const BlockGabriel &par, const std::vector<int> &pos, const std::vector< std::vector<char> > &dat,
                   std::vector< std::pair<int,int> > &bpos);


#endif // BLOCK_GABRIEL_H
