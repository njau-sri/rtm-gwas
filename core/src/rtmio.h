#ifndef RTMIO_H
#define RTMIO_H

#include "main.h"

int read_file(const string &filename, string &contents);

int write_file( const string &contents, const string &filename);

int read_genotype(const string &filename, Genotype &gt);

int write_genotype(const Genotype &gt, const string &filename);

int read_phenotype(const string &filename, Phenotype &pt);

int write_phenotype(const Phenotype &pt, const string &filename);

int read_square(const string &filename, SquareData &sd);

int write_square(const SquareData &sd, const string &filename);

#endif // RTMIO_H
