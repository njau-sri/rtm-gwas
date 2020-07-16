#ifndef PED_H
#define PED_H

#include "vcf.h"

// PLINK PED/MAP Genotype Format
//
//   white-space (' ' or '\t') delimited
//
//   PED File
//     FAM001  1  0 0  1  2  A A  G G  A C
//     FAM001  2  0 0  1  2  A A  A G  0 0
//
//   The first 6 columns are mandatory:
//     Family ID
//     Individual ID
//     Paternal ID
//     Maternal ID
//     Sex (1=male; 2=female; other=unknown)
//     Phenotype (1=unaffected; 2=affected; 0=missing disease phenotype; -9=missing quantitative phenotype)
//
//   MAP File
//     1  rs123456  0  1234555
//     1  rs234567  0  1237793
//     1  rs224534  0  1237697
//     1  rs233556  0  1337456
//
//   Exactly 4 columns:
//     chromosome (1-22, X, Y or 0 if unplaced)
//     rs# or snp identifier
//     Genetic distance (morgans)
//     Base-pair position (bp units)
//
// http://zzz.bwh.harvard.edu/plink/data.shtml

struct PedEntry
{
    std::string fid;
    std::string iid;
    std::string pid;
    std::string mid;
    std::vector<allele_t> gt;
    int sex = 1;
    int pheno = 0;
};

struct MapEntry
{
    std::string chr;
    std::string id;
    double dist = 0.0;
    int pos = 0;
};

int parse_ped_entry(const std::string &s, PedEntry &entry);

int parse_map_entry(const std::string &s, MapEntry &entry);

int read_ped(const std::string &filename, Genotype &gt);

int write_ped(const Genotype &gt, const std::string &filename);

#endif // PED_H
