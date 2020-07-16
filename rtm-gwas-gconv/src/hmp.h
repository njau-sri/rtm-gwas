#ifndef HMP_H
#define HMP_H

#include "vcf.h"

// HapMap Genotype Format
//
//   space (' ') delimited
//
//
//   rs#   alleles  chrom pos  strand  assembly#  center  protLSID  assayLSID  panelLSID  QCcode  ind1  ind2
//   snp1  A/G      1     123  NA      NA         NA      NA        NA         NA         NA      AA    AG
//   snp2  T/C      2     345  NA      NA         NA      NA        NA         NA         NA      TC    TT
//
//
//   rs#        snp identifier
//   alleles    alleles
//   chrom      chromosome identifier
//   pos        position
//   strand     SNP orientation, forward (+), reverse (-)
//   assembly#  reference sequence assembly version
//   center     genotyping center
//   protLSID   protocol identifier
//   assayLSID  assay identifier
//   panelLSID  panel identifier
//   QCcode     quality control

struct HmpEntry
{
    std::string chr;
    std::string id;
    std::vector<std::string> as;
    std::vector<allele_t> gt;
    int pos = 0;
};

int parse_hmp_header(const std::string &s, std::vector<std::string> &v);

int parse_hmp_entry(const std::string &s, HmpEntry &entry);

int read_hmp(const std::string &filename, Genotype &gt);

int write_hmp(const Genotype &gt, const std::string &filename);

#endif // HMP_H
