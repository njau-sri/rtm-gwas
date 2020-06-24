#ifndef VCF_H
#define VCF_H

#include <string>
#include <vector>

// Variant Call Format - VCF Genotype Format
//
//   Tab-delimited
//
//   ##fileformat=VCFv4.2
//   ##CHROM  POS  ID    REF  ALT  QUAL FILTER INFO FORMAT  ind1  ind2
//   1        123  snp1  A    G    .    .      .    GT      0/0   0/1
//   2        456  snp2  T    C    .    .      .    GT      0/1   0/0
//
//   CHROM   chromosome identifier
//   POS     reference position
//   ID      variant identifier
//   REF     reference base(s)
//   ALT     alternate base(s)
//   QUAL    quality score
//   FILTER  filter status
//   INFO    additional information
//   FORMAT  genotype data format
//
//   http://samtools.github.io/hts-specs/

using allele_t = unsigned char;

struct VcfEntry
{
    std::string chr;
    std::string id;
    std::vector<std::string> as;
    std::vector<allele_t> gt;
    int pos = 0;
    int ploidy = 0;
};

struct Genotype
{
    std::vector<std::string> ind;
    std::vector<std::string> loc;
    std::vector<std::string> chr;
    std::vector<int> pos;  // chromosome length <=2Gb
    std::vector< std::vector<allele_t> > dat;
    std::vector< std::vector<std::string> > allele;
    int ploidy = 0;
    enum class Source { UNKNOWN, VCF, PED, HAPMAP, GENO, SNPLDB };
    Source source = Source::UNKNOWN;
};

int parse_vcf_header(const std::string &s, std::vector<std::string> &v);

int parse_vcf_entry(const std::string &s, VcfEntry &entry);

int read_vcf(const std::string &filename, Genotype &gt);

int read_vcf_gz(const std::string &filename, Genotype &gt);

int write_vcf(const Genotype &gt, const std::string &filename, bool force_diploid = true);

int write_vcf_gz(const Genotype &gt, const std::string &filename, bool force_diploid = true);

#endif // VCF_H
