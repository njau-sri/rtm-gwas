#include <fstream>
#include <iostream>
#include <algorithm>
#include "vcf.h"
#include "split.h"


using std::size_t;


namespace {


// parse GT, return ploidy number (1 or 2), otherwise error
int parse_vcf_gt(const char *s, size_t n, int &a, int &b)
{
    auto beg = s;
    auto end = std::find(s, s + n, ':');
    auto len = static_cast<size_t>(end - beg);

    if (len == 0)
        return -1;

    if (len == 1) {
        if (*beg == '.')
            return 1;
        a = *beg - '0';
        return a < 0 || a > 9 ? -1 : 1;
    }

    if (len == 3 && (beg[1] == '/' || beg[1] == '|')) {
        if (beg[0] != '.') {
            a = beg[0] - '0';
            if (a < 0 || a > 9)
                return -1;
        }

        if (beg[2] != '.') {
            b = beg[2] - '0';
            if (b < 0 || b > 9)
                return -1;
        }

        return 2;
    }

    if (std::count(beg, end, '/') + std::count(beg, end, '|') > 1) {
        std::cerr << "ERROR: polyploidy genotype is not supported: " << std::string(beg,end) << "\n";
        return -1;
    }

    size_t pos = 0;
    while (pos < len && beg[pos] != '/' && beg[pos] != '|')
        ++pos;
    if (pos == 0 || pos == len - 1)
        return -1;

    if (pos == len) {
        a = std::stoi(std::string(beg,end));
        return a < 0 ? -1 : 1;
    }

    std::string gt(beg, beg + pos);
    if (gt != ".") {
        a = std::stoi(gt);
        if (a < 0)
            return -1;
    }

    gt.assign(beg + pos + 1, end);
    if (gt != ".") {
        b = std::stoi(gt);
        if (b < 0)
            return -1;
    }

    return 2;
}


} // namespace


int parse_vcf_header(const std::string &s, std::vector<std::string> &v)
{
    v.clear();
    split(s, "\t", v);

    if (v.size() != 8 && v.size() < 10) {
        std::cerr << "ERROR: incorrect number of columns in VCF header line: " << v.size() << "\n";
        return 1;
    }

    static const char * const z[8] = { "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO" };

    for (size_t i = 0; i < 8; ++i) {
        if (v[i] != z[i]) {
            std::cerr << "ERROR: incorrect column name in header line:" << v[i] << "\n";
            return 1;
        }
    }

    if (v.size() > 9) {
        if (v[8] != "FORMAT") {
            std::cerr << "ERROR: FORMAT is required at 9th field in VCF header line: " << v[8] << "\n";
            return 1;
        }
        v.erase(v.begin(), v.begin() + 9);
    }
    else
        v.clear();

    return 0;
}

int parse_vcf_entry(const std::string &s, VcfEntry &e)
{
    std::vector<Token> v;
    split(s, "\t", v);

    auto n = v.size();

    if (n != 8 && n < 10) {
        std::cerr << "ERROR: incorrect number of columns at VCF entry line: " << n << "\n";
        return 1;
    }

    e.chr = v[0].to_string();
    e.pos = std::stoi(v[1].to_string());

    if (v[2].size() == 1 && v[2][0] == '.')
        e.id = v[0].to_string() + "_" + v[1].to_string();
    else
        e.id = v[2].to_string();

    e.as.clear();
    e.as.push_back(v[3].to_string());
    split(v[4].to_string(), ",", e.as);

    if (n == 8)
        return 0;

    if (v[8].size() < 2 || v[8][0] != 'G' || v[8][1] != 'T') {
        std::cerr << "ERROR: GT is required and must be the first sub-field of FORMAT: " << v[8].to_string() << "\n";
        return 1;
    }

    auto na = static_cast<int>(e.as.size());
    if (na > 255) {
        std::cerr << "ERROR: exceed the maximum number of alleles (255): " << na << "\n";
        return 1;
    }

    e.gt.clear();

    for (size_t i = 9; i < n; ++i) {
        int a = -9, b = -9;
        int info = parse_vcf_gt(v[i].data(), v[i].size(), a, b);

        if ((info != 1 && info != 2) || a >= na || b >= na) {
            std::cerr << "ERROR: invalid genotype data: " << v[i].to_string() << "\n";
            return 1;
        }

        if (i == 9)
            e.ploidy = info;

        if (info != e.ploidy) {
            std::cerr << "ERROR: ploidy doesn't match: " << v[i].to_string() << "\n";
            return 1;
        }

        e.gt.push_back(a < 0 ? 0 : static_cast<allele_t>(a+1));
        if (e.ploidy == 2)
            e.gt.push_back(b < 0 ? 0 : static_cast<allele_t>(b+1));
    }

    return 0;
}

int read_vcf(const std::string &filename, Genotype &gt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    size_t ln = 0;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;

        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        if (line.compare(0, 2, "##") == 0)
            continue;

        if (line.compare(0, 1, "#") == 0) {
            if (parse_vcf_header(line, gt.ind) != 0)
                return 1;
            break;
        }

        std::cerr << "ERROR: VCF header line is required\n";
        return 1;
    }

    VcfEntry e;

    for (std::string line; std::getline(ifs, line);) {
        ++ln;

        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        if (parse_vcf_entry(line, e) != 0)
            return 1;

        if (gt.ploidy <= 0)
            gt.ploidy = e.ploidy;

        if (e.ploidy != gt.ploidy) {
            std::cerr << "ERROR: ploidy doesn't match at line " << ln << "\n";
            return 1;
        }

        if (e.gt.size() != static_cast<size_t>(e.ploidy) * gt.ind.size()) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << "\n";
            return 1;
        }

        gt.loc.push_back(e.id);
        gt.chr.push_back(e.chr);
        gt.pos.push_back(e.pos);
        gt.allele.push_back(e.as);
        gt.dat.push_back(e.gt);
    }

    return 0;
}

int write_vcf(const Genotype & gt, const std::string & filename, bool force_diploid)
{
    std::ofstream ofs(filename);
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    std::vector<std::string> codes(256);
    codes[0] = ".";
    for (size_t i = 1; i < 256; ++i)
        codes[i] = std::to_string(i-1);

    auto n = gt.ind.size();

    ofs << "##fileformat=VCFv4.2\n";
    ofs << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    if ( ! gt.ind.empty() ) {
        ofs << "\tFORMAT";
        for (size_t i = 0; i < n; ++i)
            ofs << "\t" << gt.ind[i];
    }
    ofs << "\n";

    std::string line;
    auto m = gt.loc.size();
    bool haploid = gt.ploidy != 2;

    for (size_t j = 0; j < m; ++j) {
        ofs << gt.chr[j] << "\t" << gt.pos[j] << "\t" << gt.loc[j] << "\t";

        if ( gt.allele[j].empty() )
            ofs << ".\t.";
        else {
            auto na = gt.allele[j].size();
            ofs << gt.allele[j][0] << "\t";
            if (na == 1) {
                ofs << ".";
            }
            else {
                ofs << gt.allele[j][1];
                for (size_t k = 2; k < na; ++k)
                    ofs << "," << gt.allele[j][k];
            }
        }

        ofs << "\t.\t.\t.";

        if ( gt.ind.empty() ) {
            ofs << "\n";
            continue;
        }

        line.clear();
        line.append("\tGT");

        if ( haploid ) {
            for (size_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i];
                line.append("\t").append(codes[a]);
                if ( force_diploid )
                    line.append("/").append(codes[a]);
            }
        }
        else {
            for (size_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i*2], b = gt.dat[j][i*2+1];
                line.append("\t").append(codes[a]).append("/").append(codes[b]);
            }
        }

        ofs << line << "\n";
    }

    return 0;
}
