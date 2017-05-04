#include <limits>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "vcf.h"
#include "split.h"
#include "number.h"

namespace {

int parse_vcf_header(const string &str, vector<string> &ind)
{
    vector<string> vs;
    split([](char c) { return c == '\t'; }, str.begin(), str.end(), vs);

    if (vs.size() != 8 && vs.size() < 10) {
        std::cerr << "ERROR: incorrect number of columns in header line: " << vs.size() << "\n";
        return 1;
    }

    if (vs[0] != "#CHROM" || vs[1] != "POS" || vs[2] != "ID" || vs[3] != "REF" ||
        vs[4] != "ALT" || vs[5] != "QUAL" || vs[6] != "FILTER" || vs[7] != "INFO") {
        std::cerr << "ERROR: incorrect column names in header line: "
                  << vs[0] << "\t" << vs[1] << "\t" << vs[2] << "\t" << vs[3] << "\t"
                  << vs[4] << "\t" << vs[5] << "\t" << vs[6] << "\t" << vs[7] << "\n";
        return 1;
    }

    if (vs.size() > 9) {
        if (vs[8] != "FORMAT") {
            std::cerr << "ERROR: FORMAT is required at 9th field in header line: " << vs[8] << "\n";
            return 1;
        }
        ind.assign(vs.begin() + 9, vs.end());
    }

    return 0;
}

int parse_vcf_gt(const Token &t, int &a, int &b)
{
    auto beg = t.data();
    auto end = std::find(t.data(), t.data() + t.size(), ':');
    size_t len = end - beg;

    a = b = -9;

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
        std::cerr << "ERROR: unsuppored polyploidy genotype: " << string(beg,end) << "\n";
        return -1;
    }

    size_t pos = 0;
    while (pos < len && beg[pos] != '/' && beg[pos] != '|')
        ++pos;
    if (pos == 0 || pos == len - 1)
        return -1;

    bool ok = false;

    if (pos == len) {
        a = number<int>(string(beg,end), &ok);
        return ! ok || a < 0 ? -1 : 1;
    }

    string gt(beg, beg + pos);
    if (gt != ".") {
        a = number<int>(gt, &ok);
        if ( ! ok || a < 0 )
            return -1;
    }

    gt.assign(beg + pos + 1, end);
    if (gt != ".") {
        b = number<int>(gt, &ok);
        if ( ! ok || b < 0 )
            return -1;
    }

    return 2;
}

} // namespace

int read_vcf(const string &filename, Genotype &gt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    string ver;

    size_t ln = 0;
    for (string line; std::getline(ifs,line); ) {
        ++ln;
        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        if (ln == 1 && line.compare(0,13,"##fileformat=") == 0) {
            ver = line.substr(13);
            if (ver != "VCFv4.0" && ver != "VCFv4.1" && ver != "VCFv4.2" && ver != "VCFv4.3") {
                std::cerr << "ERROR: unsupported VCF file version: " << line << "\n";
                return 1;
            }
        }

        if (line.compare(0,7,"##INFO=") == 0)
            continue;

        if (line.compare(0,9,"##FILTER=") == 0)
            continue;

        if (line.compare(0,9,"##FORMAT=") == 0)
            continue;

        if (line.compare(0,9,"##contig=") == 0)
            continue;

        if (line.compare(0,2,"##") == 0)
            continue;

        if (line.compare(0,1,"#") == 0) {
            if ( parse_vcf_header(line, gt.ind) )
                return 1;
            break;
        }

        std::cerr << "ERROR: #CHROME header line is required\n";
        return 1;
    }

    auto delim = [](char c) { return c == '\t'; };

    int ploidy = -1;
    size_t ncols = gt.ind.empty() ? 8 : gt.ind.size() + 9;

    for (string line; std::getline(ifs,line); ) {
        ++ln;
        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        vector<Token> vt;
        split(delim, line.begin(), line.end(), vt);

        if (vt.size() != ncols) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " ("
                      << vt.size() << " != " << ncols << "): " << filename << "\n";
            return 1;
        }

        gt.chr.push_back(string(vt[0]));
        gt.pos.push_back(number<int>(string(vt[1])));

        if (vt[2] == Token(".",1))
            gt.loc.push_back(string(vt[0]) + "_" + string(vt[1]));
        else
            gt.loc.push_back(string(vt[2]));

        if ( gt.ind.empty() )
            continue;

        auto format = string(vt[8]);
        if (format.compare(0,2,"GT") != 0) {
            std::cerr << "ERROR: GT is required and must be the first sub-field of FORMAT: "
                      << format << "\n";
            return 1;
        }

        auto ref = string(vt[3]);
        auto alt = string(vt[4]);
        std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
        std::transform(alt.begin(), alt.end(), alt.begin(), ::toupper);

        vector<string> as;
        as.push_back(ref);
        split([](char c) { return c == ','; }, alt.begin(), alt.end(), as);
        gt.allele.push_back(as);
        int na = as.size();

        vector<allele_t> v;
        for (size_t i = 9; i < ncols; ++i) {
            int a = -9, b = -9;
            int info = parse_vcf_gt(vt[i], a, b);
            if (info < 0 || a >= na || b >= na) {
                std::cerr << "ERROR: invalid genotype data: " << string(vt[i]) << "\n";
                return 1;
            }

            if (ploidy == -1)
                ploidy = info;

            if (info != ploidy) {
                std::cerr << "ERROR: inconsistent ploidy (" << ploidy
                          << ") genotype: " << string(vt[i]) << "\n";
                return 1;
            }

            v.push_back(a < 0 ? 0 : a + 1);
            if (ploidy == 2)
                v.push_back(b < 0 ? 0 : b + 1);
        }

        gt.dat.push_back(v);
    }

    gt.dist.assign(gt.loc.size(), 0.0);
    gt.ploidy = ploidy;

    return 0;
}

int write_vcf(const Genotype &gt, const string &filename, bool diploid)
{
    std::ofstream ofs(filename);
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto m = gt.loc.size(), n = gt.ind.size();
    bool haploid = gt.ploidy == 1;

    ofs << "##fileformat=VCFv4.2\n";

    ofs << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    if ( ! gt.ind.empty() ) {
        ofs << "\tFORMAT";
        for (size_t i = 0; i < n; ++i)
            ofs << "\t" << gt.ind[i];
    }
    ofs << "\n";

    int p = std::numeric_limits<allele_t>::max() + 1;
    vector<string> gs(p);
    gs[0] = ".";
    for (int i = 1; i < p; ++i)
        gs[i] = std::to_string(i-1);

    string line;

    for (size_t j = 0; j < m; ++j) {
        line.clear();

        line.append(gt.chr[j]).append("\t");
        line.append(std::to_string(gt.pos[j])).append("\t");
        line.append(gt.loc[j]).append("\t");

        if ( ! gt.allele[j].empty() ) {
            line.append(gt.allele[j][0]).append("\t");
            auto na = gt.allele[j].size();
            if (na > 1) {
                for (size_t k = 1; k < na; ++k)
                    line.append(gt.allele[j][k]).append(",");
                line.pop_back();
            }
            else
                line.append(".");
        }
        else
            line.append(".\t.");

        line.append("\t.\t.\t.");

        if ( gt.ind.empty() ) {
            ofs << line << "\n";
            continue;
        }

        line.append("\tGT");

        if ( haploid ) {
            for (size_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i];
                line.append("\t").append(gs[a]);
                if ( diploid )
                    line.append("/").append(gs[a]);
            }
        }
        else {
            for (size_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i*2], b = gt.dat[j][i*2+1];
                line.append("\t").append(gs[a]).append("/").append(gs[b]);
            }
        }

        ofs << line << "\n";
    }

    return 0;
}
