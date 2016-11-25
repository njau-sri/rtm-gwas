#include <fstream>
#include <iostream>
#include "plinkio.h"
#include "strsplit.h"
#include "util.h"

// PLINK, http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml

namespace {

int read_plink_map(const string &filename, Genotype &gt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    auto delim = [](char c) { return c == ' ' || c == '\t' || c == '\r'; };

    size_t ln = 0;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<string> vs;
        strsplit(delim, line.begin(), line.end(), vs);
        if ( vs.empty() ) {
            std::cerr << "WARN: skipping empty line: " << ln << "\n";
            continue;
        }

        if (vs.size() != 4) {
            std::cerr << "ERROR: expected 4 columns at line " << ln << ": " << filename << "\n";
            return 1;
        }

        gt.chr.push_back(vs[0]);
        gt.loc.push_back(vs[1]);
        gt.dist.push_back(number<double>(vs[2]));
        gt.pos.push_back(number<int>(vs[3]));
    }

    return 0;
}

int read_plink_ped(const string &filename, Genotype &gt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    auto delim = [](char c) { return c == ' ' || c == '\t' || c == '\r'; };

    size_t ln = 0;
    auto m = gt.loc.size();
    vector< vector<allele_t> > dat;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<Token> vt;
        strsplit(delim, line.begin(), line.end(), vt);
        if ( vt.empty() ) {
            std::cerr << "WARN: skipping empty line: " << ln << "\n";
            continue;
        }

        if (vt.size() != 6 + 2*m) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vt.size() << " != "
                      << 6 + 2*m << "): " << filename << "\n";
            return 1;
        }

        gt.ind.push_back(string(vt[0]) + "_" + string(vt[1]));

        vector<allele_t> v;
        for (auto itr = vt.begin() + 6; itr != vt.end(); ++itr) {
            if (itr->size() != 1) {
                std::cerr << "ERROR: invalid PLINK genotype code: " << string(*itr) << "\n";
                return 1;
            }
            v.push_back((*itr)[0]);
        }

        dat.push_back(v);
    }

    auto n = dat.size();
    const allele_t missing = '0';

    for (size_t j = 0; j < m; ++j) {
        vector<allele_t> v;
        for (size_t i = 0; i < n; ++i) {
            v.push_back(dat[i][j*2]);
            v.push_back(dat[i][j*2+1]);
        }

        auto u = v;
        std::sort(u.begin(), u.end());
        u.erase(std::unique(u.begin(), std::remove(u.begin(), u.end(), missing)), u.end());
        if (u.size() > 2) {
            std::cerr << "ERROR: exceed the maximum number of alleles: " << u.size() << "\n";
            return 1;
        }

        vector<string> allele;
        for (auto a : u)
            allele.emplace_back(1,a);
        gt.allele.push_back(allele);

        for (auto &a : v)
            a = a == missing ? 0 : index(u,a)+1;

        gt.dat.push_back(v);
    }

    gt.ploidy = 2;

    return 0;
}

} // namespace

bool is_compatible_plink(const Genotype &gt)
{
    for (auto &v : gt.allele) {
        for (auto &e : v) {
            if (e.size() != 1)
                return false;
            auto a = e[0];
            if (a != 'A' && a != 'C' && a != 'G' && a != 'T' &&
                a != '1' && a != '2' && a != '3' && a != '4')
                return false;
        }
    }

    return true;
}

int read_plink(const string &prefix, Genotype &gt)
{
    int info = read_plink_map(prefix + ".map", gt);

    if (info == 0)
        info = read_plink_ped(prefix + ".ped", gt);

    return info;
}

int write_plink(const Genotype &gt, const string &prefix)
{
    std::ofstream pedf(prefix + ".ped");
    if ( ! pedf ) {
        std::cerr << "ERROR: can't open file for writing: " << prefix << ".ped\n";
        return 1;
    }

    std::ofstream mapf(prefix + ".map");
    if ( ! mapf ) {
        std::cerr << "ERROR: can't open file for writing: " << prefix << ".map\n";
        return 1;
    }

    auto m = gt.loc.size(), n = gt.ind.size();
    bool haploid = gt.ploidy == 1;

    string line;

    for (size_t i = 0; i < n; ++i) {
        line.clear();
        line.append(gt.ind[i]).append(" 0 0 0 0");
        auto k1 = haploid ? i : i*2;
        auto k2 = haploid ? i : i*2+1;
        for (size_t j = 0; j < m; ++j) {
            auto a = gt.dat[j][k1];
            auto b = gt.dat[j][k2];
            if (a == 0 || b == 0)
                line.append(" 0 0");
            else {
                line.push_back(' ');
                line.append(gt.allele[j][a-1]);
                line.push_back(' ');
                line.append(gt.allele[j][b-1]);
            }
        }
        pedf << i+1 << " " << line << "\n";
    }

    for (size_t j = 0; j < m; ++j)
        mapf << gt.chr[j] << " " << gt.loc[j] << " " << gt.dist[j] << " " << gt.pos[j] << "\n";

    return 0;
}
