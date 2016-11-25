#include <fstream>
#include <numeric>
#include <iostream>
#include "hapmapio.h"
#include "strsplit.h"
#include "util.h"

bool is_compatible_hapmap(const Genotype &gt)
{
    for (auto &v : gt.allele) {
        for (auto &e : v) {
            if (e.size() != 1)
                return false;
            auto a = e[0];
            if (a != 'A' && a != 'C' && a != 'G' && a != 'T' && a != 'I' && a != 'D')
                return false;
        }
    }

    return true;
}

int read_hapmap(const string &filename, Genotype &gt)
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

        if (vs.size() < 11) {
            std::cerr << "ERROR: expected at least 11 columns at line " << ln << ": " << filename << "\n";
            return 1;
        }

        gt.ind.assign(vs.begin() + 11, vs.end());

        break;
    }

    auto n = gt.ind.size();
    const allele_t missing = 'N';

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<Token> vt;
        strsplit(delim, line.begin(), line.end(), vt);
        if ( vt.empty() ) {
            std::cerr << "WARN: skipping empty line: " << ln << "\n";
            continue;
        }

        if (vt.size() != 11 + n) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vt.size() << " != "
                      << 11 + n << "): " << filename << "\n";
            return 1;
        }

        gt.loc.push_back(string(vt[0]));
        gt.chr.push_back(string(vt[2]));
        gt.pos.push_back(number<int>(string(vt[3])));

        vector<allele_t> v;
        for (auto itr = vt.begin() + 11; itr != vt.end(); ++itr) {
            auto newsize = v.size() + 2;
            if (itr->size() == 2) {
                auto a = (*itr)[0], b = (*itr)[1];
                if (a == 'A' || a == 'C' || a == 'G' || a == 'T' || a == 'I' || a == 'D' || a == 'N')
                    v.push_back(a);
                if (b == 'A' || b == 'C' || b == 'G' || b == 'T' || b == 'I' || b == 'D' || b == 'N')
                    v.push_back(b);
            }
            if (v.size() != newsize) {
                std::cerr << "ERROR: invalid HapMap genotype code: " << string(*itr) << "\n";
                return 1;
            }
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

    gt.dist.assign(gt.loc.size(), 0.0);
    gt.ploidy = 2;

    return 0;
}

int write_hapmap(const Genotype &gt, const string &filename)
{
    std::ofstream ofs(filename);
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto m = gt.loc.size(), n = gt.ind.size();
    bool haploid = gt.ploidy == 1;

    ofs << "rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode";
    for (size_t i = 0; i < n; ++i)
        ofs << " " << gt.ind[i];
    ofs << "\n";

    string line;

    for (size_t j = 0; j < m; ++j) {
        line.clear();
        line.append("NA NA NA NA NA NA NA");
        for (size_t i = 0; i < n; ++i) {
            auto a = haploid ? gt.dat[j][i] : gt.dat[j][i*2];
            auto b = haploid ? gt.dat[j][i] : gt.dat[j][i*2+1];
            if (a == 0 || b == 0)
                line.append(" NN");
            else {
                line.push_back(' ');
                line.append(gt.allele[j][a-1]).append(gt.allele[j][b-1]);
            }
        }

        string allele;
        for (auto e : gt.allele[j]) {
            allele.push_back(e[0]);
            allele.push_back('/');
        }
        if (allele.size() == 2)
            allele.push_back(allele[0]);
        else if ( allele.empty() )
            allele.append("N/N");

        ofs << gt.loc[j] << " " << allele.substr(0,3) << " " << gt.chr[j] << " " << gt.pos[j] << " " << line << "\n";
    }

    return 0;
}
