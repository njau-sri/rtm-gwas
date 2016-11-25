#include <fstream>
#include <iostream>
#include "rtmio.h"
#include "strsplit.h"
#include "util.h"

namespace {

// http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html

char encode_iupac(char a, char b)
{
    if (a == b && (a == 'A' || a == 'C' || a == 'G' || a == 'T'))
        return a;

    if (a > b)
        std::swap(a, b);

    if (a == 'A') {
        if (b == 'C') return 'M';
        if (b == 'G') return 'R';
        if (b == 'T') return 'W';
    }

    if (a == 'C') {
        if (b == 'G') return 'S';
        if (b == 'T') return 'Y';
    }

    if (a == 'G' && b == 'T')
        return 'K';

    return 'N';
}

std::pair<char,char> decode_iupac(char c)
{
    char a = 0, b = 0;

    if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
        a = b = c;
    else if (c == 'W')
        a = 'A', b = 'T';
    else if (c == 'S')
        a = 'C', b = 'G';
    else if (c == 'M')
        a = 'A', b = 'C';
    else if (c == 'K')
        a = 'G', b = 'T';
    else if (c == 'R')
        a = 'A', b = 'G';
    else if (c == 'Y')
        a = 'C', b = 'T';

    return {a, b};
}

bool is_compatible_iupac(const Genotype &gt)
{
    for (auto &v : gt.allele) {
        for (auto &e : v) {
            if (e.size() != 1)
                return false;
            auto a = e[0];
            if (a != 'A' && a != 'C' && a != 'G' && a != 'T')
                return false;
        }
    }

    return true;
}

bool is_iupac(const vector< vector<allele_t> > &dat)
{
    bool ret = false;

    for (auto &v : dat) {
        for (auto &a : v) {
            switch (a) {
            case 'A': case 'C': case 'G': case 'T': case 0:
                break;
            case 'W': case 'S': case 'M': case 'K': case 'R': case 'Y':
                ret = true;
            default:
                return false;
            }
        }
    }

    return ret;
}

int read_genotype_char(const string &filename, Genotype &gt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    size_t ln = 0;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<string> vs;
        strsplit([](char c) { return c == ' ' || c == '\t' || c == '\r'; }, line.begin(), line.end(), vs);
        if ( vs.empty() ) {
            std::cerr << "WARN: skipping empty line: " << ln << "\n";
            continue;
        }

        if (vs.size() < 4) {
            std::cerr << "ERROR: expected at least 4 columns at line " << ln << ": " << filename << "\n";
            return 1;
        }

        gt.ind.assign(vs.begin() + 4, vs.end());

        break;
    }

    auto delim = [](char c) { return c == ' ' || c == '\t' ||  c == ':' || c == '\r'; };

    size_t ploidy = 0;
    auto n = gt.ind.size();

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<Token> vt;
        strsplit(delim, line.begin(), line.end(), vt);
        if ( vt.empty() ) {
            std::cerr << "WARN: skipping empty line: " << ln << "\n";
            continue;
        }

        if (ploidy == 0) {
            ploidy = vt.size() > 4 + n ? (vt.size() - 4) / n : 1;
            if (ploidy > 2) {
                std::cerr << "ERROR: unsupported polyploidy genotype data: " << ploidy << "\n";
                return 1;
            }
        }

        if (vt.size() != 4 + ploidy * n) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vt.size() << " != "
                      << 4 + ploidy * n << "): " << filename << "\n";
            return 1;
        }

        gt.loc.push_back(string(vt[0]));
        gt.chr.push_back(string(vt[1]));
        gt.pos.push_back(number<int>(string(vt[2])));
        gt.dist.push_back(number<double>(string(vt[3])));

        vector<allele_t> v;
        for (auto itr = vt.begin() + 4; itr != vt.end(); ++itr) {
            if (itr->size() == 1)
                v.push_back((*itr)[0]);
            else
                return -1;
        }

        for (auto &e : v) {
            if (e == 'N' || e == '-' || e == '.' || e == '?')
                e = 0;
        }

        gt.dat.push_back(v);
    }

    bool iupac = ploidy == 1 && is_iupac(gt.dat);

    for (auto &v : gt.dat) {
        if ( iupac ) {
            vector<allele_t> w;
            w.reserve( v.size() * 2 );
            for (auto a : v) {
                auto p = decode_iupac(a);
                w.push_back(p.first);
                w.push_back(p.second);
            }
            v.swap(w);
        }

        auto u = v;
        std::sort(u.begin(), u.end());
        u.erase(std::unique(u.begin(), std::remove(u.begin(), u.end(), allele_t(0))), u.end());

        vector<string> allele;
        for (auto a : u)
            allele.emplace_back(1,a);
        gt.allele.push_back(allele);

        for (auto &a : v)
            a = a == 0 ? 0 : index(u,a)+1;
    }

    gt.ploidy = ploidy;

    return 0;
}

int read_genotype_string(const string &filename, Genotype &gt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    size_t ln = 0;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<string> vs;
        strsplit([](char c) { return c == ' ' || c == '\t' || c == '\r'; }, line.begin(), line.end(), vs);
        if ( vs.empty() ) {
            std::cerr << "WARN: skipping empty line: " << ln << "\n";
            continue;
        }

        if (vs.size() < 4) {
            std::cerr << "ERROR: expected at least 4 columns at line " << ln << ": " << filename << "\n";
            return 1;
        }

        gt.ind.assign(vs.begin() + 4, vs.end());

        break;
    }

    auto delim = [](char c) { return c == ' ' || c == '\t' || c == ':' || c == '\r'; };

    size_t ploidy = 0;
    auto n = gt.ind.size();
    const Token missing("?", 1);

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<Token> vt;
        strsplit(delim, line.begin(), line.end(), vt);
        if ( vt.empty() ) {
            std::cerr << "WARN: skipping empty line: " << ln << "\n";
            continue;
        }

        if (ploidy == 0) {
            ploidy = vt.size() > 4 + n ? (vt.size() - 4) / n : 1;
            if (ploidy > 2) {
                std::cerr << "ERROR: unsupported polyploidy genotype data: " << ploidy << "\n";
                return 1;
            }
        }

        if (vt.size() != 4 + ploidy * n) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vt.size() << " != "
                      << 4 + ploidy * n << "): " << filename << "\n";
            return 1;
        }

        gt.loc.push_back(string(vt[0]));
        gt.chr.push_back(string(vt[1]));
        gt.pos.push_back(number<int>(string(vt[2])));
        gt.dist.push_back(number<double>(string(vt[3])));

        vector<Token> u(vt.begin() + 4, vt.end());
        std::sort(u.begin(), u.end());
        u.erase(std::unique(u.begin(), std::remove(u.begin(), u.end(), missing)), u.end());

        if (u.size() > std::numeric_limits<allele_t>::max()) {
            std::cerr << "ERROR: exceed the maximum number of alleles: " << u.size() << "\n";
            return 1;
        }

        vector<string> allele;
        for (auto &e : u)
            allele.push_back(string(e));
        gt.allele.push_back(allele);

        vector<allele_t> v;
        for (auto itr = vt.begin() + 4; itr != vt.end(); ++itr) {
            if (*itr == missing)
                v.push_back(0);
            else
                v.push_back(index(u,*itr)+1);
        }

        gt.dat.push_back(v);
    }

    gt.ploidy = ploidy;

    return 0;
}

} // namespace

int read_file(const string &filename, string &contents)
{
    std::ifstream ifs(filename, std::ios::binary);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    ifs.seekg(0, std::ios::end);
    auto n = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    contents.resize(n);
    if ( ! contents.empty() )
        ifs.read(&contents[0], n);

    return 0;
}

int write_file( const string &contents, const string &filename)
{
    std::ofstream ofs(filename);

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    ofs << contents;

    return 0;
}

int read_genotype(const string &filename, Genotype &gt)
{
    int info = read_genotype_char(filename, gt);

    if (info < 0) {
        gt.loc.clear(); gt.chr.clear(); gt.pos.clear();
        gt.dist.clear(); gt.dat.clear(); gt.allele.clear();
        info = read_genotype_string(filename, gt);
    }

    return info;
}

int write_genotype(const Genotype &gt, const string &filename)
{
    std::ofstream ofs(filename);
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto m = gt.loc.size(), n = gt.ind.size();
    bool haploid = gt.ploidy == 1;
    bool iupac = is_compatible_iupac(gt);
    const string missing = "?";

    string line;

    ofs << "Locus\tChromosome\tPosition\tDistance";
    for (size_t i = 0; i < n; ++i)
        ofs << "\t" << gt.ind[i];
    ofs << "\n";

    for (size_t j = 0; j < m; ++j) {
        ofs << gt.loc[j] << "\t" << gt.chr[j] << "\t" << gt.pos[j] << "\t" << gt.dist[j];

        line.clear();

        for (size_t i = 0; i < n; ++i) {
            line.push_back('\t');
            if ( haploid ) {
                auto a = gt.dat[j][i];
                line.append(a == 0 ? missing : gt.allele[j][a-1]);
            }
            else {
                auto a = gt.dat[j][i*2];
                auto b = gt.dat[j][i*2+1];
                if ( iupac ) {
                    if (a == 0 || b == 0)
                        line.push_back('N');
                    else
                        line.push_back( encode_iupac(gt.allele[j][a-1][0], gt.allele[j][b-1][0]) );
                }
                else {
                    line.append(a == 0 ? missing : gt.allele[j][a-1]);
                    line.push_back(':');
                    line.append(b == 0 ? missing : gt.allele[j][b-1]);
                }
            }
        }

        ofs << line << "\n";
    }

    return 0;
}

int read_phenotype(const string &filename, Phenotype &pt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    auto delim = [](char c) { return c == ' ' || c == '\t' || c == '\r'; };

    size_t ln = 0;
    vector<string> colnames;

    for (string line; std::getline(ifs,line); ) {
        ++ln;
        strsplit(delim, line.begin(), line.end(), colnames);
        if ( ! colnames.empty() )
            break;
        std::cerr << "WARN: skipping empty line: " << ln << "\n";
    }

    vector<size_t> jphe;
    size_t jenv = 0, jblk = 0;

    for (size_t j = 1; j < colnames.size(); ++j) {
        if (colnames[j] == "_ENV_")
            jenv = j;
        else if (colnames[j] == "_BLK_")
            jblk = j;
        else
            jphe.push_back(j);
    }

    for (auto j : jphe)
        pt.phe.push_back(colnames[j]);

    vector< vector<double> > dat;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<string> vs;
        strsplit(delim, line.begin(), line.end(), vs);
        if ( vs.empty() ) {
            std::cerr << "WARN: skipping empty line: " << ln << "\n";
            continue;
        }

        if (vs.size() != colnames.size()) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vs.size() << "!="
                      << colnames.size() << "): " << filename << "\n";
            return 1;
        }

        pt.ind.push_back(vs[0]);

        if (jenv > 0)
            pt.env.push_back(vs[jenv]);

        if (jblk > 0)
            pt.blk.push_back(vs[jblk]);

        vector<double> v;

        for (auto j : jphe) {
            if (vs[j] == "?" || vs[j] == "NA" || vs[j] == ".")
                v.push_back(std::numeric_limits<double>::quiet_NaN());
            else
                v.push_back(number<double>(vs[j]));
        }

        dat.push_back(v);
    }

    auto m = pt.phe.size(), n = pt.ind.size();

    for (size_t j = 0; j < m; ++j) {
        vector<double> v(n);
        for (size_t i = 0; i < n; ++i)
            v[i] = dat[i][j];
        pt.dat.push_back(v);
    }

    return 0;
}

int write_phenotype(const Phenotype &pt, const string &filename)
{
    std::ofstream ofs(filename);

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto m = pt.phe.size();
    auto n = pt.ind.size();

    ofs << "Indiv";
    for (size_t j = 0; j < m; ++j)
        ofs << "\t" << pt.phe[j];
    ofs << "\n";

    for (size_t i = 0; i < n; ++i) {
        ofs << pt.ind[i];
        for (size_t j = 0; j < m; ++j)
            ofs << "\t" << pt.dat[j][i];
        ofs << "\n";
    }

    return 0;
}

int read_square(const string &filename, SquareData &sd)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    auto delim = [](char c) { return c == ' ' || c == '\t' || c == '\r'; };

    size_t ln = 0, cc = 0;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<string> vs;
        strsplit(delim, line.begin(), line.end(), vs);
        if ( vs.empty() ) {
            std::cerr << "WARN: skipping empty line: " << ln << "\n";
            continue;
        }

        if (cc == 0)
            cc = vs.size();

        if (vs.size() != cc) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vs.size() << "!="
                      << cc << "): " << filename << "\n";
            return 1;
        }

        if (sd.ind.size() >= cc) {
            std::cerr << "ERROR: data must be square: " << filename << "\n";
            return 1;
        }

        sd.ind.push_back(vs[0]);

        vector<double> v;
        for (size_t i = 1; i < cc; ++i)
            v.push_back(number<double>(vs[i]));

        sd.dat.push_back(v);
    }

    return 0;
}

int write_square(const SquareData &sd, const string &filename)
{
    std::ofstream ofs(filename);

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto n = sd.ind.size();

    for (size_t i = 0; i < n; ++i) {
        ofs << sd.ind[i];
        for (size_t j = 0; j < n; ++j)
            ofs << "\t" << sd.dat[i][j];
        ofs << "\n";
    }

    return 0;
}
