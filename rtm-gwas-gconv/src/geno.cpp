#include "geno.h"

#include <limits>
#include <algorithm>

#include "cfile.h"
#include "print.h"
#include "stringutil.h"
#include "vectorutil.h"

namespace {

    // http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html
    //
    //   W --- A/T
    //   S --- C/G
    //   M --- A/C
    //   K --- G/T
    //   R --- A/G
    //   Y --- C/T

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

    std::pair<char, char> decode_iupac(char c)
    {
        char a = 0, b = 0;

        if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
            a = b = c;
        else if (c == 'W') {
            a = 'A';
            b = 'T';
        }
        else if (c == 'S') {
            a = 'C';
            b = 'G';
        }
        else if (c == 'M') {
            a = 'A';
            b = 'C';
        }
        else if (c == 'K') {
            a = 'G';
            b = 'T';
        }
        else if (c == 'R') {
            a = 'A';
            b = 'G';
        }
        else if (c == 'Y') {
            a = 'C';
            b = 'T';
        }

        return { a, b };
    }

    bool is_iupac(const std::vector< std::vector<allele_t> > &dat)
    {
        bool ret = false;

        for (auto &v : dat) {
            for (auto &a : v) {
                switch (a) {
                case 'A': case 'C': case 'G': case 'T': case 0:
                    break;
                case 'W': case 'S': case 'M': case 'K': case 'R': case 'Y':
                    ret = true;
                    break;
                default:
                    return false;
                }
            }
        }

        return ret;
    }

    int check_compat_iupac(const Genotype &gt)
    {
        if (gt.source == Genotype::Source::PED)
            return 0;

        for (auto &v : gt.allele) {
            for (auto &e : v) {
                if (length(e) != 1)
                    return 1;
                char a = e[0];
                if (a != 'A' && a != 'C' && a != 'G' && a != 'T')
                    return 2;
            }
        }

        return 0;
    }

    int check_homozygous(const Genotype &gt)
    {
        if (gt.ploidy != 2)
            return 0;

        isize_t n = length(gt.ind);

        for (auto &v : gt.dat) {
            for (isize_t i = 0; i < n; ++i) {
                if (v[i*2] != v[i*2+1])
                    return 1;
            }
        }

        return 0;
    }

    int read_genotype_char(const std::string &filename, Genotype &gt)
    {
        CFileLineReader file(filename);

        if (!file) {
            eprint("ERROR: can't open file for reading: %s\n", filename);
            return 1;
        }

        for (std::string line; file.read(line); ) {
            auto vs = split(line, " \t");
            if (vs.empty())
                continue;

            if (length(vs) < 3) {
                eprint("ERROR: expected at least 3 columns at the first line\n");
                return 1;
            }

            gt.ind.assign(vs.begin() + 3, vs.end());
            break;
        }

        isize_t ploidy = 0;
        isize_t n = length(gt.ind);

        for (std::string line; file.read(line); ) {
            std::vector<Token> vt;
            split(line, " \t/:", vt);
            if (vt.empty())
                continue;
            isize_t nvt = length(vt);
            if (ploidy == 0) {
                ploidy = nvt > (3 + n) ? (nvt - 3) / n : 1;
                if (ploidy > 2) {
                    eprint("ERROR: polyploidy (%td) genotype is not supported: %s\n", ploidy, vt[0].str());
                    return 1;
                }
            }

            if (nvt != 3 + ploidy * n) {
                eprint("ERROR: column count doesn't match (%td != %td): %s\n", nvt, 3 + ploidy * n, vt[0].str());
                return 1;
            }

            gt.loc.push_back(vt[0].str());
            gt.chr.push_back(vt[1].str());
            gt.pos.push_back(std::stoi(vt[2].str()));

            std::vector<allele_t> v;
            for (auto itr = vt.begin() + 3; itr != vt.end(); ++itr) {
                if (length(*itr) == 1)
                    v.push_back((allele_t) ((*itr)[0]));
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
            if (iupac) {
                std::vector<allele_t> w;
                w.reserve(length(v) * 2);
                for (auto a : v) {
                    auto p = decode_iupac(static_cast<char>(a));
                    w.push_back((allele_t) p.first);
                    w.push_back((allele_t) p.second);
                }
                v.swap(w);
            }

            auto u = v;
            std::sort(u.begin(), u.end());
            u.erase(std::unique(u.begin(), std::remove(u.begin(), u.end(), (allele_t) 0)), u.end());

            std::vector<std::string> allele;
            for (auto a : u)
                allele.emplace_back(1, a);
            gt.allele.push_back(allele);

            for (auto &a : v)
                a = a == 0 ? 0 : (allele_t) (index(u, a) + 1);
        }

        gt.ploidy = iupac ? 2 : (int) ploidy;

        return 0;
    }

    int read_genotype_string(const std::string &filename, Genotype &gt)
    {
        CFileLineReader file(filename);

        if (!file) {
            eprint("ERROR: can't open file for reading: %s\n", filename);
            return 1;
        }

        for (std::string line; file.read(line); ) {
            std::vector<std::string> vs;
            split(line, " \t", vs);
            if (vs.empty())
                continue;

            if (length(vs) < 3) {
                eprint("ERROR: expected at least 3 columns at the first line\n");
                return 1;
            }

            gt.ind.assign(vs.begin() + 3, vs.end());
            break;
        }

        isize_t ploidy = 0;
        isize_t n = length(gt.ind);
        static const Token mis1("?", 1);
        static const Token mis2(".", 1);
        static const Token mis3("-", 1);
        static const Token mis4("N", 1);

        for (std::string line; file.read(line); ) {
            std::vector<Token> vt;
            split(line, " \t/:", vt);
            if (vt.empty())
                continue;
            isize_t nvt = length(vt);
            if (ploidy == 0) {
                ploidy = nvt > (3 + n) ? (nvt - 3) / n : 1;
                if (ploidy > 2) {
                    eprint("ERROR: polyploidy (%td) genotype is not supported: %s\n", ploidy, vt[0].str());
                    return 1;
                }
            }

            if (length(vt) != 3 + ploidy * n) {
                eprint("ERROR: column count doesn't match at line (%td != %td): %s\n", nvt, 3 + ploidy * n, vt[0].str());
                return 1;
            }

            gt.loc.push_back(vt[0].str());
            gt.chr.push_back(vt[1].str());
            gt.pos.push_back(std::stoi(vt[2].str()));

            std::vector<Token> u(vt.begin() + 3, vt.end());
            std::sort(u.begin(), u.end());
            auto endmis = std::remove(u.begin(), u.end(), mis1);
            endmis = std::remove(u.begin(), endmis, mis2);
            endmis = std::remove(u.begin(), endmis, mis3);
            endmis = std::remove(u.begin(), endmis, mis4);
            u.erase(std::unique(u.begin(), endmis), u.end());

            isize_t na = length(u);
            if (na > (isize_t) std::numeric_limits<allele_t>::max()) {
                eprint("ERROR: exceed the maximum number of alleles: %td\n", na);
                return 1;
            }

            std::vector<std::string> allele;
            for (auto &e : u)
                allele.push_back(e.str());
            gt.allele.push_back(allele);

            std::vector<allele_t> v;
            for (auto itr = vt.begin() + 3; itr != vt.end(); ++itr) {
                if (*itr == mis1 || *itr == mis2 || *itr == mis3 || *itr == mis4)
                    v.push_back(0);
                else
                    v.push_back((allele_t) (index(u, *itr) + 1));
            }

            gt.dat.push_back(v);
        }

        gt.ploidy = (int) ploidy;

        return 0;
    }

} // namespace

int read_geno(const std::string &filename, Genotype &gt)
{
    int info = read_genotype_char(filename, gt);

    if (info < 0) {
        gt.loc.clear();
        gt.chr.clear();
        gt.pos.clear();
        gt.dat.clear();
        gt.allele.clear();
        info = read_genotype_string(filename, gt);
    }

    gt.source = Genotype::Source::GENO;

    return info;
}

int write_geno(const Genotype &gt, const std::string &filename)
{
    CFile file(filename, "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s\n", filename);
        return 1;
    }

    isize_t m = length(gt.loc);
    isize_t n = length(gt.ind);
    bool haploid = gt.ploidy != 2;
    bool iupac = check_compat_iupac(gt) == 0;
    bool homo = check_homozygous(gt) == 0;
    static const std::string missing = iupac ? "N" : "?";

    std::string line;

    fprint(file, "Locus\tChromosome\tPosition");
    for (isize_t i = 0; i < n; ++i)
        fprint(file, "\t%s", gt.ind[i]);
    fprint(file, "\n");

    for (isize_t j = 0; j < m; ++j) {
        line.clear();

        for (isize_t i = 0; i < n; ++i) {
            line.push_back('\t');
            if (haploid) {
                allele_t a = gt.dat[j][i];
                line.append(a == 0 ? missing : gt.allele[j][a-1]);
            }
            else {
                allele_t a = gt.dat[j][i*2];
                allele_t b = gt.dat[j][i*2+1];
                if (iupac) {
                    if (a == 0 || b == 0)
                        line.append(missing);
                    else
                        line.push_back(encode_iupac(gt.allele[j][a-1][0], gt.allele[j][b-1][0]));
                }
                else {
                    line.append(a == 0 ? missing : gt.allele[j][a-1]);
                    if (!homo) {
                        line.push_back('/');
                        line.append(b == 0 ? missing : gt.allele[j][b-1]);
                    }
                }
            }
        }

        fprint(file, "%s\t%s\t%d", gt.loc[j], gt.chr[j], gt.pos[j]);
        fprint(file, "%s\n", line);
    }

    return 0;
}
