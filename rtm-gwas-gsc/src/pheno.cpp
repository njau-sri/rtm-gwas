#include "pheno.h"

#include <cmath>
#include <limits>

#include "cfile.h"
#include "print.h"
#include "stringutil.h"
#include "vectorutil.h"

namespace {

    const std::string kENV = "_ENV_";
    const std::string kBLK = "_BLK_";

} // namespace

int read_pheno(const std::string &filename, Phenotype &pt)
{
    CFileLineReader file(filename);

    if (!file) {
        eprint("ERROR: can't open file for reading: %s\n", filename);
        return 1;
    }

    isize_t ln = 0;
    std::vector<std::string> header;

    for (std::string line; file.read(line); ) {
        ++ln;
        split(line, " \t", header);
        if (!header.empty())
            break;
    }

    if (count(header, kENV) > 1) {
        eprint("ERROR: multiple keywords %s is not allowed\n", kENV);
        return 1;
    }

    if (count(header, kBLK) > 1) {
        eprint("ERROR: multiple keywords %s is not allowed\n", kBLK);
        return 1;
    }

    std::vector<isize_t> jphe;
    isize_t jenv = -1, jblk = -1;
    isize_t nhead = length(header);

    for (isize_t j = 1; j < nhead; ++j) {
        if (header[j] == kENV)
            jenv = j;
        else if (header[j] == kBLK)
            jblk = j;
        else
            jphe.push_back(j);
    }

    for (auto j : jphe)
        pt.phe.push_back(header[j]);

    std::vector< std::vector<double> > dat;

    for (std::string line; file.read(line); ) {
        ++ln;

        auto vs = split(line, " \t");
        if (vs.empty())
            continue;

        if (length(vs) != nhead) {
            eprint("ERROR: column count doesn't match at line %td: %td != %td\n", ln, length(vs), nhead);
            return 1;
        }

        pt.ind.push_back(vs[0]);

        if (jenv != -1)
            pt.env.push_back(vs[jenv]);

        if (jblk != -1)
            pt.blk.push_back(vs[jblk]);

        std::vector<double> v;

        for (auto j : jphe) {
            if (vs[j] == "?" || vs[j] == "NA" || vs[j] == ".")
                v.push_back(std::numeric_limits<double>::quiet_NaN());
            else
                v.push_back(std::stod(vs[j]));
        }

        dat.push_back(v);
    }

    isize_t m = length(pt.phe);
    isize_t n = length(pt.ind);

    for (isize_t j = 0; j < m; ++j) {
        std::vector<double> v(n);
        for (isize_t i = 0; i < n; ++i)
            v[i] = dat[i][j];
        pt.dat.push_back(v);
    }

    if (!pt.env.empty() && length(unique(pt.env)) < 2) {
        pt.env.clear();
        eprint("WARNING: ignoring invalid %s factor\n", kENV);
    }

    if (!pt.blk.empty() && length(unique(pt.blk)) < 2) {
        pt.blk.clear();
        eprint("WARNING: ignoring invalid %s factor\n", kBLK);
    }

    return 0;
}

int write_pheno(const Phenotype &pt, const std::string &filename)
{
    CFile file(filename, "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s\n", filename);
        return 1;
    }

    fprint(file, "Indiv");
    if (!pt.env.empty())
        fprint(file, "\t%s", kENV);
    if (!pt.blk.empty())
        fprint(file, "\t%s", kBLK);
    for (auto &e : pt.phe)
        fprint(file, "\t%s", e);
    fprint(file, "\n");

    isize_t m = length(pt.phe);
    isize_t n = length(pt.ind);

    for (isize_t i = 0; i < n; ++i) {
        fprint(file, "%s", pt.ind[i]);
        if (!pt.env.empty())
            fprint(file, "\t%s", pt.env[i]);
        if (!pt.blk.empty())
            fprint(file, "\t%s", pt.blk[i]);
        for (isize_t j = 0; j < m; ++j)
            fprint(file, "\t%g", pt.dat[j][i]);
        fprint(file, "\n");
    }

    return 0;
}

int read_covar(const std::string &filename, Covariate &ct)
{
    CFileLineReader file(filename);

    if (!file) {
        eprint("ERROR: can't open file for reading: %s\n", filename);
        return 1;
    }

    isize_t ln = 0;
    std::vector<std::string> header;

    for (std::string line; file.read(line); ) {
        ++ln;
        split(line, " \t", header);
        if (!header.empty())
            break;
    }

    isize_t nhead = length(header);
    for (isize_t j = 1; j < nhead; ++j)
        ct.phe.push_back(header[j]);

    std::vector< std::vector<double> > dat;

    for (std::string line; file.read(line); ) {
        ++ln;
        auto vs = split(line, " \t");
        if (vs.empty())
            continue;

        if (length(vs) != nhead) {
            eprint("ERROR: column count doesn't match at line %td: %td != %td\n", ln, length(vs), nhead);
            return 1;
        }

        ct.ind.push_back(vs[0]);
        vs.erase(vs.begin());

        std::vector<double> v;
        for (auto &e : vs) {
            v.push_back(std::stod(e));
            if ( ! std::isfinite(v.back()) ) {
                eprint("ERROR: only finite value is allowed: %s\n", e);
                return 1;
            }
        }
        dat.push_back(v);
    }

    isize_t m = length(ct.phe);
    isize_t n = length(ct.ind);

    for (isize_t j = 0; j < m; ++j) {
        std::vector<double> v(n);
        for (isize_t i = 0; i < n; ++i)
            v[i] = dat[i][j];
        ct.dat.push_back(v);
    }

    return 0;
}

int write_covar(const Covariate &ct, const std::string &filename)
{
    CFile file(filename, "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s\n", filename);
        return 1;
    }

    fprint(file, "Indiv");
    for (auto &e : ct.phe)
        fprint(file, "\t%s", e);
    fprint(file, "\n");

    isize_t m = length(ct.phe);
    isize_t n = length(ct.ind);

    for (isize_t i = 0; i < n; ++i) {
        fprint(file, "%s", ct.ind[i]);
        for (isize_t j = 0; j < m; ++j)
            fprint(file, "\t%g", ct.dat[j][i]);
        fprint(file, "\n");
    }

    return 0;
}

int read_square(const std::string &filename, SquareData &sd)
{
    CFileLineReader file(filename);

    if (!file) {
        eprint("ERROR: can't open file for reading: %s\n", filename);
        return 1;
    }

    isize_t ln = 0, cc = 0;

    for (std::string line; file.read(line); ) {
        ++ln;

        auto vs = split(line, " \t");
        if (vs.empty())
            continue;

        if (cc < 1)
            cc = length(vs);

        if (length(vs) != cc) {
            eprint("ERROR: column count doesn't match at line %td: %td != %td\n", ln, length(vs), cc);
            return 1;
        }

        sd.ind.push_back(vs[0]);
        vs.erase(vs.begin());

        std::vector<double> v;

        for (auto &e : vs) {
            v.push_back(std::stod(e));
            if ( ! std::isfinite(v.back()) ) {
                eprint("ERROR: only finite value is allowed: %s\n", e);
                return 1;
            }
        }

        sd.dat.push_back(v);
    }

    if (cc > 0 && cc != length(sd.ind) + 1) {
        eprint("ERROR: data is not square in file: %s\n", filename);
        return 1;
    }

    return 0;
}

int write_square(const SquareData &sd, const std::string &filename)
{
    CFile file(filename, "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s\n", filename);
        return 1;
    }

    isize_t n = length(sd.ind);

    for (isize_t i = 0; i < n; ++i) {
        fprint(file, "%s", sd.ind[i]);
        for (isize_t j = 0; j < n; ++j)
            fprint(file, "\t%g", sd.dat[j][i]);
        fprint(file, "\n");
    }

    return 0;
}

int read_qtl_effect(const std::string &filename, QtlEffect &qe)
{
    CFileLineReader file(filename);

    if (!file) {
        eprint("ERROR: can't open file for reading: %s\n", filename);
        return 1;
    }

    int ln = 0;
    std::string curr;

    for (std::string line; file.read(line); ) {
        ++ln;

        auto vs = split(line, " \t");

        if (vs.empty())
            continue;

        if (length(vs) == 1 && vs[0].find('>') == 0) {
            curr = vs[0].substr(1);
            continue;
        }

        if (length(vs) != 3)
            continue;

        if (vs[0].find("_ENV_") != std::string::npos)
            continue;

        if (vs[0].find("_BLK_") != std::string::npos)
            continue;

        if (!curr.empty()) {
            std::vector<std::string> as;
            split(vs[1], ":/|", as);
            if (length(as) == 1 || (length(as) == 2 && as[0] == as[1])) {
                qe.phe.push_back(curr);
                qe.qtl.push_back(vs[0]);
                qe.allele.push_back(as[0]);
                qe.effect.push_back(std::stod(vs[2]));
            }
        }
    }

    return 0;
}

int read_map(const std::string &filename, GeneticMap &gm)
{
    CFileLineReader file(filename);

    if (!file) {
        eprint("ERROR: can't open file for reading: %s\n", filename);
        return 1;
    }

    int ln = 0;

    for (std::string line; file.read(line); ) {
        ++ln;

        auto vs = split(line, " \t");

        if (vs.empty())
            continue;

        if (length(vs) < 3) {
            eprint("ERROR: expected at least three columns at line %d\n", ln);
            return 1;
        }

        gm.loc.push_back(vs[0]);
        gm.chr.push_back(vs[1]);

        double a = std::stod(vs[2]);
        gm.pos.push_back(a);

        if (!std::isfinite(a) || a < 0.0) {
            eprint("ERROR: invalid map position: %s\n", line);
            return 1;
        }
    }

    return 0;
}
