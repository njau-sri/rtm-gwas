#include "ped.h"

#include <cstring>
#include <algorithm>

#include "cfile.h"
#include "print.h"
#include "stringutil.h"
#include "vectorutil.h"

namespace {

    // fid_iid -> fid iid, otherwise, fid = 1,2,... iid = ind
    void parse_fid_iid(const std::vector<std::string> &ind,
                       std::vector<std::string> &fid,
                       std::vector<std::string> &iid)
    {
        bool decode = true;
        for (auto &e : ind) {
            if (e.front() == '_' || e.back() == '_' || e.find('_') == std::string::npos) {
                decode = false;
                break;
            }
        }

        if (decode) {
            fid.clear();
            iid.clear();
            for (auto &e : ind) {
                auto v = split(e, "_");
                if (length(v) != 2) {
                    decode = false;
                    break;
                }
                fid.push_back(v[0]);
                iid.push_back(v[1]);
            }
        }

        if (!decode) {
            fid.clear();
            iid = ind;
            isize_t n = length(ind);
            for (isize_t i = 0; i < n; ++i)
                fid.push_back(std::to_string(i+1));
        }
    }

    int parse_ped_gt(const char *s, isize_t n, char &a)
    {
        if (n != 1) {
            eprint("ERROR: PED genotype must be represented by 1 character: %s\n", std::string(s, n));
            return 1;
        }

        a = s[0];

        if (a == '1') a = 'A';
        else if (a == '2') a = 'C';
        else if (a == '3') a = 'G';
        else if (a == '4') a = 'T';
        else if (a == '0') a = 'N';
        else if (a != 'A' && a != 'C' && a != 'G' && a != 'T') {
            eprint("ERROR: invalid PED genotype code: %s\n", std::string(s, n));
            return 2;
        }

        return 0;
    }

    int check_compat_ped(const Genotype &gt)
    {
        if (gt.source == Genotype::Source::SNPLDB)
            return 1;

        static const char gs[] = { 'A','C','G','T','1','2','3','4' };

        for (auto &v : gt.allele) {
            if (v.size() > 2)
                return 1;
            for (auto &e : v) {
                if (e.size() != 1)
                    return 2;
                if (std::memchr(gs, e[0], sizeof gs) == NULL)
                    return 3;
            }
        }

        return 0;
    }

    int read_ped_map(const std::string &filename, Genotype &gt)
    {
        CFileLineReader file(filename);

        if (!file) {
            eprint("ERROR: can't open file for reading: %s\n", filename);
            return 1;
        }

        MapEntry entry;

        for (std::string line; file.read(line); ) {
            if (parse_map_entry(line, entry) != 0)
                return 1;
            gt.loc.push_back(entry.id);
            gt.chr.push_back(entry.chr);
            gt.pos.push_back(entry.pos);
        }

        return 0;
    }

    int write_ped_map(const Genotype &gt, const std::string &filename)
    {
        CFile file(filename, "w");

        if (!file) {
            eprint("ERROR: can't open file for writing: %s\n", filename);
            return 1;
        }

        isize_t m = length(gt.loc);

        for (isize_t j = 0; j < m; ++j)
            fprint(file, "%s %s 0 %d\n", gt.chr[j], gt.loc[j], gt.pos[j]);

        return 0;
    }

} // namespace

int parse_ped_entry(const std::string &s, PedEntry &entry)
{
    std::vector<Token> v;
    split(s, " \t", v);

    if (length(v) < 6) {
        eprint("ERROR: incorrect number of columns at PLINK PED entry line: %td\n", length(v));
        return 1;
    }

    entry.fid = v[0].str();
    entry.iid = v[1].str();
    entry.pid = v[2].str();
    entry.mid = v[3].str();
    entry.sex = std::stoi(v[4].str());
    entry.pheno = std::stoi(v[5].str());

    entry.gt.clear();
    for (auto itr = v.begin() + 6; itr != v.end(); ++itr) {
        char a;
        if (parse_ped_gt(itr->data(), length(*itr), a) != 0)
            return 1;
        entry.gt.push_back(static_cast<allele_t>(a));
    }

    return 0;
}

int parse_map_entry(const std::string &s, MapEntry &entry)
{
    auto v = split(s, " \t");

    if (length(v) != 4) {
        eprint("ERROR: expected 4 columns at PLINK MAP entry line: %td\n", length(v));
        return 1;
    }

    entry.chr = v[0];
    entry.id = v[1];
    entry.dist = std::stod(v[2]);
    entry.pos = std::stoi(v[3]);

    return 0;
}

int read_ped(const std::string &filename, Genotype &gt)
{
    if (read_ped_map(filename + ".map", gt) != 0)
        return 1;

    CFileLineReader file(filename + ".ped");

    if (!file) {
        eprint("ERROR: can't open file for reading: %s.ped\n", filename);
        return 1;
    }

    PedEntry entry;

    std::vector<std::string> fid, iid;
    std::vector< std::vector<allele_t> > dat;

    for (std::string line; file.read(line); ) {
        if (parse_ped_entry(line, entry) != 0)
            return 1;

        if (length(entry.gt) != 2 * length(gt.loc)) {
            eprint("ERROR: column count doesn't match at %s %s\n", entry.fid, entry.iid);
            return 1;
        }

        fid.push_back(entry.fid);
        iid.push_back(entry.iid);
        dat.push_back(entry.gt);
    }

    if (has_duplicate(iid)) {
        isize_t n = length(iid);
        for (isize_t i = 0; i < n; ++i)
            fid[i].append("_").append(iid[i]);
        gt.ind.swap(fid);
    }
    else
        gt.ind.swap(iid);

    isize_t m = length(gt.loc);
    isize_t n = length(dat);
    std::vector<std::string> as;

    for (isize_t j = 0; j < m; ++j) {
        std::vector<allele_t> v;
        for (isize_t i = 0; i < n; ++i) {
            v.push_back(dat[i][j*2]);
            v.push_back(dat[i][j*2+1]);
        }

        auto z = v;
        std::sort(z.begin(), z.end());
        z.erase(std::unique(z.begin(), std::remove(z.begin(), z.end(), allele_t('N'))), z.end());

        isize_t na = length(z);
        if (na > 2) {
            eprint("ERROR: PED variant must be bi-allelic: %c", z[0]);
            for (isize_t i = 1; i < na; ++i)
                eprint("/%c", z[i]);
            eprint("\n");
            return 1;
        }

        as.clear();
        allele_t mis = 'N', ref = 'N', alt = 'N';
        if (!z.empty()) {
            ref = z[0];
            as.emplace_back(1, z[0]);
        }
        if (na == 2) {
            alt = z[1];
            as.emplace_back(1, z[1]);
        }

        for (auto &a : v) {
            if (a == mis)
                a = 0;
            else if (a == ref)
                a = 1;
            else if (a == alt)
                a = 2;
        }

        gt.dat.push_back(v);
        gt.allele.push_back(as);
    }

    gt.ploidy = 2;
    gt.source = Genotype::Source::PED;

    return 0;
}

int write_ped(const Genotype &gt, const std::string &filename)
{
    int info = check_compat_ped(gt);

    if (info != 0) {
        eprint("ERROR: genotype data is not compatible with PED format: %d\n", info);
        return 1;
    }

    if (write_ped_map(gt, filename + ".map") != 0)
        return 2;

    CFile file(filename + ".ped", "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s.ped\n", filename);
        return 3;
    }

    isize_t m = length(gt.loc);
    isize_t n = length(gt.ind);
    bool haploid = gt.ploidy != 2;

    std::vector<std::string> fid, iid;
    parse_fid_iid(gt.ind, fid, iid);

    std::string line;

    for (isize_t i = 0; i < n; ++i) {
        line.clear();
        isize_t k1 = haploid ? i : i * 2;
        isize_t k2 = haploid ? i : i * 2 + 1;
        for (isize_t j = 0; j < m; ++j) {
            auto a = gt.dat[j][k1];
            auto b = gt.dat[j][k2];
            if (a && b) {
                line.push_back(' ');
                line.append(gt.allele[j][a-1]);
                line.push_back(' ');
                line.append(gt.allele[j][b-1]);
            }
            else
                line.append(" 0 0");
        }
        fprint(file, "%s %s 0 0 1 0", fid[i], iid[i]);
        fprint(file, "%s\n", line);
    }

    return 0;
}
