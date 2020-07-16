#include "hmp.h"

#include <cstring>
#include <algorithm>

#include "cfile.h"
#include "print.h"
#include "stringutil.h"
#include "vectorutil.h"

namespace {

    int parse_hmp_gt(const char *s, isize_t n, char &a, char &b)
    {
        static const char gs[] = { 'A','C','G','T','I','D','N' };

        if (n != 2) {
            eprint("ERROR: HapMap genotype must be represented by 2 characters: %s\n", std::string(s, n));
            return 1;
        }

        a = s[0];
        b = s[1];

        if (std::memchr(gs, a, sizeof gs) == NULL ||
            std::memchr(gs, b, sizeof gs) == NULL) {
            eprint("ERROR: invalid HapMap genotype: %c%c\n", a, b);
            return 2;
        }

        return 0;
    }

    int check_compat_hmp(const Genotype &gt)
    {
        if (gt.source == Genotype::Source::HAPMAP)
            return 0;

        if (gt.source == Genotype::Source::PED)
            return 0;

        if (gt.source == Genotype::Source::SNPLDB)
            return 1;

        static const char gs[] = { 'A','C','G','T','-','I','D' };

        for (auto &v : gt.allele) {
            if (v.size() > 2)
                return 2;
            for (auto &e : v) {
                if (e.size() == 1 && std::memchr(gs, e[0], sizeof gs) == NULL)
                    return 3;
                if (e.find_first_not_of("ACGT") != std::string::npos)
                    return 4;
            }
        }

        return 0;
    }

} // namespace


int parse_hmp_header(const std::string &s, std::vector<std::string> &v)
{
    v.clear();
    split(s, " \t", v);

    if (v.size() < 11) {
        eprint("ERROR: incorrect number of columns in HapMap header: %td\n", length(v));
        return 1;
    }

    v.erase(v.begin(), v.begin() + 11);

    return 0;
}

int parse_hmp_entry(const std::string &s, HmpEntry &entry)
{
    std::vector<Token> v;
    split(s, " \t", v);

    if (length(v) < 11) {
        eprint("ERROR: incorrect number of columns in HapMap entry: %td\n", length(v));
        return 1;
    }

    entry.id = v[0].str();
    entry.chr = v[2].str();
    entry.pos = std::stoi(v[3].str());

    entry.gt.clear();
    for (auto itr = v.begin() + 11; itr != v.end(); ++itr) {
        char a = 0, b = 0;
        if (parse_hmp_gt(itr->data(), length(*itr), a, b) != 0)
            return 1;
        entry.gt.push_back((allele_t) a);
        entry.gt.push_back((allele_t) b);
    }

    // A, C, G, T, -, AG, ...
    entry.as.clear();
    split(v[1].str(), "/", entry.as);

    if (length(entry.as) != 2 || entry.as[0] == "N" || entry.as[1] == "N") {
        auto z = entry.gt;
        std::sort(z.begin(), z.end());
        z.erase(std::unique(z.begin(), std::remove(z.begin(), z.end(), 'N')), z.end());

        if (length(z) > 2) {
            eprint("ERROR: HapMap variant must be bi-allelic:");
            for (auto a : z)
                eprint(" %c", a);
            eprint("\n");
            return 1;
        }

        // A, C, G, T, I, D
        entry.as.clear();
        if (!z.empty())
            entry.as.emplace_back(1, z[0]);
        if (length(z) == 2)
            entry.as.emplace_back(1, z[1]);
    }

    if (entry.as.empty()) {
        std::fill(entry.gt.begin(), entry.gt.end(), (allele_t) 0);
        return 0;
    }

    if (length(entry.as) == 1) {
        std::fill(entry.gt.begin(), entry.gt.end(), (allele_t) 1);
        return 0;
    }

    allele_t ref = 'N', alt = 'N';
    if (entry.as[0] == "A" || entry.as[0] == "C" || entry.as[0] == "G" || entry.as[0] == "T")
        ref = (allele_t) entry.as[0][0];
    if (entry.as[1] == "A" || entry.as[1] == "C" || entry.as[1] == "G" || entry.as[1] == "T")
        alt = (allele_t) entry.as[1][0];

    static const allele_t mis = 'N', ins = 'I', del = 'D';

    for (auto &a : entry.gt) {
        if (a == mis)
            a = (allele_t) 0;
        else if (a == ref)
            a = (allele_t) 1;
        else if (a == alt)
            a = (allele_t) 2;
        else if (a == ins)
            a = entry.as[0] == "-" ? 2 : 1;
        else if (a == del)
            a = entry.as[0] == "-" ? 1 : 2;
        else {
            eprint("ERROR: inconsistent allele code: %s, %s, %c\n", entry.id, v[1].str(), a);
            return 1;
        }
    }

    return 0;
}

int read_hmp(const std::string &filename, Genotype &gt)
{
    CFileLineReader file(filename);

    if (!file) {
        eprint("ERROR: can't open file for reading: %s\n", filename);
        return 1;
    }

    for (std::string line; file.read(line); ) {
        if (parse_hmp_header(line, gt.ind) != 0)
            return 1;
        break;
    }

    HmpEntry entry;

    for (std::string line; file.read(line); ) {
        if (parse_hmp_entry(line, entry) != 0)
            return 1;

        if (length(entry.gt) != 2 * length(gt.ind)) {
            eprint("ERROR: column count doesn't match at %s\n", entry.id);
            return 1;
        }

        gt.loc.push_back(entry.id);
        gt.chr.push_back(entry.chr);
        gt.pos.push_back(entry.pos);
        gt.allele.push_back(entry.as);
        gt.dat.push_back(entry.gt);
    }

    gt.ploidy = 2;
    gt.source = Genotype::Source::HAPMAP;

    return 0;
}

int write_hmp(const Genotype &gt, const std::string &filename)
{
    int info = check_compat_hmp(gt);

    if (info != 0) {
        eprint("ERROR: genotype data is not compatible with HapMap format: %d\n", info);
        return 1;
    }

    CFile file(filename, "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s\n", filename);
        return 1;
    }

    isize_t m = length(gt.loc);
    isize_t n = length(gt.ind);
    bool haploid = gt.ploidy != 2;

    std::vector<char> indel(m, 0);
    if (gt.source == Genotype::Source::HAPMAP) {
        for (isize_t j = 0; j < m; ++j) {
            auto& as = gt.allele[j];
            if (length(as) == 1 && (as[0] == "-" || length(as[0]) > 1))
                indel[j] = 1;
            if (length(as) == 2 && (as[0] == "-" || as[1] == "-"))
                indel[j] = 1;
        }
    }

    fprint(file, "rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode");
    for (isize_t i = 0; i < n; ++i)
        fprint(file, " %s", gt.ind[i]);
    fprint(file, "\n");

    std::string line;

    for (isize_t j = 0; j < m; ++j) {
        line.clear();

        if (indel[j]) {
            char ref = 'D', alt = 'I';
            if (gt.allele[j][0] != "-")
                std::swap(ref, alt);
            for (isize_t i = 0; i < n; ++i) {
                auto a = haploid ? gt.dat[j][i] : gt.dat[j][i*2];
                auto b = haploid ? gt.dat[j][i] : gt.dat[j][i*2+1];
                if (a && b) {
                    line.push_back(' ');
                    line.push_back(a == 1 ? ref : alt);
                    line.push_back(b == 1 ? ref : alt);
                }
                else
                    line.append(" NN");
            }
        }
        else {
            for (isize_t i = 0; i < n; ++i) {
                auto a = haploid ? gt.dat[j][i] : gt.dat[j][i*2];
                auto b = haploid ? gt.dat[j][i] : gt.dat[j][i*2+1];
                if (a && b) {
                    line.push_back(' ');
                    line.append(gt.allele[j][a-1]).append(gt.allele[j][b-1]);
                }
                else
                    line.append(" NN");
            }
        }

        fprint(file, "%s ", gt.loc[j]);

        if (gt.allele[j].empty())
            fprint(file, "N/N");
        else {
            fprint(file, "%s", gt.allele[j][0]);
            if (length(gt.allele[j]) == 1)
                fprint(file, "/N");
            else
                fprint(file, "/%s", gt.allele[j][1]);
        }

        fprint(file, " %s %d + NA NA NA NA NA NA", gt.chr[j], gt.pos[j]);
        fprint(file, "%s\n", line);
    }

    return 0;
}
