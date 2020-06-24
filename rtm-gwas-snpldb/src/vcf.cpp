#include "vcf.h"

#include <algorithm>

#include "cfile.h"
#include "print.h"
#include "stringutil.h"
#include "vectorutil.h"

#ifdef HAVE_ZLIB
#include "gzipfile.h"
#endif // HAVE_ZLIB

namespace {

    // parse GT, return ploidy number (1 or 2), otherwise error
    int parse_vcf_gt(const char *s, isize_t n, int &a, int &b)
    {
        auto beg = s;
        auto end = std::find(s, s + n, ':');
        isize_t len = end - beg;

        a = b = -1;

        if (len == 0)
            return -1;

        // haploid, single char: i, [0,9]
        if (len == 1) {
            if (*beg == '.')
                return 1;
            a = *beg - '0';
            return a < 0 || a > 9 ? -2 : 1;
        }

        // diploid, single char: i/j or i|j, [0,9]
        if (len == 3 && (beg[1] == '/' || beg[1] == '|')) {
            if (beg[0] != '.') {
                a = beg[0] - '0';
                if (a < 0 || a > 9)
                    return -3;
            }

            if (beg[2] != '.') {
                b = beg[2] - '0';
                if (b < 0 || b > 9)
                    return -4;
            }

            return 2;
        }

        // multiple char

        isize_t nsep = std::count(beg, end, '/') + std::count(beg, end, '|');

        if (nsep > 1) {
            eprint("ERROR: polyploidy genotype is not supported: %s\n", std::string(beg, end));
            return -5;
        }

        if (nsep == 0) {
            a = std::stoi(std::string(beg, end));
            return a < 0 ? -6 : 1;
        }

        isize_t pos = 0;
        while (pos < len && beg[pos] != '/' && beg[pos] != '|')
            ++pos;
        if (pos == 0 || pos == len - 1)
            return -7;

        std::string gt;
        gt.reserve(len);

        gt.assign(beg, beg + pos);
        if (gt != ".") {
            a = std::stoi(gt);
            if (a < 0)
                return -8;
        }

        gt.assign(beg + pos + 1, end);
        if (gt != ".") {
            b = std::stoi(gt);
            if (b < 0)
                return -9;
        }

        return 2;
    }

} // namespace

int parse_vcf_header(const std::string &s, std::vector<std::string> &v)
{
    v.clear();
    split(s, "\t", v);

    isize_t n = length(v);

    if (n != 8 && n < 10) {
        eprint("ERROR: incorrect number of columns in VCF header: %td\n", length(v));
        return 1;
    }

    static const char * const z[8] = { "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO" };

    for (int i = 0; i < 8; ++i) {
        if (v[i] != z[i]) {
            eprint("ERROR: incorrect column name in VCF header: %s\n", v[i]);
            return 1;
        }
    }

    if (n > 9) {
        if (v[8] != "FORMAT") {
            eprint("ERROR: FORMAT is required at 9th field in VCF header: %s\n", v[8]);
            return 1;
        }
        v.erase(v.begin(), v.begin() + 9);
    }
    else
        v.clear();

    return 0;
}

int parse_vcf_entry(const std::string &s, VcfEntry &entry)
{
    std::vector<Token> v;
    split(s, "\t", v);

    isize_t n = length(v);

    if (n != 8 && n < 10) {
        eprint("ERROR: incorrect number of columns at VCF entry: %td\n", n);
        return 1;
    }

    entry.chr = v[0].str();
    entry.pos = std::stoi(v[1].str());

    if (length(v[2]) == 1 && v[2][0] == '.')
        entry.id = v[0].str() + "_" + v[1].str();
    else
        entry.id = v[2].str();

    entry.as.clear();
    entry.as.push_back(v[3].str());
    split(v[4].str(), ",", entry.as);

    if (n == 8)
        return 0;

    if (length(v[8]) < 2 || v[8][0] != 'G' || v[8][1] != 'T') {
        eprint("ERROR: GT is required and must be the first sub-field of FORMAT: %s\n", v[8].str());
        return 1;
    }

    isize_t na = length(entry.as);
    if (na > 255) {
        eprint("ERROR: exceed the maximum number of alleles (255): %td\n", na);
        return 1;
    }

    entry.gt.clear();

    for (isize_t i = 9; i < n; ++i) {
        int a = -1, b = -1;
        int info = parse_vcf_gt(v[i].data(), length(v[i]), a, b);

        if ((info != 1 && info != 2) || a >= na || b >= na) {
            eprint("ERROR: invalid genotype data: %s %d %s %s %s\n",
                   entry.chr, entry.pos, v[3].str(), v[4].str(), v[i].str());
            return 1;
        }

        if (i == 9)
            entry.ploidy = info;

        if (info != entry.ploidy) {
            eprint("ERROR: genotype ploidy doesn't match (%d != %d): %s\n", entry.ploidy, info, v[i].str());
            return 1;
        }

        entry.gt.push_back(a < 0 ? 0 : (allele_t) (a+1));
        if (entry.ploidy == 2)
            entry.gt.push_back(b < 0 ? 0 : (allele_t) (b+1));
    }

    return 0;
}

int read_vcf(const std::string &filename, Genotype &gt)
{
    CFileLineReader file(filename);

    if (!file) {
        eprint("ERROR: can't open file for reading: %s\n", filename);
        return 1;
    }

    isize_t ln = 0;

    for (std::string line; file.read(line); ) {
        ++ln;

        if (line.compare(0, 2, "##") == 0)
            continue;

        if (line.compare(0, 1, "#") == 0) {
            if (parse_vcf_header(line, gt.ind) != 0)
                return 1;
            break;
        }

        eprint("ERROR: VCF header line is required\n");
        return 1;
    }

    VcfEntry entry;

    for (std::string line; file.read(line); ) {
        ++ln;

        if (parse_vcf_entry(line, entry) != 0)
            return 1;

        if (gt.ploidy != 1 && gt.ploidy != 2)
            gt.ploidy = entry.ploidy;

        if (entry.ploidy > gt.ploidy) {
            eprint("ERROR: genotype ploidy doesn't match at line %td\n", ln);
            return 1;
        }

        if (entry.ploidy < gt.ploidy) {
            std::vector<allele_t> v;
            v.reserve(length(entry.gt) * 2);
            for (auto a : entry.gt) {
                v.push_back(a);
                v.push_back(a);
            }
            entry.gt.swap(v);
            entry.ploidy = gt.ploidy;
        }

        if (length(entry.gt) != entry.ploidy * length(gt.ind)) {
            eprint("ERROR: column count doesn't match at line %td\n", ln);
            return 1;
        }

        gt.loc.push_back(entry.id);
        gt.chr.push_back(entry.chr);
        gt.pos.push_back(entry.pos);
        gt.allele.push_back(entry.as);
        gt.dat.push_back(entry.gt);
    }

    gt.source = Genotype::Source::VCF;

    return 0;
}

int read_vcf_gz(const std::string &filename, Genotype &gt)
{
#ifdef HAVE_ZLIB

    GzipFileLineReader file(filename);

    if (!file) {
        eprint("ERROR: can't open file for reading: %s\n", filename);
        return 1;
    }

    isize_t ln = 0;

    for (std::string line; file.read(line); ) {
        ++ln;

        if (line.compare(0, 2, "##") == 0)
            continue;

        if (line.compare(0, 1, "#") == 0) {
            if (parse_vcf_header(line, gt.ind) != 0)
                return 1;
            break;
        }

        eprint("ERROR: VCF header line is required\n");
        return 1;
    }

    VcfEntry entry;

    for (std::string line; file.read(line); ) {
        ++ln;

        if (parse_vcf_entry(line, entry) != 0)
            return 1;

        if (gt.ploidy != 1 && gt.ploidy != 2)
            gt.ploidy = entry.ploidy;

        if (entry.ploidy > gt.ploidy) {
            eprint("ERROR: genotype ploidy doesn't match at line %td\n", ln);
            return 1;
        }

        if (entry.ploidy < gt.ploidy) {
            std::vector<allele_t> v;
            v.reserve(length(entry.gt) * 2);
            for (auto a : entry.gt) {
                v.push_back(a);
                v.push_back(a);
            }
            entry.gt.swap(v);
            entry.ploidy = gt.ploidy;
        }

        if (length(entry.gt) != entry.ploidy * length(gt.ind)) {
            eprint("ERROR: column count doesn't match at line %td\n", ln);
            return 1;
        }

        gt.loc.push_back(entry.id);
        gt.chr.push_back(entry.chr);
        gt.pos.push_back(entry.pos);
        gt.allele.push_back(entry.as);
        gt.dat.push_back(entry.gt);
    }

    gt.source = Genotype::Source::VCF;

    return 0;

#else

    eprint("ERROR: Gzip compressed VCF file (.vcf.gz) is not supported: %s\n", filename);
    return 1;

#endif // HAVE_ZLIB
}

int write_vcf(const Genotype &gt, const std::string &filename, bool force_diploid)
{
    CFile file(filename, "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s\n", filename);
        return 1;
    }

    std::vector<std::string> codes(256);
    codes[0] = ".";
    for (int i = 1; i < 256; ++i)
        codes[i] = std::to_string(i-1);

    isize_t n = length(gt.ind);

    fprint(file, "##fileformat=VCFv4.2\n");
    fprint(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if (!gt.ind.empty()) {
        fprint(file, "\tFORMAT");
        for (isize_t i = 0; i < n; ++i)
            fprint(file, "\t%s", gt.ind[i]);
    }
    fprint(file, "\n");

    std::string line;
    isize_t m = length(gt.loc);
    bool haploid = gt.ploidy != 2;

    for (isize_t j = 0; j < m; ++j) {
        fprint(file, "%s\t%d\t%s\t", gt.chr[j], gt.pos[j], gt.loc[j]);

        if (gt.allele[j].empty())
            fprint(file, ".\t.");
        else {
            isize_t na = length(gt.allele[j]);
            fprint(file, "%s\t", gt.allele[j][0]);
            if (na == 1)
                fprint(file, ".");
            else {
                fprint(file, "%s", gt.allele[j][1]);
                for (isize_t k = 2; k < na; ++k)
                    fprint(file, ",%s", gt.allele[j][k]);
            }
        }

        fprint(file, "\t.\t.\t.");

        if (gt.ind.empty()) {
            fprint(file, "\n");
            continue;
        }

        line.clear();
        line.append("\tGT");

        if (haploid) {
            for (isize_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i];
                line.append("\t").append(codes[a]);
                if (force_diploid)
                    line.append("/").append(codes[a]);
            }
        }
        else {
            for (isize_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i*2], b = gt.dat[j][i*2+1];
                line.append("\t").append(codes[a]).append("/").append(codes[b]);
            }
        }

        fprint(file, "%s\n", line);
    }

    return 0;
}

int write_vcf_gz(const Genotype &gt, const std::string &filename, bool force_diploid)
{
#ifdef HAVE_ZLIB

    GzipFile file(filename, "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s\n", filename);
        return 1;
    }

    std::vector<std::string> codes(256);
    codes[0] = ".";
    for (int i = 1; i < 256; ++i)
        codes[i] = std::to_string(i-1);

    isize_t n = length(gt.ind);

    gzputs(file, "##fileformat=VCFv4.2\n");
    gzputs(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if (!gt.ind.empty()) {
        gzputs(file, "\tFORMAT");
        for (isize_t i = 0; i < n; ++i) {
            gzputc(file, '\t');
            gzputs(file, gt.ind[i].c_str());
        }
    }
    gzputc(file, '\n');

    std::string line;
    isize_t m = length(gt.loc);
    bool haploid = gt.ploidy != 2;

    for (isize_t j = 0; j < m; ++j) {
        gzprintf(file, "%s\t%d\t%s\t", gt.chr[j].c_str(), gt.pos[j], gt.loc[j].c_str());

        if (gt.allele[j].empty())
            gzputs(file, ".\t.");
        else {
            isize_t na = length(gt.allele[j]);
            gzputs(file, gt.allele[j][0].c_str());
            gzputc(file, '\t');
            if (na == 1) {
                gzputc(file, '.');
            }
            else {
                gzputs(file, gt.allele[j][1].c_str());
                for (isize_t k = 2; k < na; ++k) {
                    gzputc(file, ',');
                    gzputs(file, gt.allele[j][k].c_str());
                }
            }
        }

        gzputs(file, "\t.\t.\t.");

        if (gt.ind.empty()) {
            gzputc(file, '\n');
            continue;
        }

        line.clear();
        line.append("\tGT");

        if (haploid) {
            for (isize_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i];
                line.append("\t").append(codes[a]);
                if (force_diploid)
                    line.append("/").append(codes[a]);
            }
        }
        else {
            for (isize_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i*2], b = gt.dat[j][i*2+1];
                line.append("\t").append(codes[a]).append("/").append(codes[b]);
            }
        }

        gzputs(file, line.c_str());
        gzputc(file, '\n');
    }

    return 0;

#else

    eprint("ERROR: Gzip compressed VCF file (.vcf.gz) is not supported: %s\n", filename);
    return 1;

#endif // HAVE_ZLIB
}
