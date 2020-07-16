#include "rtmgwasld.h"

#include <cmath>
#include <limits>
#include <algorithm>

#include "types.h"
#include "print.h"
#include "cmdline.h"
#include "vcf.h"
#include "cfile.h"
#include "stringutil.h"
#include "vectorutil.h"
#include "popgen.h"

#ifndef RTM_GWAS_VERSION
#define RTM_GWAS_VERSION "unknown"
#endif // RTM_GWAS_VERSION

namespace {

    double kNaN = std::numeric_limits<double>::quiet_NaN();

    struct Parameter
    {
        std::string vcf;
        std::string out;
        std::string loc;
        double rsq = 0.5;
        int maxdist = 500000;
    } par;

    void show_help()
    {
        const char *help =
            "Usage: rtm-gwas-ld [options...]\n"
            "Options:\n"
            "  --vcf <>            SNP genotype data file in VCF format\n"
            "  --out <ld.out>      output file prefix\n"
            "  --maxdist <500000>  maximum inter-locus distance (bp)\n"
            "  --with-loc <>       locus list (locus names) file\n"
            "  --min-r2 <0.5>      minimum r2 threshold for --with-loc output\n"
            "\n";
        eprint(help);
    }

    int parse_cmdline(int argc, char *argv[])
    {
        CmdLine cmd;

        cmd.add("--vcf", "");
        cmd.add("--out", "ld.out");
        cmd.add("--maxdist", "500000");
        cmd.add("--with-loc", "");
        cmd.add("--min-r2", "0.5");

        if (cmd.parse(argc, argv) != 0)
            return 1;

        par.vcf = cmd.get("--vcf");
        par.out = cmd.get("--out");
        par.maxdist = std::stoi(cmd.get("--maxdist"));
        par.loc = cmd.get("--with-loc");
        par.rsq = std::stod(cmd.get("--min-r2"));

        return 0;
    }

    void calc_dprime_rsq_hap_snp(const std::vector<allele_t> &x, const std::vector<allele_t> &y, double &dprime, double &rsq)
    {
        dprime = rsq = kNaN;

        static const allele_t a = 1, b = 1;

        isize_t n = length(x);
        isize_t nt = 0, na = 0, nb = 0, nab = 0;

        for (isize_t i = 0; i < n; ++i) {
            if (x[i] && y[i]) {
                ++nt;
                if (x[i] == a) {
                    ++na;
                    if (y[i] == b)
                        ++nab;
                }
                if (y[i] == b)
                    ++nb;
            }
        }

        if (nt < 2)
            return;

        double fnt = nt;

        calc_dprime_rsq(na/fnt, nb/fnt, nab/fnt, dprime, rsq);

        dprime = std::fabs(dprime);
    }

    void calc_dprime_rsq_dip_snp(const std::vector<allele_t> &x, const std::vector<allele_t> &y, double &dprime, double &rsq)
    {
        dprime = rsq = kNaN;

        isize_t nt = 0;
        int g[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

        isize_t n = length(x) / 2;

        for (isize_t i = 0; i < n; ++i) {
            isize_t i1 = i*2;
            isize_t i2 = i*2+1;
            if (x[i1] && x[i2] && y[i1] && y[i2]) {
                int a = x[i1] + x[i2] - 2;
                int b = y[i1] + y[i2] - 2;
                ++nt;
                ++g[a][b];
            }
        }

        if (nt < 2)
            return;

        double p[4];
        est_hap_freq(1000, 1e-5, g, p);

        calc_dprime_rsq(p[0] + p[1], p[0] + p[2], p[0], dprime, rsq);

        dprime = std::fabs(dprime);
    }

    void calc_dprime_rsq_hap(const std::vector<allele_t> &x, const std::vector<allele_t> &y, double &dprime, double &rsq)
    {
        dprime = rsq = kNaN;

        isize_t n = length(x);

        allele_t xmax = * std::max_element(x.begin(), x.end());
        allele_t ymax = * std::max_element(y.begin(), y.end());

        if (xmax < 2 || ymax < 2)
            return;

        if (xmax == 2 && ymax == 2) {
            calc_dprime_rsq_hap_snp(x, y, dprime, rsq);
            return;
        }

        std::vector<isize_t> idx;
        idx.reserve(n);
        for (isize_t i = 0; i < n; ++i) {
            if (x[i] && y[i])
                idx.push_back(i);
        }

        isize_t nt = length(idx);
        if (nt < 2)
            return;

        dprime = rsq = 0.0;

        for (allele_t a = 1; a <= xmax; ++a) {
            for (allele_t b = 1; b <= ymax; ++b) {
                isize_t na = 0, nb = 0, nab = 0;

                for (auto i : idx) {
                    if (x[i] == a) {
                        ++na;
                        if (y[i] == b)
                            ++nab;
                    }
                    if (y[i] == b)
                        ++nb;
                }

                if (na == 0 || nb == 0)
                    continue;

                double fnt = nt;
                double pa = na / fnt;
                double pb = nb / fnt;
                double pab = nab / fnt;

                double dprime_ab = 0.0, rsq_ab = 0.0;
                calc_dprime_rsq(pa, pb, pab, dprime_ab, rsq_ab);

                // Farnir F, Coppieters W, Arranz JJ, et al. Extensive genome-wide linkage disequilibrium in cattle.
                //    Genome Res. 2000, 10(2):220-227. doi:10.1101/gr.10.2.220

                pab = pa * pb;
                dprime += pab * std::fabs(dprime_ab);
                rsq += pab * rsq_ab;
            }
        }
    }

    // assuming phased
    void calc_dprime_rsq_dip(const std::vector<allele_t> &x, const std::vector<allele_t> &y, double &dprime, double &rsq)
    {
        dprime = rsq = kNaN;

        isize_t n = length(x) / 2;

        allele_t xmax = * std::max_element(x.begin(), x.end());
        allele_t ymax = * std::max_element(y.begin(), y.end());

        if (xmax < 2 || ymax < 2)
            return;

        if (xmax == 2 && ymax == 2) {
            calc_dprime_rsq_dip_snp(x, y, dprime, rsq);
            return;
        }

        std::vector<isize_t> idx;
        idx.reserve(n);
        for (isize_t i = 0; i < n; ++i) {
            isize_t i1 = i*2;
            isize_t i2 = i*2+1;
            if (x[i1] && x[i2] && y[i1] && y[i2])
                idx.push_back(i);
        }

        isize_t nt = length(idx);
        if (nt < 2)
            return;

        dprime = rsq = 0.0;

        for (allele_t a = 1; a <= xmax; ++a) {
            for (allele_t b = 1; b <= ymax; ++b) {
                int na = 0, nb = 0, nab = 0;

                for (auto i : idx) {
                    isize_t i1 = i*2;
                    isize_t i2 = i*2+1;
                    if (x[i1] == a) {
                        ++na;
                        if (y[i1] == b)
                            ++nab;
                    }

                    if (x[i2] == a) {
                        ++na;
                        if (y[i2] == b)
                            ++nab;
                    }

                    if (y[i1] == b)
                        ++nb;

                    if (y[i2] == b)
                        ++nb;
                }

                if (na == 0 || nb == 0)
                    continue;

                double fnt = 2*nt;
                double pa = na / fnt;
                double pb = nb / fnt;
                double pab = nab / fnt;

                double dprime_ab = 0.0, rsq_ab = 0.0;
                calc_dprime_rsq(pa, pb, pab, dprime_ab, rsq_ab);

                pab = pa * pb;
                dprime += pab * std::fabs(dprime_ab);
                rsq += pab * rsq_ab;
            }
        }
    }

    int read_locus_list(const std::string &filename, std::vector<std::string> &loc)
    {
        std::string text;
        if (read_all_text(filename, text) != 0)
            return 1;

        split(text, " \t\r\n", loc);

        return 0;
    }

    int ld_with_loc(const Genotype &gt)
    {
        std::vector<std::string> loc;
        if (read_locus_list(par.loc, loc) != 0)
            return 1;

        CFile file(par.out + ".with", "w");

        if (!file) {
            eprint("ERROR: can't open file: %s.with\n", par.out);
            return 1;
        }

        fprint(file, "Locus1\tChromosome1\tPosition1\tLocus2\tChromosome2\tPosition2\tDPrime\tRSquare\n");

        isize_t m = length(gt.loc);
        bool diploid = gt.ploidy == 2;

        for (auto curr_loc : loc) {
            isize_t i = index(gt.loc, curr_loc);
            if (i == -1) {
                eprint("WARNING: can't find locus: %\n", curr_loc);
                continue;
            }

            eprint("INFO: locus %s\n", curr_loc);

            for (isize_t j = 0; j < m; ++j) {
                if (j == i)
                    continue;

                isize_t n1 = length(gt.allele[i]);
                isize_t n2 = length(gt.allele[j]);

                if (n1 < 2 || n2 < 2)
                    continue;

                double dprime = kNaN;
                double rsq = kNaN;

                if (n1 == 2 && n2 == 2) {
                    if (diploid)
                        calc_dprime_rsq_dip_snp(gt.dat[i], gt.dat[j], dprime, rsq);
                    else
                        calc_dprime_rsq_hap_snp(gt.dat[i], gt.dat[j], dprime, rsq);
                }
                else {
                    if (diploid)
                        calc_dprime_rsq_dip(gt.dat[i], gt.dat[j], dprime, rsq);
                    else
                        calc_dprime_rsq_hap(gt.dat[i], gt.dat[j], dprime, rsq);
                }

                if (!std::isfinite(dprime) || !std::isfinite(rsq))
                    continue;

                if (rsq < par.rsq)
                    continue;

                fprint(file, "%s\t%s\t%d\t%s\t%s\t%d\t%g\t%g\n",
                       gt.loc[i], gt.chr[i], gt.pos[i], gt.loc[j], gt.chr[j], gt.pos[j], dprime, rsq);
            }
        }

        eprint("INFO: LD has finished successfully\n");

        return 0;
    }

} // namespace

int rtm_gwas_ld(int argc, char *argv[])
{
    eprint("RTM-GWAS %s LD (Built on %s %s)\n", RTM_GWAS_VERSION, __DATE__, __TIME__);

    if (argc < 2) {
        show_help();
        return 1;
    }

    if (parse_cmdline(argc, argv) != 0)
        return 1;

    Genotype gt;
    eprint("INFO: reading genotype file...\n");
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    eprint("INFO: %td individuals, %td loci\n", length(gt.ind), length(gt.loc));

    if (gt.ploidy != 1 && gt.ploidy != 2) {
        eprint("ERROR: unsupported ploidy of genotype data: %d\n", gt.ploidy);
        return 1;
    }

    if (!par.loc.empty())
        return ld_with_loc(gt);

    CFile file(par.out + ".all", "w");

    if (!file) {
        eprint("ERROR: can't open file: %s.all\n", par.out);
        return 1;
    }

    fprint(file, "Chromosome\tLocus1\tPosition1\tLocus2\tPosition2\tDistance\tDPrime\tRSquare\n");

    int n = par.maxdist / 1000;
    std::vector<isize_t> count(n, 0);
    std::vector<double> rsq_mean(n, 0.0);
    std::vector<double> dprime_mean(n, 0.0);

    isize_t m = length(gt.loc);

    for (isize_t i = 0; i < m; ++i) {
        for (isize_t j = i + 1; j < m; ++j) {
            if (gt.chr[j] != gt.chr[i])
                continue;

            int dist = std::abs(gt.pos[j] - gt.pos[i]);
            if (dist > par.maxdist || dist == 0)
                continue;

            isize_t n1 = length(gt.allele[i]);
            isize_t n2 = length(gt.allele[j]);

            if (n1 < 2 || n2 < 2)
                continue;

            double dprime = kNaN;
            double rsq = kNaN;

            if (n1 == 2 && n2 == 2) {
                if (gt.ploidy == 1)
                    calc_dprime_rsq_hap_snp(gt.dat[i], gt.dat[j], dprime, rsq);
                else
                    calc_dprime_rsq_dip_snp(gt.dat[i], gt.dat[j], dprime, rsq);
            }
            else {
                if (gt.ploidy == 1)
                    calc_dprime_rsq_hap(gt.dat[i], gt.dat[j], dprime, rsq);
                else
                    calc_dprime_rsq_dip(gt.dat[i], gt.dat[j], dprime, rsq);
            }

            if (!std::isfinite(dprime) || !std::isfinite(rsq))
                continue;

            int k = dist / 1000;
            if (dist % 1000 == 0)
                --k;

            ++count[k];
            rsq_mean[k] += (rsq - rsq_mean[k]) / count[k];
            dprime_mean[k] += (dprime - dprime_mean[k]) / count[k];

            fprint(file, "%s\t%s\t%d\t%s\t%d\t%d\t%g\t%g\n",
                   gt.chr[i], gt.loc[i], gt.pos[i], gt.loc[j], gt.pos[j], dist, dprime, rsq);
        }
    }

    CFile file2(par.out + ".sum", "w");

    if (!file2) {
        eprint("ERROR: can't open file: %s.sum\n", par.out);
        return 1;
    }

    fprint(file2, "<=dist\tpairs\tr2\tdp\n");

    for (int i = 0; i < n; ++i) {
        if (count[i] < 1)
            continue;
        fprint(file2, "%d\t%td\t%g\t%g\n", (i+1)*1000, count[i], rsq_mean[i], dprime_mean[i]);
    }

    eprint("INFO: LD has finished successfully\n");

    return 0;
}
