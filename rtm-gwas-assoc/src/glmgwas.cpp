#include "glmgwas.h"

#include <limits>
#include <numeric>
#include <functional>
#include "print.h"
#include "cmdline.h"
#include "openmp.h"
#include "vectorutil.h"
#include "statutil.h"
#include "lmfit.h"
#include "pheno.h"
#include "cfile.h"

#ifndef GLM_GWAS_VERSION
#define GLM_GWAS_VERSION "unknown"
#endif // GLM_GWAS_VERSION

namespace {

    double kNaN = std::numeric_limits<double>::quiet_NaN();

    struct Parameter
    {
        std::string vcf;
        std::string pheno;
        std::string covar;
        std::string out;
        int thread = 0;
    } par;

    void show_help()
    {
        const char *help =
            "Usage: glm-gwas [options...]\n"
            "Options:\n"
            "  --vcf <>         genotype file in VCF format\n"
            "  --pheno <>       phenotype data file\n"
            "  --covar <>       covariate data file\n"
            "  --out <glm.out>  output file prefix\n"
            "  --thread <0>     set the number of threads\n"
            "\n";
        eprint(help);
    }

    void parse_cmdline(int argc, char *argv[])
    {
        CmdLine cmd;

        cmd.add("--vcf", "");
        cmd.add("--pheno", "");
        cmd.add("--covar", "");
        cmd.add("--out", "glm.out");
        cmd.add("--thread", "0");

        cmd.parse(argc, argv);

        par.vcf = cmd.get("--vcf");
        par.pheno = cmd.get("--pheno");
        par.covar = cmd.get("--covar");
        par.out = cmd.get("--out");
        par.thread = std::stoi(cmd.get("--thread"));
    }

    void filter_ind(const std::vector<isize_t> &idx, Phenotype &pt)
    {
        subset(pt.ind, idx).swap(pt.ind);

        if (!pt.env.empty())
            subset(pt.env, idx).swap(pt.env);

        if (!pt.blk.empty())
            subset(pt.blk, idx).swap(pt.blk);

        for (auto &e : pt.dat)
            subset(e, idx).swap(e);
    }

    void filter_ind(const std::vector<isize_t> &idx, Covariate &ct)
    {
        subset(ct.ind, idx).swap(ct.ind);

        for (auto &e : ct.dat)
            subset(e, idx).swap(e);
    }

} // namespace

void intersect(Genotype &gt, Phenotype &pt, Covariate &ct, std::vector<isize_t> &gi)
{
    bool require = gt.ind != pt.ind;
    bool covar = !ct.phe.empty() && !ct.ind.empty();

    if (!require && covar)
        require = gt.ind != ct.ind;

    if (!require) {
        gi.resize(gt.ind.size());
        std::iota(gi.begin(), gi.end(), 0);
        return;
    }

    eprint("INFO: performing data intersection by individual...\n");

    auto ind = intersect(gt.ind, pt.ind);
    if (covar)
        ind = intersect(ind, ct.ind);

    gi.clear();
    std::vector<isize_t> pi, ci;
    isize_t n = length(pt.ind);

    for (isize_t i = 0; i < n; ++i) {
        if (std::binary_search(ind.begin(), ind.end(), pt.ind[i])) {
            pi.push_back(i);
            gi.push_back(index(gt.ind, pt.ind[i]));
            if (covar)
                ci.push_back(index(ct.ind, pt.ind[i]));
        }
    }

    filter_ind(pi, pt);

    if (covar)
        filter_ind(ci, ct);

    eprint("INFO: there are %td individuals after intersection\n", length(ind));

    if (length(gt.ind) - length(ind) > 50)
        eprint("WARNING: the number of common individuals is much less than genotype\n");
}

void design(const Phenotype &pt, std::vector< std::vector<double> > &ac, std::vector< std::vector<double> > &ic)
{
    std::vector< std::vector<double> > xenv, xblk;

    if (!pt.env.empty()) {
        design2(grpidx(pt.env), xenv);
        ac.insert(ac.end(), xenv.begin(), xenv.end());
        ic.insert(ic.end(), xenv.begin(), xenv.end());
    }

    if (!pt.blk.empty()) {
        design2(grpidx(pt.blk), xblk);
        if (xenv.empty())
            ac.insert(ac.end(), xblk.begin(), xblk.end());
        else {
            design3(grpidx(pt.env), xenv);
            for (auto &x : xenv) {
                for (auto y : xblk) {
                    std::transform(x.begin(), x.end(), y.begin(), y.begin(), std::multiplies<double>());
                    ac.push_back(y);
                }
            }
        }
    }
}

int assoc_glm(const Genotype &gt, const std::vector<isize_t> &gi,
              const std::vector<double> &y,
              const std::vector< std::vector<double> > &ac,
              const std::vector< std::vector<double> > &ic,
              std::vector< std::vector<double> > &out)
{
    isize_t n = length(y);
    isize_t m = length(gt.dat);

    out.assign(ic.empty() ? 2 : 6, std::vector<double>(m, kNaN));

    double sst = css(n, y.data());

    isize_t p0 = 1 + length(ac);
    std::vector<double> x0(n, 1);
    for (auto &e : ac)
        x0.insert(x0.end(), e.begin(), e.end());

    LmFit lm;
    int info = lmfit(n, p0, y.data(), x0.data(), &lm);

    if (info != 0) {
        eprint("ERROR: lmfit failed in assoc_glm (%s:%d): %d\n", __FILE__, __LINE__, info);
        return 1;
    }

    double dfe0 = lm.dfe;
    double sse0 = lm.sse;

    bool haploid = gt.ploidy != 2;
    std::vector<allele_t> g1;
    std::vector< std::pair<allele_t, allele_t> > g2;

    std::vector<double> xv, yv;

    for (isize_t j = 0; j < m; ++j) {
        std::vector<isize_t> idx;
        std::vector< std::vector<double> > gx;

        if (haploid) {
            g1.clear();
            for (isize_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][gi[i]];
                if (a) {
                    g1.push_back(a);
                    idx.push_back(i);
                }
            }
            design2(grpidx(g1), gx);
        }
        else {
            g2.clear();
            for (isize_t i = 0; i < n; ++i) {
                isize_t i1 = gi[i] * 2;
                isize_t i2 = i1 + 1;
                auto a = gt.dat[j][i1];
                auto b = gt.dat[j][i2];
                if (a && b) {
                    if (a > b)
                        std::swap(a, b);
                    g2.emplace_back(a, b);
                    idx.push_back(i);
                }
            }
            design2(grpidx(g2), gx);
        }

        if (gx.empty())
            continue;

        isize_t nv = length(idx);

        double jsst = sst;
        double jdfe0 = dfe0;
        double jsse0 = sse0;

        if (nv == n) {
            yv = y;
            xv = x0;
        }
        else {
            yv = subset(y, idx);
            jsst = css(nv, yv.data());

            xv.assign(nv, 1);
            for (auto &e : ac) {
                auto z = subset(e, idx);
                xv.insert(xv.end(), z.begin(), z.end());
            }

            info = lmfit(nv, p0, yv.data(), xv.data(), &lm);

            if (info != 0) {
                eprint("ERROR: lmfit failed in assoc_glm (%s:%d): %d\n", __FILE__, __LINE__, info);
                return 2;
            }

            jdfe0 = lm.dfe;
            jsse0 = lm.sse;
        }

        isize_t p1 = p0 + length(gx);
        for (auto &e : gx)
            xv.insert(xv.end(), e.begin(), e.end());

        info = lmfit(nv, p1, yv.data(), xv.data(), &lm);

        if (info != 0) {
            eprint("ERROR: lmfit failed in assoc_glm (%s:%d): %d\n", __FILE__, __LINE__, info);
            return 3;
        }

        double jdfe1 = lm.dfe;
        double jsse1 = lm.sse;

        double dfx = jdfe0 - jdfe1;
        double ssx = jsse0 - jsse1;

        if (ic.empty()) {
            if (dfx > 0.0 && ssx > 0.0 && jdfe1 > 0.0 && jsse1 > 0.0) {
                auto f = (ssx / dfx) / (jsse1 / jdfe1);
                out[0][j] = fpval(f, dfx, jdfe1);
                out[1][j] = ssx / jsst;
            }
        }
        else {
            isize_t p2 = p1 + length(gx) * length(ic);
            for (auto &e : gx) {
                for (auto z : ic) {
                    if (nv != n)
                        z = subset(z, idx);
                    std::transform(e.begin(), e.end(), z.begin(), z.begin(), std::multiplies<double>());
                    xv.insert(xv.end(), z.begin(), z.end());
                }
            }

            info = lmfit(nv, p2, yv.data(), xv.data(), &lm);

            if (info != 0) {
                eprint("ERROR: lmfit failed in assoc_glm (%s:%d): %d\n", __FILE__, __LINE__, info);
                return 4;
            }

            double jdfe2 = lm.dfe;
            double jsse2 = lm.sse;

            double dfm = jdfe0 - jdfe2;
            double ssm = jsse0 - jsse2;

            if (dfm > 0.0 && ssm > 0.0 && jdfe2 > 0.0 && jsse2 > 0.0) {
                double f = (ssm / dfm) / (jsse2 / jdfe2);
                out[0][j] = fpval(f, dfm, jdfe2);
                out[1][j] = ssm / jsst;
            }

            if (dfx > 0.0 && ssx > 0.0 && jdfe2 > 0.0 && jsse2 > 0.0) {
                double f = (ssx / dfx) / (jsse2 / jdfe2);
                out[2][j] = fpval(f, dfx, jdfe2);
                out[3][j] = ssx / jsst;
            }

            double dfxi = jdfe1 - jdfe2;
            double ssxi = jsse1 - jsse2;

            if (dfxi > 0.0 && ssxi > 0.0 && jdfe2 > 0.0 && jsse2 > 0.0) {
                double f = (ssxi / dfxi) / (jsse2 / jdfe2);
                out[4][j] = fpval(f, dfxi, jdfe2);
                out[5][j] = ssxi / jsst;
            }
        }
    }

    return 0;
}

int assoc_glm_omp(const Genotype &gt, const std::vector<isize_t> &gi,
                  const std::vector<double> &y,
                  const std::vector< std::vector<double> > &ac,
                  const std::vector< std::vector<double> > &ic,
                  std::vector< std::vector<double> > &out)
{
    isize_t n = length(y);
    isize_t m = length(gt.dat);

    out.assign(ic.empty() ? 2 : 6, std::vector<double>(m, kNaN));

    double sst = css(n, y.data());

    isize_t p0 = 1 + length(ac);
    std::vector<double> x0(n, 1);
    for (auto &e : ac)
        x0.insert(x0.end(), e.begin(), e.end());

    LmFit lm0;
    int info = lmfit(n, p0, y.data(), x0.data(), &lm0);

    if (info != 0) {
        eprint("ERROR: lmfit failed in assoc_glm (%s:%d): %d\n", __FILE__, __LINE__, info);
        return 1;
    }

    double dfe0 = lm0.dfe;
    double sse0 = lm0.sse;

    #pragma omp parallel for
    for (isize_t j = 0; j < m; ++j) {
        std::vector<size_t> idx;
        std::vector< std::vector<double> > gx;

        if (gt.ploidy != 2) {
            std::vector<allele_t> g;
            for (isize_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][gi[i]];
                if (a) {
                    g.push_back(a);
                    idx.push_back(i);
                }
            }
            design2(grpidx(g), gx);
        }
        else {
            std::vector< std::pair<allele_t, allele_t> > g;
            for (isize_t i = 0; i < n; ++i) {
                isize_t i1 = gi[i] * 2;
                isize_t i2 = i1 + 1;
                auto a = gt.dat[j][i1];
                auto b = gt.dat[j][i2];
                if (a && b) {
                    if (a > b)
                        std::swap(a, b);
                    g.emplace_back(a, b);
                    idx.push_back(i);
                }
            }
            design2(grpidx(g), gx);
        }

        if (gx.empty())
            continue;

        isize_t nv = length(idx);

        double jsst = sst;
        double jdfe0 = dfe0;
        double jsse0 = sse0;

        std::vector<double> xv, yv, b;

        if (nv == n) {
            yv = y;
            xv = x0;
        }
        else {
            yv = subset(y, idx);
            jsst = css(nv, yv.data());

            xv.assign(nv, 1);
            for (auto &e : ac) {
                auto z = subset(e, idx);
                xv.insert(xv.end(), z.begin(), z.end());
            }

            LmFit lm;
            info = lmfit(nv, p0, yv.data(), xv.data(), &lm);

            if (info != 0) {
                eprint("ERROR: lmfit failed in assoc_glm (%s:%d): %d\n", __FILE__, __LINE__, info);
                continue;
            }

            jdfe0 = lm.dfe;
            jsse0 = lm.sse;
        }

        isize_t p1 = p0 + length(gx);
        for (auto &e : gx)
            xv.insert(xv.end(), e.begin(), e.end());

        LmFit lm;
        info = lmfit(nv, p1, yv.data(), xv.data(), &lm);

        if (info != 0) {
            eprint("ERROR: lmfit failed in assoc_glm (%s:%d): %d\n", __FILE__, __LINE__, info);
            continue;
        }

        double jdfe1 = lm.dfe;
        double jsse1 = lm.sse;

        double dfx = jdfe0 - jdfe1;
        double ssx = jsse0 - jsse1;

        if (ic.empty()) {
            if (dfx > 0.0 && ssx > 0.0 && jdfe1 > 0.0 && jsse1 > 0.0) {
                double f = (ssx / dfx) / (jsse1 / jdfe1);
                out[0][j] = fpval(f, dfx, jdfe1);
                out[1][j] = ssx / jsst;
            }
        }
        else {
            isize_t p2 = p1 + length(gx) * length(ic);
            for (auto &e : gx) {
                for (auto z : ic) {
                    if (nv != n)
                        z = subset(z, idx);
                    std::transform(e.begin(), e.end(), z.begin(), z.begin(), std::multiplies<double>());
                    xv.insert(xv.end(), z.begin(), z.end());
                }
            }

            LmFit lm2;
            info = lmfit(nv, p2, yv.data(), xv.data(), &lm2);

            if (info != 0) {
                eprint("ERROR: lmfit failed in assoc_glm (%s:%d): %d\n", __FILE__, __LINE__, info);
                continue;
            }

            double jdfe2 = lm2.dfe;
            double jsse2 = lm2.sse;

            double dfm = jdfe0 - jdfe2;
            double ssm = jsse0 - jsse2;

            if (dfm > 0 && ssm > 0 && jdfe2 > 0 && jsse2 > 0) {
                double f = (ssm / dfm) / (jsse2 / jdfe2);
                out[0][j] = fpval(f, dfm, jdfe2);
                out[1][j] = ssm / jsst;
            }

            if (dfx > 0 && ssx > 0 && jdfe2 > 0 && jsse2 > 0) {
                double f = (ssx / dfx) / (jsse2 / jdfe2);
                out[2][j] = fpval(f, dfx, jdfe2);
                out[3][j] = ssx / jsst;
            }

            double dfxi = jdfe1 - jdfe2;
            double ssxi = jsse1 - jsse2;

            if (dfxi > 0 && ssxi > 0 && jdfe2 > 0 && jsse2 > 0) {
                double f = (ssxi / dfxi) / (jsse2 / jdfe2);
                out[4][j] = fpval(f, dfxi, jdfe2);
                out[5][j] = ssxi / jsst;
            }
        }
    }

    return 0;
}

int glm_gwas(int argc, char *argv[])
{
    eprint("GLM-GWAS %s (Built on %s %s)\n", GLM_GWAS_VERSION, __DATE__, __TIME__);

    if (argc < 2) {
        show_help();
        return 1;
    }

    parse_cmdline(argc, argv);

    if (par.thread > 0)
        call_omp_set_num_threads(par.thread);

    Genotype gt;

    eprint("INFO: reading genotype file...\n");
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    eprint("INFO: %td individuals, %td loci\n", length(gt.ind), length(gt.loc));

    Phenotype pt;

    eprint("INFO: reading phenotype data file...\n");
    if (read_pheno(par.pheno, pt) != 0)
        return 1;
    eprint("INFO: %td observations, %td traits\n", length(pt.ind), length(pt.phe));

    Covariate ct;

    if (!par.covar.empty()) {
        eprint("INFO: reading covariate data file...\n");
        if (read_covar(par.covar, ct) != 0)
            return 1;
        eprint("INFO: %td individuals, %td covariates\n", length(ct.ind), length(ct.phe));
    }

    std::vector<isize_t> gi;
    intersect(gt, pt, ct, gi);

    if (gi.empty()) {
        eprint("ERROR: no valid observations are found\n");
        return 1;
    }

    std::vector< std::vector<double> > ac, ic;
    design(pt, ac, ic);

    ac.insert(ac.end(), ct.dat.begin(), ct.dat.end());

    isize_t nt = length(pt.phe);
    std::vector< std::vector<double> > out;

    for (isize_t t = 0; t < nt; ++t) {
        std::vector<char> keep;
        for (auto e : pt.dat[t])
            keep.push_back(std::isfinite(e));

        auto y = pt.dat[t];
        auto ac2 = ac;
        auto ic2 = ic;
        auto gi2 = gi;

        if (std::find(keep.begin(), keep.end(), 0) != keep.end()) {
            y = subset(y, keep);
            for (auto &e : ac2)
                e = subset(e, keep);
            for (auto &e : ic2)
                e = subset(e, keep);
            gi2 = subset(gi, keep);
        }

        std::vector< std::vector<double> > res;

        if (par.thread > 0)
            assoc_glm_omp(gt, gi2, y, ac2, ic2, res);
        else
            assoc_glm(gt, gi2, y, ac2, ic2, res);

        out.insert(out.end(), res.begin(), res.end());
    }

    CFile file(par.out, "w");

    if (!file) {
        eprint("ERROR: can't open file: %s\n", par.out);
        return 1;
    }

    fprint(file, "Locus\tChromosome\tPosition");
    for (auto &e : pt.phe) {
        fprint(file, "\t%s.model_p\t%s.model_r2", e, e);
        if (!ic.empty())
            fprint(file, "\t%s.main_p\t%s.main_r2\t%s.gxe_p\t%s.gxe_r2", e, e, e, e);
    }
    fprint(file, "\n");

    isize_t m = length(gt.loc);

    for (isize_t j = 0; j < m; ++j) {
        fprint(file, "%s\t%s\t%d", gt.loc[j], gt.chr[j], gt.pos[j]);
        for (auto &v : out)
            fprint(file, "\t%g", v[j]);
        fprint(file, "\n");
    }

    eprint("INFO: GLM-GWAS has finished successfully\n");

    return 0;
}
