#include "rtmgwasassoc.h"

#include <cmath>
#include <limits>
#include <functional>
#include "cmdline.h"
#include "print.h"
#include "stringutil.h"
#include "vectorutil.h"
#include "statutil.h"
#include "vcf.h"
#include "pheno.h"
#include "anova.h"
#include "openmp.h"
#include "cfile.h"
#include "glmgwas.h"
#include "lmfit.h"

#ifndef RTM_GWAS_VERSION
#define RTM_GWAS_VERSION "unknown"
#endif // RTM_GWAS_VERSION

namespace {

    double kNaN = std::numeric_limits<double>::quiet_NaN();

    struct Parameter
    {
        std::string vcf;
        std::string pheno;
        std::string covar;
        std::string out;
        double rsq = 0.95;
        double alpha = 0.01;
        double preselect = 0.05;
        double dfe = 0.0;
        double mse = 0.0;
        int mtc = 0;
        int thread = 0;
        bool nogxe = false;
    } par;

    void show_help()
    {
        const char *help =
            "Usage: rtm-gwas-assoc [options...]\n"
            "Options:\n"
            "  --vcf <>            genotype file in VCF format\n"
            "  --pheno <>          phenotype data file\n"
            "  --covar <>          covariate data file\n"
            "  --out <assoc.out>   output file prefix\n"
            "  --alpha <0.01>      significance level\n"
            "  --preselect <0.05>  pre-selection threshold\n"
            "  --rsq <0.95>        maximum model r-square\n"
            //"  --dfe <0>           degrees of freedom of known error\n"
            //"  --mse <0>           mean squares of known error\n"
            "  --mtc <>            multiple testing correction, BON/FDR\n"
            "  --thread <0>        set the number of threads\n"
            "  --no-gxe            exclude GxE interaction effect\n"
            "\n";
        eprint(help);
    }

    int parse_cmdline(int argc, char *argv[])
    {
        CmdLine cmd;

        cmd.add("--vcf", "");
        cmd.add("--pheno", "");
        cmd.add("--covar", "");
        cmd.add("--out", "assoc.out");
        cmd.add("--alpha", "0.01");
        cmd.add("--preselect", "0.05");
        cmd.add("--rsq", "0.95");
        cmd.add("--dfe", "0");
        cmd.add("--mse", "0");
        cmd.add("--mtc", "");
        cmd.add("--thread", "0");
        cmd.add("--no-gxe");

        if (cmd.parse(argc, argv) != 0)
            return 1;

        par.vcf = cmd.get("--vcf");
        par.pheno = cmd.get("--pheno");
        par.covar = cmd.get("--covar");
        par.out = cmd.get("--out");
        par.alpha = std::stod(cmd.get("--alpha"));
        par.preselect = std::stod(cmd.get("--preselect"));
        par.rsq = std::stod(cmd.get("--rsq"));
        par.dfe = std::stod(cmd.get("--dfe"));
        par.mse = std::stod(cmd.get("--mse"));
        par.thread = std::stoi(cmd.get("--thread"));
        par.nogxe = cmd.has("--no-gxe");

        std::string mtc = to_upper(cmd.get("--mtc"));
        if (!mtc.empty()) {
            if (mtc == "BON")
                par.mtc = 1;
            else if (mtc == "FDR")
                par.mtc = 2;
            else if (mtc == "HOLM")
                par.mtc = 3;
            else
                eprint("WARNING: invalid multiple testing correction: %s\n", mtc);
        }

        return 0;
    }

    void fit_multi_locus_model(isize_t t, const std::vector<isize_t> &loc,
                               const Genotype &gt, const std::vector<isize_t> &gi,
                               const Phenotype &pt, const Covariate &ct,
                               std::string &est, std::string &estwt, std::string &aov1, std::string &aov3)
    {
        auto y = pt.dat[t];
        isize_t n = length(y);

        std::vector<isize_t> obs;
        for (isize_t i = 0; i < n; ++i) {
            if (std::isfinite(y[i]))
                obs.push_back(i);
        }

        isize_t nv = length(obs);
        if (nv != n)
            subset(y, obs).swap(y);

        std::vector<std::string> env;
        if (!pt.env.empty())
            subset(pt.env, obs).swap(env);

        std::vector<std::string> blk;
        if (!pt.blk.empty())
            subset(pt.blk, obs).swap(blk);

        ANOVA anova;

        if (!env.empty())
            anova.add_main("_ENV_", env);

        if (!blk.empty()) {
            if (env.empty())
                anova.add_main("_BLK_", blk);
            else
                anova.add_nested("_BLK_(_ENV_)", blk, env);
        }

        if (!ct.phe.empty() && !ct.ind.empty()) {
            isize_t nc = length(ct.phe);
            for (isize_t i = 0; i < nc; ++i) {
                if (nv == n)
                    anova.add_reg(ct.phe[i], ct.dat[i]);
                else
                    anova.add_reg(ct.phe[i], subset(ct.dat[i], obs));
            }
        }

        subset(gi, obs).swap(obs);

        for (auto j : loc) {
            std::vector<std::string> gs;

            // NOTE: assuming complete genotype data

            if (gt.ploidy != 2) {
                for (auto i : obs) {
                    auto a = gt.dat[j][i];
                    if (a)
                        gs.push_back(gt.allele[j][a-1]);
                    else
                        gs.push_back(".");
                }
            }
            else {
                for (auto i : obs) {
                    auto a = gt.dat[j][i*2];
                    auto b = gt.dat[j][i*2+1];
                    if (a > b)
                        std::swap(a, b);
                    if (a && b)
                        gs.push_back(gt.allele[j][a-1] + '/' + gt.allele[j][b-1]);
                    else
                        gs.push_back("./.");
                }
            }

            anova.add_main(gt.loc[j], gs);

            if (!env.empty() && !par.nogxe)
                anova.add_crossed(gt.loc[j] + "*_ENV_", gs, env);
        }

        est = anova.solution(y).to_string();
        estwt = anova.solution_wtsum(y).to_string();

        auto aov = anova.solve1(y);
        for (auto &e : aov.p) {
            if (e == 0.0)
                e = std::numeric_limits<double>::min();
        }
        aov1 = aov.to_string();

        aov = anova.solve3(y);
        for (auto &e : aov.p) {
            if (e == 0.0)
                e = std::numeric_limits<double>::min();
        }
        aov3 = aov.to_string();
    }

} // namespace

int assoc_stepwise(const Genotype &gt, const std::vector<isize_t> &gi,
                   const std::vector<double> &y,
                   const std::vector< std::vector<double> > &ac,
                   const std::vector< std::vector<double> > &ic,
                   std::vector<double> &ps, std::vector<isize_t> &in)
{
    isize_t n = length(y);
    isize_t m = length(gt.dat);

    isize_t q = 1 + length(ac);
    std::vector<double> w(n, 1);
    for (auto &e : ac)
        w.insert(w.end(), e.begin(), e.end());

    std::vector<double> x;
    std::vector<isize_t> c;
    std::vector<isize_t> idx;

    std::vector<allele_t> g1;
    std::vector<std::pair<allele_t, allele_t> > g2;

    for (isize_t j = 0; j < m; ++j) {
        std::vector< std::vector<double> > xg;

        // NOTE: assuming complete genotype data

        if (gt.ploidy != 2) {
            g1.clear();
            for (auto i : gi)
                g1.push_back(gt.dat[j][i]);
            design2(grpidx(g1), xg);
        }
        else {
            g2.clear();
            for (auto i : gi) {
                auto a = gt.dat[j][i*2];
                auto b = gt.dat[j][i*2+1];
                if (a > b)
                    std::swap(a, b);
                g2.emplace_back(a, b);
            }
            design2(grpidx(g2), xg);
        }

        if (xg.empty())
            continue;

        idx.push_back(j);

        if (!ic.empty()) {
            isize_t ng = length(xg);
            for (isize_t k = 0; k < ng; ++k) {
                for (auto v : ic) {
                    std::transform(v.begin(), v.end(), xg[k].begin(), v.begin(), std::multiplies<double>());
                    xg.push_back(v);
                }
            }
        }

        c.push_back(length(xg));

        for (auto &e : xg)
            x.insert(x.end(), e.begin(), e.end());
    }

    isize_t p = length(idx);

    StepwiseFit out;
    out.par.sle = out.par.sls = par.alpha;
    out.par.dfe = par.dfe;
    out.par.mse = par.mse;
    out.par.rsq = par.rsq;
    out.par.thread = par.thread;
    out.par.mtc = par.mtc;

    stepwisefit(n, q, p, y.data(), w.data(), x.data(), c.data(), &out);

    subset(idx, out.in).swap(in);

    ps.assign(m, kNaN);
    for (isize_t i = 0; i < p; ++i)
        ps[idx[i]] = out.ps[i];

    return 0;
}

int rtm_gwas_assoc(int argc, char *argv[])
{
    eprint("RTM-GWAS %s ASSOC (Built on %s %s)\n", RTM_GWAS_VERSION, __DATE__, __TIME__);

    if (argc < 2) {
        show_help();
        return 1;
    }

    if (parse_cmdline(argc, argv) != 0)
        return 1;

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

    if (par.nogxe)
        ic.clear();

    ac.insert(ac.end(), ct.dat.begin(), ct.dat.end());

    CFile locfile(par.out + ".loc", "w");
    if (!locfile) {
        eprint("ERROR: can't open file: %s.loc\n", par.out);
        return 1;
    }

    CFile estfile(par.out + ".est", "w");
    if (!estfile) {
        eprint("ERROR: can't open file: %s.est\n", par.out);
        return 1;
    }

    CFile estwtfile(par.out + ".est.wt", "w");
    if (!estwtfile) {
        eprint("ERROR: can't open file: %s.est.wt\n", par.out);
        return 1;
    }

    CFile aov1file(par.out + ".aov1", "w");
    if (!aov1file) {
        eprint("ERROR: can't open file: %s.aov1\n", par.out);
        return 1;
    }

    CFile aov3file(par.out + ".aov3", "w");
    if (!aov3file) {
        eprint("ERROR: can't open file: %s.aov3\n", par.out);
        return 1;
    }

    isize_t m = length(gt.loc);
    isize_t nt = length(pt.phe);
    std::vector< std::vector<double> > ps;

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

        std::vector<isize_t> loc;
        std::vector<double> p1;

        if (par.preselect <= 0.0 || par.preselect >= 1.0)
            assoc_stepwise(gt, gi2, y, ac2, ic2, p1, loc);
        else {
            std::vector< std::vector<double> > glm;
            if (par.thread > 0)
                assoc_glm_omp(gt, gi2, y, ac2, ic2, glm);
            else
                assoc_glm(gt, gi2, y, ac2, ic2, glm);

            std::vector<isize_t> idx;

            if (ic.empty()) {
                p1 = glm[0];
                for (isize_t j = 0; j < m; ++j) {
                    double pval = glm[0][j];
                    if (std::isfinite(pval) && pval <= par.preselect)
                        idx.push_back(j);
                }
            }
            else {
                p1 = glm[2];
                for (isize_t j = 0; j < m; ++j) {
                    double pval = glm[2][j];
                    if (std::isfinite(pval) && pval <= par.preselect) {
                        idx.push_back(j);
                        continue;
                    }
                    pval = glm[4][j];
                    if (std::isfinite(pval) && pval <= par.preselect)
                        idx.push_back(j);
                }
            }

            eprint("INFO: after pre-selection, there are %td loci left\n", length(idx));

            if (!idx.empty()) {
                Genotype gtsub;
                gtsub.ploidy = gt.ploidy;
                subset(gt.dat, idx).swap(gtsub.dat);
                subset(gt.loc, idx).swap(gtsub.loc);

                std::vector<double> p2;
                assoc_stepwise(gtsub, gi2, y, ac2, ic2, p2, loc);

                subset(idx, loc).swap(loc);

                for (isize_t k = 0; k < length(idx); ++k)
                    p1[idx[k]] = p2[k];
            }
        }

        ps.push_back(p1);

        fprint(locfile, ">%s\n", pt.phe[t]);
        for (auto j : loc)
            fprint(locfile, "%s\n", gt.loc[j]);
        fprint(locfile, "\n");

        std::string est, estwt, aov1, aov3;
        fit_multi_locus_model(t, loc, gt, gi, pt, ct, est, estwt, aov1, aov3);

        fprint(estfile, ">%s\n%s\n", pt.phe[t], est);
        fprint(estwtfile, ">%s\n%s\n", pt.phe[t], estwt);
        fprint(aov1file, ">%s\n%s\n", pt.phe[t], aov1);
        fprint(aov3file, ">%s\n%s\n", pt.phe[t], aov3);
    }

    CFile pvalfile(par.out + ".pval", "w");

    if (!pvalfile) {
        eprint("ERROR: can't open file: %s.pval\n", par.out);
        return 1;
    }

    fprint(pvalfile, "Locus\tChromosome\tPosition");
    for (auto &e : pt.phe)
        fprint(pvalfile, "\t%s", e);
    fprint(pvalfile, "\n");

    for (isize_t j = 0; j < m; ++j) {
        fprint(pvalfile, "%s\t%s\t%d", gt.loc[j], gt.chr[j], gt.pos[j]);
        for (isize_t t = 0; t < nt; ++t)
            fprint(pvalfile, "\t%g", ps[t][j]);
        fprint(pvalfile, "\n");
    }

    return  0;
}
