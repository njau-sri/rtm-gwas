#include <cmath>
#include <fstream>
#include <iostream>
#include <functional>
#include "version.h"
#include "cmdline.h"
#include "util.h"
#include "vcf.h"
#include "pheno.h"
#include "statsutil.h"
#include "lsfit.h"
#include "stepreg.h"
#include "anova.h"


using std::ptrdiff_t;
using std::size_t;


namespace {


struct Parameter
{
    std::string vcf;
    std::string pheno;
    std::string covar;
    std::string out;
    double rsq = 0.95;
    double alpha = 0.01;
    double preselect = 0.05;
    int mtc = 0;
    bool nogxe = false;
    bool openmp = false;
} par ;


void filter_ind(const std::vector<size_t> &idx, Phenotype &pt)
{
    subset(pt.ind,idx).swap(pt.ind);

    if ( ! pt.env.empty() )
        subset(pt.env,idx).swap(pt.env);

    if ( ! pt.blk.empty() )
        subset(pt.blk,idx).swap(pt.blk);

    for (auto &e : pt.dat)
        subset(e,idx).swap(e);
}

void filter_ind(const std::vector<size_t> &idx, Covariate &ct)
{
    subset(ct.ind,idx).swap(ct.ind);

    for (auto &e : ct.dat)
        subset(e,idx).swap(e);
}

void merge(Genotype &gt, Phenotype &pt, Covariate &ct, std::vector<size_t> &gi)
{
    bool docovar = ! ct.phe.empty() && ! ct.ind.empty();

    bool domerge = gt.ind != pt.ind;
    if ( ! domerge && docovar )
        domerge = gt.ind != ct.ind;

    if ( ! domerge ) {
        gi.resize(gt.ind.size());
        std::iota(gi.begin(), gi.end(), 0);
        return;
    }

    std::cerr << "INFO: performing data intersection by individual...\n";

    auto ind = intersect(gt.ind, pt.ind);
    if ( docovar )
        ind = intersect(ind, ct.ind);

    gi.clear();
    std::vector<size_t> pi, ci;
    size_t n = pt.ind.size();

    for (size_t i = 0; i < n; ++i) {
        if ( std::binary_search(ind.begin(), ind.end(), pt.ind[i]) ) {
            pi.push_back(i);
            gi.push_back( index(gt.ind, pt.ind[i]) );
            if ( docovar )
                ci.push_back( index(ct.ind, pt.ind[i]) );
        }
    }

    filter_ind(pi, pt);

    if ( docovar )
        filter_ind(ci, ct);

    std::cerr << "INFO: there are " << ind.size() << " individuals after intersection\n";
}

void parse_envblk(const Phenotype &pt, std::vector< std::vector<double> > &ac, std::vector< std::vector<double> > &ic)
{
    std::vector< std::vector<double> > xenv, xblk;

    if ( ! pt.env.empty() ) {
        idummy1(factor(pt.env), xenv);
        ac.insert(ac.end(), xenv.begin(), xenv.end());
        ic.insert(ic.end(), xenv.begin(), xenv.end());
    }

    if ( ! pt.blk.empty() ) {
        idummy1(factor(pt.blk), xblk);
        if ( xenv.empty() )
            ac.insert(ac.end(), xblk.begin(), xblk.end());
        else {
            idummy3(factor(pt.env), xenv);
            for (auto &e : xenv) {
                for (auto v : xblk) {
                    std::transform(e.begin(), e.end(), v.begin(), v.begin(), std::multiplies<double>());
                    ac.push_back(v);
                }
            }
        }
    }
}

int assoc_glm(const Genotype &gt, const std::vector<size_t> &gi, const std::vector<double> &y,
              const std::vector< std::vector<double> > &ac,
              const std::vector< std::vector<double> > &ic,
              std::vector< std::vector<double> > &res)
{
    auto n = y.size();
    auto m = gt.dat.size();

    // modelP, modelR2, mainP, mainR2, intP, intR2
    res.assign(ic.empty() ? 2 : 6, std::vector<double>(m, std::numeric_limits<double>::quiet_NaN()));

    double sst = calc_css(y);

    std::vector<double> x0(n, 1);
    for (auto &v : ac)
        x0.insert(x0.end(), v.begin(), v.end());

    std::vector<double> b;
    double dfe0 = 0, sse0 = 0;
    lsfit(y, x0, b, dfe0, sse0);

    std::vector<allele_t> g1;
    std::vector< std::pair<allele_t,allele_t> > g2;

    std::vector<double> x, y2;

    for (size_t j = 0; j < m; ++j) {
        std::vector<size_t> idx;
        std::vector< std::vector<double> > x1;

        if (gt.ploidy == 1) {
            g1.clear();
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii];
                if ( a ) {
                    g1.push_back(a);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g1), x1);
        }
        else {
            g2.clear();
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii*2];
                auto b = gt.dat[j][ii*2+1];
                if ( a && b ) {
                    if (a > b)
                        std::swap(a, b);
                    g2.emplace_back(a, b);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g2), x1);
        }

        if ( x1.empty() )
            continue;

        auto n2 = idx.size();

        auto jsst = sst;
        auto jdfe0 = dfe0;
        auto jsse0 = sse0;

        if (n2 == n) {
            y2 = y;
            x = x0;
        }
        else {
            y2 = subset(y, idx);

            jsst = calc_css(y2);

            x.assign(n2, 1);
            for (auto &e : ac) {
                auto z = subset(e, idx);
                x.insert(x.end(), z.begin(), z.end());
            }

            lsfit(y2, x, b, jdfe0, jsse0);
        }

        for (auto &e : x1)
            x.insert(x.end(), e.begin(), e.end());

        double jdfe1 = 0, jsse1 = 0;
        lsfit(y2, x, b, jdfe1, jsse1);

        auto dfx = jdfe0 - jdfe1;
        auto ssx = jsse0 - jsse1;

        if ( ic.empty() ) {
            if (dfx > 0 && ssx > 0 && jdfe1 > 0 && jsse1 > 0) {
                auto f = (ssx / dfx) / (jsse1 / jdfe1);
                res[0][j] = fpval(f, dfx, jdfe1);
                res[1][j] = ssx / jsst;
            }
        }
        else {
            for (auto &e : x1) {
                for (auto z : ic) {
                    if (n2 != n)
                        z = subset(z, idx);
                    std::transform(e.begin(), e.end(), z.begin(), z.begin(), std::multiplies<double>());
                    x.insert(x.end(), z.begin(), z.end());
                }
            }

            double jdfe2 = 0, jsse2 = 0;
            lsfit(y2, x, b, jdfe2, jsse2);

            auto dfm = jdfe0 - jdfe2;
            auto ssm = jsse0 - jsse2;

            if (dfm > 0 && ssm > 0 && jdfe2 > 0 && jsse2 > 0) {
                auto f = (ssm / dfm) / (jsse2 / jdfe2);
                res[0][j] = fpval(f, dfm, jdfe2);
                res[1][j] = ssm / jsst;
            }

            if (dfx > 0 && ssx > 0 && jdfe2 > 0 && jsse2 > 0) {
                auto f = (ssx / dfx) / (jsse2 / jdfe2);
                res[2][j] = fpval(f, dfx, jdfe2);
                res[3][j] = ssx / jsst;
            }

            auto dfxi = jdfe1 - jdfe2;
            auto ssxi = jsse1 - jsse2;

            if (dfxi > 0 && ssxi > 0 && jdfe2 > 0 && jsse2 > 0) {
                auto f = (ssxi / dfxi) / (jsse2 / jdfe2);
                res[4][j] = fpval(f, dfxi, jdfe2);
                res[5][j] = ssxi / jsst;
            }
        }
    }

    return 0;
}

int assoc_glm_omp(const Genotype &gt, const std::vector<size_t> &gi, const std::vector<double> &y,
                  const std::vector< std::vector<double> > &ac,
                  const std::vector< std::vector<double> > &ic,
                  std::vector< std::vector<double> > &res)
{
    auto n = y.size();
    auto m = gt.dat.size();

    // modelP, modelR2, mainP, mainR2, intP, intR2
    res.assign(ic.empty() ? 2 : 6, std::vector<double>(m, std::numeric_limits<double>::quiet_NaN()));

    double sst = calc_css(y);

    std::vector<double> x0(n, 1);
    for (auto &v : ac)
        x0.insert(x0.end(), v.begin(), v.end());

    std::vector<double> b0;
    double dfe0 = 0, sse0 = 0;
    lsfit(y, x0, b0, dfe0, sse0);

    // in earlier OpenMP specifications (<3.0), unsigned integer is not allowed in loop construct
    auto m2 = static_cast<ptrdiff_t>(m);

    #pragma omp parallel for
    for (ptrdiff_t j2 = 0; j2 < m2; ++j2) {
        auto j = static_cast<size_t>(j2);

        std::vector<size_t> idx;
        std::vector< std::vector<double> > x1;

        if (gt.ploidy == 1) {
            std::vector<allele_t> g;
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii];
                if ( a ) {
                    g.push_back(a);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g), x1);
        }
        else {
            std::vector< std::pair<allele_t,allele_t> > g;
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii*2];
                auto b = gt.dat[j][ii*2+1];
                if ( a && b ) {
                    if (a > b)
                        std::swap(a, b);
                    g.emplace_back(a, b);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g), x1);
        }

        if ( x1.empty() )
            continue;

        auto n2 = idx.size();

        auto jsst = sst;
        auto jdfe0 = dfe0;
        auto jsse0 = sse0;

        std::vector<double> x, y2, b;

        if (n2 == n) {
            y2 = y;
            x = x0;
        }
        else {
            y2 = subset(y, idx);

            jsst = calc_css(y2);

            x.assign(n2, 1);
            for (auto &e : ac) {
                auto z = subset(e, idx);
                x.insert(x.end(), z.begin(), z.end());
            }

            lsfit(y2, x, b, jdfe0, jsse0);
        }

        for (auto &e : x1)
            x.insert(x.end(), e.begin(), e.end());

        double jdfe1 = 0, jsse1 = 0;
        lsfit(y2, x, b, jdfe1, jsse1);

        auto dfx = jdfe0 - jdfe1;
        auto ssx = jsse0 - jsse1;

        if ( ic.empty() ) {
            if (dfx > 0 && ssx > 0 && jdfe1 > 0 && jsse1 > 0) {
                auto f = (ssx / dfx) / (jsse1 / jdfe1);
                res[0][j] = fpval(f, dfx, jdfe1);
                res[1][j] = ssx / jsst;
            }
        }
        else {
            for (auto &e : x1) {
                for (auto z : ic) {
                    if (n2 != n)
                        z = subset(z, idx);
                    std::transform(e.begin(), e.end(), z.begin(), z.begin(), std::multiplies<double>());
                    x.insert(x.end(), z.begin(), z.end());
                }
            }

            double jdfe2 = 0, jsse2 = 0;
            lsfit(y2, x, b, jdfe2, jsse2);

            auto dfm = jdfe0 - jdfe2;
            auto ssm = jsse0 - jsse2;

            if (dfm > 0 && ssm > 0 && jdfe2 > 0 && jsse2 > 0) {
                auto f = (ssm / dfm) / (jsse2 / jdfe2);
                res[0][j] = fpval(f, dfm, jdfe2);
                res[1][j] = ssm / jsst;
            }

            if (dfx > 0 && ssx > 0 && jdfe2 > 0 && jsse2 > 0) {
                auto f = (ssx / dfx) / (jsse2 / jdfe2);
                res[2][j] = fpval(f, dfx, jdfe2);
                res[3][j] = ssx / jsst;
            }

            auto dfxi = jdfe1 - jdfe2;
            auto ssxi = jsse1 - jsse2;

            if (dfxi > 0 && ssxi > 0 && jdfe2 > 0 && jsse2 > 0) {
                auto f = (ssxi / dfxi) / (jsse2 / jdfe2);
                res[4][j] = fpval(f, dfxi, jdfe2);
                res[5][j] = ssxi / jsst;
            }
        }
    }

    return 0;
}

int assoc_stepwise(const Genotype &gt, const std::vector<size_t> &gi, const std::vector<double> &y,
                   const std::vector< std::vector<double> > &ac,
                   const std::vector< std::vector<double> > &ic,
                   std::vector<double> &ps, std::vector<size_t> &in)
{
    auto n = y.size();
    auto m = gt.dat.size();

    std::vector<double> x0(n,1);
    for (auto &e : ac)
        x0.insert(x0.end(), e.begin(), e.end());

    std::vector<size_t> idx;
    std::vector<allele_t> g1;
    std::vector<std::pair<allele_t, allele_t> > g2;

    std::vector<size_t> c1;
    std::vector<double> x1;

    for (size_t j = 0; j < m; ++j) {
        std::vector< std::vector<double> > xg;

        // !!! missing genotype is not allowed !!!

        if (gt.ploidy == 1) {
            g1.clear();
            for (auto i : gi)
                g1.push_back(gt.dat[j][i]);
            idummy1(factor(g1), xg);
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
            idummy1(factor(g2), xg);
        }

        if (xg.empty())
            continue;

        idx.push_back(j);

        if ( ! ic.empty() ) {
            auto ng = xg.size();
            for (size_t k = 0; k < ng; ++k) {
                for (auto v : ic) {
                    std::transform(v.begin(), v.end(), xg[k].begin(), v.begin(), std::multiplies<double>());
                    xg.push_back(v);
                }
            }
        }

        c1.push_back( xg.size() );

        for (auto &e : xg)
            x1.insert(x1.end(), e.begin(), e.end());
    }

    StepReg sr;
    sr.sle = sr.sls = par.alpha;
    sr.rsq = par.rsq;
    sr.mtc = par.mtc;

    std::vector<double> px;

    if ( par.openmp )
        sr.forward_omp(x0, x1, c1, y, px, in);
    else
        sr.forward(x0, x1, c1, y, px, in);

    sr.backward(x0, x1, c1, y, px, in);

    subset(idx,in).swap(in);

    ps.assign(m, std::numeric_limits<double>::quiet_NaN());

    auto nx = idx.size();
    for (size_t i = 0; i < nx; ++i)
        ps[idx[i]] = px[i];

    return 0;
}

int fit_multi_locus_model(size_t t, const std::vector<size_t> &loc,
                          const Genotype &gt, const std::vector<size_t> &gi,
                          const Phenotype &pt, const Covariate &ct,
                          std::string &est, std::string &aov1, std::string &aov3)
{
    auto y = pt.dat[t];
    auto n = y.size();

    std::vector<size_t> obs;
    for (size_t i = 0; i < n; ++i) {
        if (std::isfinite(y[i]))
            obs.push_back(i);
    }

    auto nv = obs.size();
    if (nv != n)
        subset(y,obs).swap(y);

    std::vector<std::string> env;
    if ( ! pt.env.empty() )
        subset(pt.env,obs).swap(env);

    std::vector<std::string> blk;
    if ( ! pt.blk.empty() )
        subset(pt.blk,obs).swap(blk);

    ANOVA aov;

    if ( ! env.empty() )
        aov.add_main("_ENV_", env);

    if ( ! blk.empty() ) {
        if ( env.empty() )
            aov.add_main("_BLK_", blk);
        else
            aov.add_nested("_BLK_(_ENV_)", blk, env);
    }

    if ( ! ct.phe.empty() && ! ct.ind.empty() ) {
        auto m = ct.phe.size();
        for (size_t i = 0; i < m; ++i) {
            if (nv == n)
                aov.add_reg(ct.phe[i], ct.dat[i]);
            else
                aov.add_reg(ct.phe[i], subset(ct.dat[i],obs));
        }
    }

    subset(gi,obs).swap(obs);

    for (auto j : loc) {
        std::vector<std::string> gs;
        if (gt.ploidy == 1) {
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

        aov.add_main(gt.loc[j], gs);

        if ( ! env.empty() && ! par.nogxe )
            aov.add_crossed(gt.loc[j] + "*_ENV_", gs, env);
    }

    auto sol = aov.solution(y);
    auto at1 = aov.solve1(y);
    auto at3 = aov.solve3(y);

    for (auto &e : at1.p) {
        if (e == 0.0)
            e = std::numeric_limits<double>::min();
    }

    for (auto &e : at3.p) {
        if (e == 0.0)
            e = std::numeric_limits<double>::min();
    }

    est = sol.to_string();
    aov1 = at1.to_string();
    aov3 = at3.to_string();

    return 0;
}

} // namespace


int rtm_gwas_assoc(int argc, char *argv[])
{
    std::cerr << "RTM-GWAS " RTM_GWAS_VERSION_STRING " ASSOC (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd;

    cmd.add("--vcf", "VCF file", "");
    cmd.add("--pheno", "phenotype file", "");
    cmd.add("--covar", "covariate file", "");
    cmd.add("--out", "output file", "assoc.out");
    cmd.add("--alpha", "significance level", "0.01");
    cmd.add("--preselect", "pre-selection threshold", "0.05");
    cmd.add("--mtc", "multiple testing correction, BON/FDR", "");
    cmd.add("--rsq", "maximum model r-square", "0.95");

    cmd.add("--no-gxe", "ignore GxE interaction effect");
    cmd.add("--openmp", "enable OpenMP multithreading");

    cmd.parse(argc, argv);

    if (argc < 2) {
        cmd.show();
        return 1;
    }

    par.vcf = cmd.get("--vcf");
    par.pheno = cmd.get("--pheno");
    par.covar = cmd.get("--covar");
    par.out = cmd.get("--out");
    par.alpha = std::stod(cmd.get("--alpha"));
    par.preselect = std::stod(cmd.get("--preselect"));
    par.rsq = std::stod(cmd.get("--rsq"));

    par.nogxe = cmd.has("--no-gxe");
    par.openmp = cmd.has("--openmp");

    if ( ! cmd.get("--mtc").empty() ) {
        auto mtc = cmd.get("--mtc");
        std::transform(mtc.begin(), mtc.end(), mtc.begin(), ::toupper);
        if (mtc == "BON")
            par.mtc = 1;
        else if (mtc == "FDR")
            par.mtc = 2;
        else if (mtc == "HOLM")
            par.mtc = 3;
        else
            par.mtc = 1;
    }

    Genotype gt;
    Phenotype pt;
    Covariate ct;

    std::cerr << "INFO: reading genotype file...\n";
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    std::cerr << "INFO: reading phenotype file...\n";
    if (read_pheno(par.pheno, pt) != 0)
        return 1;
    std::cerr << "INFO: " << pt.ind.size() << " observations, " << pt.phe.size() << " traits\n";

    if ( ! par.covar.empty() ) {
        std::cerr << "INFO: reading covariate file...\n";
        if (read_covar(par.covar, ct) != 0)
            return 1;
        std::cerr << "INFO: " << ct.ind.size() << " individuals, " << ct.phe.size() << " covariates\n";
    }

    std::vector<size_t> gi;
    merge(gt, pt, ct, gi);

    if ( gi.empty() ) {
        std::cerr << "ERROR: no valid observations are found\n";
        return 1;
    }

    std::vector< std::vector<double> > ac, ic;
    parse_envblk(pt, ac, ic);

    if (par.nogxe)
        ic.clear();

    ac.insert(ac.end(), ct.dat.begin(), ct.dat.end());

    std::ofstream ofs_loc(par.out + ".loc");
    if ( ! ofs_loc ) {
        std::cerr << "ERROR: can't open file: " << par.out << ".loc" << "\n";
        return 1;
    }

    std::ofstream ofs_est(par.out + ".est");
    if (!ofs_est) {
        std::cerr << "ERROR: can't open file: " << par.out << ".est" << "\n";
        return 1;
    }

    std::ofstream ofs_aov1(par.out + ".aov1");
    if (!ofs_aov1) {
        std::cerr << "ERROR: can't open file: " << par.out << ".aov1" << "\n";
        return 1;
    }

    std::ofstream ofs_aov3(par.out + ".aov3");
    if (!ofs_aov3) {
        std::cerr << "ERROR: can't open file: " << par.out << ".aov3" << "\n";
        return 1;
    }

    auto m = gt.loc.size();
    auto nt = pt.phe.size();
    std::vector< std::vector<double> > ps;

    for (size_t t = 0; t < nt; ++t) {
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

        std::vector<size_t> loc;
        std::vector<double> p1;

        if (par.preselect <= 0 || par.preselect >= 1)
            assoc_stepwise(gt, gi2, y, ac2, ic2, p1, loc);
        else {
            std::vector< std::vector<double> > glm;
            if ( par.openmp )
                assoc_glm_omp(gt, gi2, y, ac2, ic2, glm);
            else
                assoc_glm(gt, gi2, y, ac2, ic2, glm);

            std::vector<size_t> idx;

            if ( ic.empty() ) {
                p1 = glm[0];
                for (size_t j = 0; j < m; ++j) {
                    auto pval = glm[0][j];
                    if (std::isfinite(pval) && pval <= par.preselect)
                        idx.push_back(j);
                }
            }
            else {
                p1 = glm[2];
                for (size_t j = 0; j < m; ++j) {
                    auto pval = glm[2][j];
                    if (std::isfinite(pval) && pval <= par.preselect) {
                        idx.push_back(j);
                        continue;
                    }
                    pval = glm[4][j];
                    if (std::isfinite(pval) && pval <= par.preselect)
                        idx.push_back(j);
                }
            }

            std::cerr << "INFO: after pre-selection, there are " << idx.size() << " loci\n";

            if ( ! idx.empty() ) {
                Genotype gt2;
                gt2.ploidy = gt.ploidy;
                subset(gt.dat,idx).swap(gt2.dat);
                subset(gt.loc,idx).swap(gt2.loc);

                std::vector<double> p2;
                assoc_stepwise(gt2, gi2, y, ac2, ic2, p2, loc);

                subset(idx,loc).swap(loc);

                for (size_t k = 0; k < idx.size(); ++k)
                    p1[idx[k]] = p2[k];
            }
        }

        ps.push_back(p1);

        ofs_loc << ">" << pt.phe[t] << "\n";
        for (auto j : loc)
            ofs_loc << gt.loc[j] << "\n";
        ofs_loc << "\n";

        std::string est, aov1, aov3;
        fit_multi_locus_model(t, loc, gt, gi, pt, ct, est, aov1, aov3);

        ofs_est << ">" << pt.phe[t] << "\n" << est << "\n";
        ofs_aov1 << ">" << pt.phe[t] << "\n" << aov1 << "\n";
        ofs_aov3 << ">" << pt.phe[t] << "\n" << aov3 << "\n";
    }

    std::ofstream ofs_ps(par.out + ".pval");
    if ( ! ofs_ps ) {
        std::cerr << "ERROR: can't open file: " << par.out << ".pval" << "\n";
        return 1;
    }

    ofs_ps << "Locus\tChromosome\tPosition";
    for (auto &e : pt.phe)
        ofs_ps << "\t" << e;
    ofs_ps << "\n";

    for (size_t j = 0; j < m; ++j) {
        ofs_ps << gt.loc[j] << "\t" << gt.chr[j] << "\t" << gt.pos[j];
        for (size_t t = 0; t < nt; ++t)
            ofs_ps << "\t" << ps[t][j];
        ofs_ps << "\n";
    }

    return  0;
}
