#include <cmath>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include "appassoc.h"
#include "cmdline.h"
#include "util.h"
#include "stat.h"
#include "plinkio.h"
#include "hapmapio.h"
#include "vcfio.h"
#include "rtmio.h"
#include "glm.h"
#include "assoclm.h"
#include "assoclmm.h"
#include "assocrtm.h"

namespace {

void filter_ind(const vector<size_t> &idx, Genotype &gt)
{
    bool haploid = gt.ploidy == 1;

    if ( haploid ) {
        for (auto &v : gt.dat)
            subset(v,idx).swap(v);
    }
    else {
        vector<size_t> idx2;
        for (auto i : idx) {
            idx2.push_back(i*2);
            idx2.push_back(i*2+1);
        }
        for (auto &v : gt.dat)
            subset(v,idx2).swap(v);
    }

    subset(gt.ind,idx).swap(gt.ind);
}

void filter_ind(const vector<size_t> &idx, Phenotype &pt)
{
    subset(pt.ind,idx).swap(pt.ind);

    if ( ! pt.env.empty() )
        subset(pt.env,idx).swap(pt.env);

    if ( ! pt.blk.empty() )
        subset(pt.blk,idx).swap(pt.blk);

    for (auto &v : pt.dat)
        subset(v,idx).swap(v);
}

void filter_ind(const vector<size_t> &idx, SquareData &sd)
{
    subset(sd.ind,idx).swap(sd.ind);
    subset(sd.dat,idx).swap(sd.dat);
    for (auto &v : sd.dat)
        subset(v,idx).swap(v);
}

} // namespace

int AppAssoc::run(int argc, char *argv[])
{
    auto cmd = std::make_shared<CmdLine>("assoc [options]");

    cmd->add("--vcf", "VCF file", "");
    cmd->add("--ped", "PLINK PED/MAP file prefix", "");
    cmd->add("--hmp", "HapMap file",  "");
    cmd->add("--geno", "genotype data file", "");
    cmd->add("--pheno", "phenotype data file", "");
    cmd->add("--covar", "covariate data file", "");
    cmd->add("--kin", "kinship data file", "");
    cmd->add("--out", "output file", "appassoc.out");
    cmd->add("--method", "association analysis method, RTM/LM/LMM", "RTM");
    cmd->add("--alpha", "RTM significance level", "0.05");
    cmd->add("--preselect", "RTM pre-selection threshold", "0.05");
    cmd->add("--multtest", "RTM multiple testing correction, BON/FDR",  "");
    cmd->add("--maxrsq", "RTM maximum model R-Square", "0.99");
    cmd->add("--maxqtl", "RTM maximum number of QTLs", "0");
    cmd->add("--sstype", "RTM sum of squares type", "1");

    cmd->add("--no-gxe", "ignore GxE interaction effect");

    if (argc < 2) {
        cmd->help();
        return 1;
    }

    cmd->parse(argc, argv);

    m_par.vcf = cmd->get("--vcf");
    m_par.ped = cmd->get("--ped");
    m_par.hmp = cmd->get("--hmp");
    m_par.geno = cmd->get("--geno");
    m_par.pheno = cmd->get("--pheno");
    m_par.covar = cmd->get("--covar");
    m_par.kin = cmd->get("--kin");
    m_par.out = cmd->get("--out");
    m_par.method = cmd->get("--method");
    m_par.multtest = cmd->get("--multtest");
    m_par.alpha = number<double>(cmd->get("--alpha"));
    m_par.preselect = number<double>(cmd->get("--preselect"));
    m_par.maxrsq = number<double>(cmd->get("--maxrsq"));
    m_par.maxqtl = number<int>(cmd->get("--maxqtl"));
    m_par.sstype = number<int>(cmd->get("--sstype"));

    m_par.no_gxe = cmd->has("--no-gxe");

    cmd.reset();

    std::transform(m_par.method.begin(), m_par.method.end(), m_par.method.begin(), ::toupper);
    std::transform(m_par.multtest.begin(), m_par.multtest.end(), m_par.multtest.begin(), ::toupper);

    if (m_par.method != "LMM" && !m_par.kin.empty()) {
        std::cerr << "WARNING: unused argument --kin " << m_par.kin << "\n";
        m_par.kin.clear();
    }

    int info = perform();

    return info;
}

int AppAssoc::perform()
{
    load_genotype();

    load_phenotype();

    load_covariate();

    load_kinship();

    merge_data();

    if ( m_obs.empty() ) {
        std::cerr << "ERROR: no valid observations were found\n";
        return 1;
    }

    parse_factor();

    if (m_par.method == "LM") {
        LM();
    }
    else if (m_par.method == "LMM") {
        if ( m_kin.ind.empty() ) {
            std::cerr << "ERROR: kinship data is requied for LMM method.\n";
            return 1;
        }
        LMM();
    }
    else if (m_par.method == "RTM") {
        RTM();
    }
    else {
        std::cerr << "ERROR: invalid analysis method: " << m_par.method << "\n";
        return 1;
    }

    return 0;
}

void AppAssoc::load_genotype()
{
    if ( m_par.vcf.empty() && m_par.ped.empty() && m_par.hmp.empty() && m_par.geno.empty() )
        return;

    std::cerr << "INFO: reading genotype file...\n";

    int info = 0;

    if ( ! m_par.vcf.empty() )
        info = read_vcf(m_par.vcf, m_gt);
    else if ( ! m_par.ped.empty() )
        info = read_plink(m_par.ped, m_gt);
    else if ( ! m_par.hmp.empty() )
        info = read_hapmap(m_par.hmp, m_gt);
    else if ( ! m_par.geno.empty() )
        info = read_genotype(m_par.geno, m_gt);

    if (info != 0) {
        m_gt.loc.clear();
        m_gt.ind.clear();
        m_gt.dat.clear();
    }

    std::cerr << "INFO: " << m_gt.ind.size() << " individuals and " << m_gt.loc.size() << " loci were observed\n";
}

void AppAssoc::load_phenotype()
{
    if ( m_par.pheno.empty() )
        return;

    std::cerr << "INFO: reading phenotype file...\n";

    int ret = read_phenotype(m_par.pheno, m_pt);

    if (ret != 0) {
        m_pt.phe.clear();
        m_pt.ind.clear();
    }

    std::cerr << "INFO: " << m_pt.ind.size() << " observations and " << m_pt.phe.size() << " phenotypes were observed\n";
}

void AppAssoc::load_covariate()
{
    if ( m_par.covar.empty() )
        return;

    std::cerr << "INFO: reading covariate file...\n";

    int ret = read_phenotype(m_par.covar, m_ct);

    if (ret != 0) {
        m_ct.phe.clear();
        m_ct.ind.clear();
    }

    std::cerr << "INFO: " << m_ct.ind.size() << " observations and " << m_ct.phe.size() << " covariates were observed\n";
}

void AppAssoc::load_kinship()
{
    if ( m_par.kin.empty() )
        return;

    std::cerr << "INFO: reading kinship file...\n";

    int info = read_square(m_par.kin, m_kin);

    if (info != 0) {
        m_kin.ind.clear();
        m_kin.dat.clear();
    }

    std::cerr << "INFO: " << m_kin.ind.size() << " individuals were observed\n";
}

void AppAssoc::merge_data()
{
    bool has_kin = !m_kin.ind.empty();
    bool has_covar = !m_ct.phe.empty() && !m_ct.ind.empty();

    bool require_merge = m_gt.ind != m_pt.ind;

    if (!require_merge && has_covar)
        require_merge = m_gt.ind != m_ct.ind;

    if (!require_merge && has_kin)
        require_merge = m_gt.ind != m_kin.ind;

    if ( ! require_merge ) {
        m_obs.resize(m_gt.ind.size());
        std::iota(m_obs.begin(), m_obs.end(), size_t(0));
        return;
    }

    std::cerr << "INFO: merging data set by individual...\n";

    if (has_kin && m_gt.ind != m_kin.ind) {
        auto ind = intersect(m_gt.ind, m_kin.ind);
        vector<size_t> gi, ki;
        for (auto &e : ind) {
            gi.push_back( index(m_gt.ind, e) );
            ki.push_back( index(m_kin.ind, e) );
        }
        filter_ind(gi, m_gt);
        filter_ind(ki, m_kin);
    }

    auto ind = intersect(m_gt.ind, m_pt.ind);
    if ( has_covar )
        ind = intersect(ind, m_ct.ind);

    vector<size_t> gi, pi, ci;
    for (auto itr = m_pt.ind.begin(); itr != m_pt.ind.end(); ++itr) {
        if ( std::binary_search(ind.begin(), ind.end(), *itr) ) {
            pi.push_back( itr - m_pt.ind.begin() );
            gi.push_back( index(m_gt.ind, *itr) );
            if ( has_covar )
                ci.push_back( index(m_ct.ind, *itr) );
        }
    }

    m_obs.swap(gi);
    filter_ind(pi, m_pt);
    if ( has_covar )
        filter_ind(ci, m_ct);

    std::cerr << "INFO: after merging data set, there are " << ind.size() << " individuals\n";
}

void AppAssoc::parse_factor()
{
    vector< vector<double> > xenv, xblk;

    if ( ! m_pt.env.empty() ) {
        design(1, factor(m_pt.env), xenv);
        m_addcov.insert(m_addcov.end(), xenv.begin(), xenv.end());
        if ( ! m_par.no_gxe )
            m_intcov.insert(m_intcov.end(), xenv.begin(), xenv.end());
    }

    if ( ! m_pt.blk.empty() ) {
        design(1, factor(m_pt.blk), xblk);
        if (xenv.empty()) {
            m_addcov.insert(m_addcov.end(), xblk.begin(), xblk.end());
        }
        else {
            design(3, factor(m_pt.env), xenv);
            for (auto &e : xenv) {
                for (auto v : xblk) {
                    std::transform(e.begin(), e.end(), v.begin(), v.begin(), std::multiplies<double>());
                    m_addcov.push_back(v);
                }
            }
        }
    }
}

void AppAssoc::LM()
{
    auto m = m_gt.loc.size();
    auto t = m_pt.phe.size();
    bool haploid = m_gt.ploidy == 1;

    std::ofstream ofs(m_par.out + ".lm.txt");
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << m_par.out << ".lm.txt" << "\n";
        return;
    }

    ofs << "Locus\tDF\tP\tR2";
    if ( ! m_intcov.empty() )
        ofs << "\tDFi\tPi\tR2i";
    ofs << "\n\n";

    for (size_t i = 0; i < t; ++i) {
        std::cerr << "INFO: processing phenotype: " << m_pt.phe[i] << "\n";

        vector<bool> nonan;
        for (auto e : m_pt.dat[i])
            nonan.push_back( std::isfinite(e) );

        auto obs = m_obs;
        auto y = m_pt.dat[i];
        auto addcov = m_addcov;
        auto intcov = m_intcov;
        addcov.insert(addcov.end(), m_ct.dat.begin(), m_ct.dat.end());

        if ( contain(nonan, false) ) {
            obs = subset(obs, nonan);
            y = subset(y, nonan);
            for (auto &e : addcov)
                e = subset(e, nonan);
            for (auto &e : intcov)
                e = subset(e, nonan);
        }

        if (y.size() < 10) {
            std::cerr << "ERROR: not enough valid observations: " << y.size() << "\n";
            continue;
        }

        vector<double> result;
        assoc_LM(haploid, m_gt.dat, obs, y, addcov, intcov, result);

        vector<double> result_ps(result.begin() + m, result.begin() + 2*m);
        std::sort(result_ps.begin(), result_ps.end());

        auto fdr = fdr_threshold(m_par.alpha, result_ps);
        ofs << ">" << m_pt.phe[i] << " BON:" << m_par.alpha / m << "|FDR:" << fdr << "\n";

        for (size_t j = 0; j < m; ++j) {
            ofs << m_gt.loc[j] << "\t" << result[j] << "\t" << result[m+j] << "\t" << result[2*m+j];
            if ( ! intcov.empty() )
                ofs << "\t" << result[3*m+j] << "\t" << result[4*m+j] << "\t" << result[5*m+j];
            ofs << "\n";
        }
        ofs << "\n";
    }
}

void AppAssoc::LMM()
{
    auto m = m_gt.loc.size();
    auto t = m_pt.phe.size();
    bool haploid = m_gt.ploidy == 1;

    std::ofstream ofs(m_par.out + ".lmm.txt");
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << m_par.out << ".lmm.txt" << "\n";
        return;
    }

    ofs << "Locus\tDF\tP\n\n";

    for (size_t i = 0; i < t; ++i) {
        std::cerr << "INFO: processing phenotype: " << m_pt.phe[i] << "\n";

        vector<bool> nonan;
        for (auto e : m_pt.dat[i])
            nonan.push_back( std::isfinite(e) );

        auto obs = m_obs;
        auto y = m_pt.dat[i];
        auto addcov = m_addcov;
        auto intcov = m_intcov;
        addcov.insert(addcov.end(), m_ct.dat.begin(), m_ct.dat.end());

        if ( contain(nonan, false) ) {
            obs = subset(obs, nonan);
            y = subset(y, nonan);
            for (auto &e : addcov)
                e = subset(e, nonan);
            for (auto &e : intcov)
                e = subset(e, nonan);
        }

        if (y.size() < 10) {
            std::cerr << "ERROR: not enough valid observations: " << y.size() << "\n";
            continue;
        }

        vector<double> result;
        assoc_LMM(haploid, m_gt.dat, obs, y, addcov, intcov, m_kin.dat, result);

        vector<double> result_ps(result.begin() + m, result.begin() + 2*m);
        std::sort(result_ps.begin(), result_ps.end());

        auto fdr = fdr_threshold(m_par.alpha, result_ps);
        ofs << ">" << m_pt.phe[i] << " BON:" << m_par.alpha / m << "|FDR:" << fdr << "\n";

        for (size_t j = 0; j < m; ++j)
            ofs << m_gt.loc[j] << "\t" << result[j] << "\t" << result[m+j] << "\n";
        ofs << "\n";
    }
}

void AppAssoc::RTM()
{
    auto m = m_gt.loc.size();
    auto t = m_pt.phe.size();
    bool haploid = m_gt.ploidy == 1;

    int multtest = 0;
    if (m_par.multtest == "BON")
        multtest = 1;
    else if (m_par.multtest == "FDR")
        multtest = 2;

    std::ofstream ofs(m_par.out + ".rtm.qtl.txt");
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << m_par.out << ".rtm.qtl.txt" << "\n";
        return;
    }

    std::ofstream ofs_fit(m_par.out + ".rtm.fit.txt");
    if ( ! ofs_fit )
        std::cerr << "ERROR: can't open file: " << m_par.out << ".rtm.fit.txt" << "\n";

    for (size_t i = 0; i < t; ++i) {
        std::cerr << "INFO: processing phenotype: " << m_pt.phe[i] << "\n";

        vector<bool> nonan;
        for (auto e : m_pt.dat[i])
            nonan.push_back( std::isfinite(e) );

        auto obs = m_obs;
        auto y = m_pt.dat[i];
        auto addcov = m_addcov;
        auto intcov = m_intcov;
        addcov.insert(addcov.end(), m_ct.dat.begin(), m_ct.dat.end());

        if ( contain(nonan, false) ) {
            obs = subset(obs, nonan);
            y = subset(y, nonan);
            for (auto &e : addcov)
                e = subset(e, nonan);
            for (auto &e : intcov)
                e = subset(e, nonan);
        }

        if ( y.size() < 10 ) {
            std::cerr << "WARN: not enough valid observations: " << y.size() << "\n";
            continue;
        }

        vector<size_t> result;

        if (m_par.preselect > 0.0 && m_par.preselect < 1.0) {
            vector<double> st;
            assoc_LM(haploid, m_gt.dat, obs, y, addcov, intcov, st);

            vector<size_t> idx;
            for (size_t j = 0; j < m; ++j) {
                if (st[m+j] <= m_par.preselect || ( ! intcov.empty() && st[4*m+j] <= m_par.preselect))
                    idx.push_back(j);
            }

            auto geno = subset(m_gt.dat, idx);
            std::cerr << "INFO: after pre-selection, there are " << idx.size() << " loci\n";

            assoc_RTM(multtest, m_par.maxqtl, m_par.alpha, m_par.maxrsq,
                      haploid, geno, obs, y, addcov, intcov, result);

            subset(idx, result).swap(result);
        }
        else {
            assoc_RTM(multtest, m_par.maxqtl, m_par.alpha, m_par.maxrsq,
                      haploid, m_gt.dat, obs, y, addcov, intcov, result);
        }

        ofs << ">" << m_pt.phe[i] << "\n";
        for (auto j : result)
            ofs << m_gt.loc[j] << "\n";
        ofs << "\n";

        auto s = fit(i, result);
        ofs_fit << ">" << m_pt.phe[i] << "\n" << s << "\n";
    }
}

string AppAssoc::fit(int phe, const vector<size_t> &loc) const
{
    auto y = m_pt.dat[phe];
    auto n = y.size();

    vector<size_t> obs;
    for (size_t i = 0; i < n; ++i) {
        if ( std::isfinite(y[i]) )
            obs.push_back(i);
    }

    auto nv = obs.size();
    if (nv != n)
        subset(y, obs).swap(y);

    vector<string> env;
    if ( ! m_pt.env.empty() )
        subset(m_pt.env, obs).swap(env);

    vector<string> blk;
    if ( ! m_pt.blk.empty() )
        subset(m_pt.blk, obs).swap(blk);

    GLM glm;

    if ( ! env.empty() )
        glm.add_main("_ENV_", env);

    if ( ! blk.empty() ) {
        if (env.empty())
            glm.add_main("_BLK_", blk);
        else
            glm.add_nested("_BLK_(_ENV_)", blk, env);
    }

    if ( !m_ct.phe.empty() && !m_ct.ind.empty() ) {
        auto num_covar = m_ct.phe.size();
        for (size_t i = 0; i < num_covar; ++i) {
            if (nv == n)
                glm.add_reg(m_ct.phe[i], m_ct.dat[i]);
            else
                glm.add_reg(m_ct.phe[i], subset(m_ct.dat[i], obs));
        }
    }

    subset(m_obs, obs).swap(obs);

    for (auto j : loc) {
        vector<string> gs;
        if (m_gt.ploidy == 1) {
            for (auto i : obs) {
                auto a = m_gt.dat[j][i];
                if (a == 0)
                    gs.push_back("?");
                else
                    gs.push_back(m_gt.allele[j][a-1]);
            }
        }
        else {
            for (auto i : obs) {
                auto a = m_gt.dat[j][i*2];
                auto b = m_gt.dat[j][i*2+1];
                if (a > b)
                    std::swap(a, b);
                string sab;
                if (a == 0)
                    sab.push_back('?');
                else
                    sab.append(m_gt.allele[j][a-1]);
                sab.push_back(':');
                if (b == 0)
                    sab.push_back('?');
                else
                    sab.append(m_gt.allele[j][b-1]);
                gs.push_back(sab);
            }
        }

        glm.add_main(m_gt.loc[j], gs);

        if ( ! env.empty() && ! m_par.no_gxe )
            glm.add_crossed(m_gt.loc[j] + "*_ENV_", gs, env);
    }

    auto atab = m_par.sstype == 3 ? glm.solve3(y) : glm.solve1(y);
    auto est = glm.coeff(y);

    std::ostringstream oss;

    oss << "<ANOVA>\nSource\tDF\tSS\tMS\tF\tp\n";
    for (size_t m = atab.names.size(), i = 0; i < m; ++i)
        oss << atab.names[i] << "\t" << atab.df[i] << "\t" << atab.ss[i] << "\t"
            << atab.ms[i] << "\t" << atab.f[i] << "\t" << atab.p[i] << "\n";
    oss << "Error\t" << atab.error[0] << "\t" << atab.error[1] << "\t" << atab.error[2] << "\n";
    oss << "Total\t" << atab.total[0] << "\t" << atab.total[1] << "\n";
    oss << "<SOLUTION>\nParameter\tEstimate\n";
    for (size_t m = est.params.size(), i = 0; i < m; ++i)
        oss << est.params[i] << "\t" << est.coeffs[i] << "\n";

    return oss.str();
}
