#include <memory>
#include <iostream>
#include <utility>
#include "appgsc.h"
#include "cmdline.h"
#include "plinkio.h"
#include "hapmapio.h"
#include "vcfio.h"
#include "rtmio.h"
#include "lapack.h"
#include "util.h"
#include "stat.h"

namespace {

template<typename T>
int count_shared_allele(T a, T b, T c, T d)
{
    if (a == b)
        return static_cast<int>(a == c) + (a == d);

    if (c == d)
        return static_cast<int>(c == a) + (c == b);

    return static_cast<int>(a == c) + (a == d) + (b == c) + (b == d);
}

std::pair<size_t,size_t> match_kernel_haplo(size_t n, const allele_t *x, const allele_t *y)
{
    size_t a = 0, b = 0;
    for (size_t i = 0; i < n; ++i) {
        if (x[i] != 0 && y[i] != 0) {
            ++a;
            if (x[i] == y[i])
                ++b;
        }
    }
    return std::make_pair(a, b);
}

std::pair<size_t,size_t> match_kernel_diplo(size_t n, const allele_t *x, const allele_t *y)
{
    size_t a = 0, b = 0;
    for (size_t i = 0; i < n; ++i) {
        size_t j1 = i*2, j2 = i*2+1;
        if (x[j1] != 0 && x[j2] != 0 && y[j1] != 0 && y[j2] != 0) {
            a += 2;
            b += count_shared_allele(x[j1], x[j2], y[j1], y[j2]);
        }
    }
    return std::make_pair(a, b);
}

} // namespace

int AppGSC::run(int argc, char *argv[])
{
    auto cmd = std::make_shared<CmdLine>("gsc [options]");

    cmd->add("--vcf", "VCF file", "");
    cmd->add("--ped", "PLINK PED/MAP file prefix", "");
    cmd->add("--hmp", "HapMap file", "");
    cmd->add("--geno", "genotype data file", "");
    cmd->add("--out", "output file prefix", "appgsc.out");
    cmd->add("--top", "number of eigenvectors", "10");

    if (argc < 2) {
        cmd->help();
        return 1;
    }

    cmd->parse(argc, argv);

    m_par.vcf = cmd->get("--vcf");
    m_par.ped = cmd->get("--ped");
    m_par.hmp = cmd->get("--hmp");
    m_par.geno = cmd->get("--geno");
    m_par.out = cmd->get("--out");
    m_par.top = number<int>(cmd->get("--top"));

    if (m_par.top < 1) {
        m_par.top = 10;
        std::cerr << "WARNING: invalid argument: --top " << cmd->get("--top") << "\n";
    }

    cmd.reset();

    return perform();
}

int AppGSC::perform()
{
    load_gentoype();

    if ( m_gt.loc.empty() || m_gt.ind.size() < 2 ) {
        std::cerr << "ERROR: not enough observations\n";
        return 1;
    }

    calc_gsc_matrix(m_gsc.dat);
    m_gsc.ind = m_gt.ind;
    m_evec.ind = m_gt.ind;

    vector<double> eval;
    vector< vector<double> > evec;
    eigen_decomposition(eval, evec);

    size_t top = m_par.top;
    if (top > eval.size())
        top = eval.size();

    double toteval = sum(eval);
    double cumeval = 0.0;
    std::cerr << "INFO: the top " << top << " eigenvalues are:\n";
    for (size_t i = 0; i < top; ++i) {
        m_evec.phe.push_back("ev"+std::to_string(i+1));
        m_evec.dat.push_back(evec[i]);
        cumeval += eval[i];
        std::cerr << "  " << i + 1 << " " << eval[i] << " " << 100*cumeval/toteval << "\n";
    }

    int info = write_phenotype(m_evec, m_par.out + ".evec.txt");
    if (info != 0)
        return 1;

    info = write_square(m_gsc, m_par.out + ".gsc.txt");
    if (info != 0)
        return 2;

    return 0;
}

void AppGSC::load_gentoype()
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

void AppGSC::calc_gsc_matrix(vector<vector<double> > &x) const
{
    static const size_t nb = 1000;

    size_t m = 0;
    size_t n = m_gt.ind.size();
    vector< vector<size_t> > z(n, vector<size_t>(n, 0));

    x.assign(n, vector<double>(n, 0.0));

    if (m_gt.ploidy == 1) {
        vector< vector<allele_t> > dat(n, vector<allele_t>(nb));

        for (auto &v : m_gt.dat) {
            if (m < nb) {
                for (size_t i = 0; i < n; ++i)
                    dat[i][m] = v[i];
                ++m;
                continue;
            }
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    if (j <= i)
                        continue;
                    auto p = match_kernel_haplo(m, dat[i].data(), dat[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }
            for (size_t i = 0; i < n; ++i)
                dat[i][0] = v[i];
            m = 1;
        }

        if (m > 0) {
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i + 1; j < n; ++j) {
                    auto p = match_kernel_haplo(m, dat[i].data(), dat[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }
        }
    }
    else {
        vector< vector<allele_t> > dat(n, vector<allele_t>(nb*2));

        for (auto &v : m_gt.dat) {
            if (m < nb) {
                for (size_t i = 0; i < n; ++i) {
                    dat[i][m*2] = v[i*2];
                    dat[i][m*2+1] = v[i*2+1];
                }
                ++m;
                continue;
            }
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    if (j <= i)
                        continue;
                    auto p = match_kernel_diplo(m, dat[i].data(), dat[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }
            for (size_t i = 0; i < n; ++i) {
                dat[i][0] = v[i*2];
                dat[i][1] = v[i*2+1];
            }
            m = 1;
        }

        if (m > 0) {
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i+1; j < n; ++j) {
                    auto p = match_kernel_diplo(m, dat[i].data(), dat[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }
        }
    }

    for (size_t i = 0; i < n; ++i) {
        x[i][i] = 1.0;
        for (size_t j = i+1; j < n; ++j) {
            if (z[i][j] != 0) {
                auto a = static_cast<double>(z[j][i]) / z[i][j];
                x[i][j] = x[j][i] = a;
            }
        }
    }
}

void AppGSC::eigen_decomposition(vector<double> &eval, vector<vector<double> > &evec) const
{
    vector<double> a;
    for (auto &v : m_gsc.dat)
        a.insert(a.end(), v.begin(), v.end());

    auto n = m_gsc.ind.size();
    vector<double> w(n), z(n*n);
    vector<bint> sup(2*n);
    bint m = 0;

    la_dsyevr('V', 'A', 'U', n, a.data(), n, 0.0, 0.0, 0, 0, 0.0, &m, w.data(), z.data(), n, sup.data());

    std::reverse(w.begin(), w.end());

    eval.swap(w);
    evec.clear();

    for (size_t j = 0; j < n; ++j) {
        auto jj = n - 1 - j;
        vector<double> v(n);
        for (size_t i = 0; i < n; ++i)
            v[i] = z[jj*n+i];
        evec.push_back(v);
    }
}
