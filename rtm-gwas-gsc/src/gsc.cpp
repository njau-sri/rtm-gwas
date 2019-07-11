#include <fstream>
#include <iostream>
#include <algorithm>
#include "version.h"
#include "cmdline.h"
#include "vcf.h"
#include "lapack.h"
#include "pheno.h"


using std::size_t;


namespace {


struct Parameter
{
    std::string vcf;
    std::string grm;
    std::string out;
    int top = 10;
    bool openmp = false;
} par;


template<typename T>
size_t count_shared(T a, T b, T c, T d)
{
    if (a == b)
        return (a == c ? 1 : 0) + (a == d ? 1 : 0);

    if (c == d)
        return (c == a ? 1 : 0) + (c == b ? 1 : 0);

    return (a == c ? 1 : 0) + (a == d ? 1 : 0) + (b == c ? 1 : 0) + (b == d ? 1 : 0);
}

void calc_gsc1(size_t n, const std::vector<allele_t> &x, const std::vector<allele_t> &y, size_t &a, size_t &b)
{
    for (size_t i = 0; i < n; ++i) {
        if (x[i] && y[i]) {
            ++a;
            if (x[i] == y[i])
                ++b;
        }
    }
}

void calc_gsc2(size_t n, const std::vector<allele_t> &x, const std::vector<allele_t> &y, size_t &a, size_t &b)
{
    for (size_t i = 0; i < n; ++i) {
        auto j = i * 2;
        if (x[j] && x[j+1] && y[j] && y[j+1]) {
            a += 2;
            b += count_shared(x[j], x[j+1], y[j], y[j+1]);
        }
    }
}

void calc_gsc_matrix(const Genotype &gt, std::vector< std::vector<double> > &x)
{
    static const size_t m = 10000;

    auto n = gt.ind.size();

    x.assign(n, std::vector<double>(n,1));

    std::vector< std::vector<size_t> > z(n, std::vector<size_t>(n,0));

    if (gt.ploidy == 1) {
        std::vector< std::vector<allele_t> > dat(n, std::vector<allele_t>(m));

        size_t k = 0;
        for (auto &v : gt.dat) {
            if (k < m) {
                for (size_t i = 0; i < n; ++i)
                    dat[i][k] = v[i];
                ++k;
                continue;
            }

            for (size_t i = 0; i < n; ++i)
                for (size_t j = i + 1; j < n; ++j)
                    calc_gsc1(m, dat[i], dat[j], z[i][j], z[j][i]);

            for (size_t i = 0; i < n; ++i)
                dat[i][0] = v[i];
            k = 1;
        }

        if (k > 0) {
            for (size_t i = 0; i < n; ++i)
                for (size_t j = i + 1; j < n; ++j)
                    calc_gsc1(k, dat[i], dat[j], z[i][j], z[j][i]);
        }
    }
    else {
        std::vector< std::vector<allele_t> > dat(n, std::vector<allele_t>(m * 2));

        size_t k = 0;
        for (auto &v : gt.dat) {
            if (k < m) {
                for (size_t i = 0; i < n; ++i) {
                    dat[i][k*2] = v[i*2];
                    dat[i][k*2+1] = v[i*2+1];
                }
                ++k;
                continue;
            }

            for (size_t i = 0; i < n; ++i)
                for (size_t j = i + 1; j < n; ++j)
                    calc_gsc2(m, dat[i], dat[j], z[i][j], z[j][i]);

            for (size_t i = 0; i < n; ++i) {
                dat[i][0] = v[i*2];
                dat[i][1] = v[i*2+1];
            }
            k = 1;
        }

        if (k > 0) {
            for (size_t i = 0; i < n; ++i)
                for (size_t j = i + 1; j < n; ++j)
                    calc_gsc2(k, dat[i], dat[j], z[i][j], z[j][i]);
        }
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            auto a = static_cast<double>(z[j][i]) / z[i][j];
            x[i][j] = x[j][i] = a;
        }
    }
}

void calc_gsc_matrix_omp(const Genotype &gt, std::vector< std::vector<double> > &x)
{
    static const size_t m = 10000;

    auto n = gt.ind.size();

    x.assign(n, std::vector<double>(n,1));

    std::vector< std::vector<size_t> > z(n, std::vector<size_t>(n,0));

    std::vector< std::pair<size_t,size_t> > pidx;
    pidx.reserve(n*(n-1)/2);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j)
            pidx.emplace_back(i,j);
    }

    auto p = pidx.size();

    if (gt.ploidy == 1) {
        std::vector< std::vector<allele_t> > dat(n, std::vector<allele_t>(m));

        size_t k = 0;
        for (auto &v : gt.dat) {
            if (k < m) {
                for (size_t i = 0; i < n; ++i)
                    dat[i][k] = v[i];
                ++k;
                continue;
            }

            #pragma omp parallel for
            for (size_t l = 0; l < p; ++l) {
                auto i = pidx[l].first;
                auto j = pidx[l].second;
                calc_gsc1(m, dat[i], dat[j], z[i][j], z[j][i]);
            }

            for (size_t i = 0; i < n; ++i)
                dat[i][0] = v[i];
            k = 1;
        }

        if (k > 0) {
            #pragma omp parallel for
            for (size_t l = 0; l < p; ++l) {
                auto i = pidx[l].first;
                auto j = pidx[l].second;
                calc_gsc1(k, dat[i], dat[j], z[i][j], z[j][i]);
            }
        }
    }
    else {
        std::vector< std::vector<allele_t> > dat(n, std::vector<allele_t>(m*2));

        size_t k = 0;
        for (auto &v : gt.dat) {
            if (k < m) {
                for (size_t i = 0; i < n; ++i) {
                    dat[i][k*2] = v[i*2];
                    dat[i][k*2+1] = v[i*2+1];
                }
                ++k;
                continue;
            }

            #pragma omp parallel for
            for (size_t l = 0; l < p; ++l) {
                auto i = pidx[l].first;
                auto j = pidx[l].second;
                calc_gsc2(m, dat[i], dat[j], z[i][j], z[j][i]);
            }

            for (size_t i = 0; i < n; ++i) {
                dat[i][0] = v[i*2];
                dat[i][1] = v[i*2+1];
            }
            k = 1;
        }

        if (k > 0) {
            #pragma omp parallel for
            for (size_t l = 0; l < p; ++l) {
                auto i = pidx[l].first;
                auto j = pidx[l].second;
                calc_gsc2(k, dat[i], dat[j], z[i][j], z[j][i]);
            }
        }
    }

    #pragma omp parallel for
    for (size_t l = 0; l < p; ++l) {
        auto i = pidx[l].first;
        auto j = pidx[l].second;
        auto a = static_cast<double>(z[j][i]) / z[i][j];
        x[i][j] = x[j][i] = a;
    }
}

int eigen(const std::vector< std::vector<double> > &mat,
          std::vector<double> &eval,
          std::vector< std::vector<double> > &evec)
{
    std::vector<double> a;
    for (auto &v : mat)
        a.insert(a.end(), v.begin(), v.end());

    auto n = mat.size();
    bint m = 0;
    std::vector<bint> sup(n*2);
    std::vector<double> w(n), z(n*n);

    int info = C_dsyevr('V', 'A', 'U', n, a.data(), n, 0, 0, 0, 0, 0, &m, w.data(), z.data(), n, sup.data());

    if (info != 0)
        return 1;

    std::reverse(w.begin(), w.end());

    eval.swap(w);
    evec.clear();

    for (size_t j = 0; j < n; ++j) {
        auto k = n - 1 - j;
        std::vector<double> v(n);
        for (size_t i = 0; i < n; ++i)
            v[i] = z[k*n+i];
        evec.push_back(v);
    }

    return 0;
}


} // namespace


int rtm_gwas_gsc(int argc, char *argv[])
{
    std::cerr << "RTM-GWAS " RTM_GWAS_VERSION_STRING " GSC (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd;

    cmd.add("--vcf", "VCF file", "");
    cmd.add("--grm", "genetic relationship matrix file", "");
    cmd.add("--out", "output file", "gsc.out");
    cmd.add("--top", "number of eigenvectors", "10");
    cmd.add("--openmp", "enable OpenMP multithreading");

    cmd.parse(argc, argv);

    if (argc < 2) {
        cmd.show();
        return 1;
    }

    par.vcf = cmd.get("--vcf");
    par.grm = cmd.get("--grm");
    par.out = cmd.get("--out");
    par.top = std::stoi(cmd.get("--top"));
    par.openmp = cmd.has("--openmp");

    Genotype gt;
    SquareData sd;

    if ( par.grm.empty() ) {
        std::cerr << "INFO: reading genotype file...\n";
        if (read_vcf(par.vcf, gt) != 0)
            return 1;
        std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

        if ( par.openmp )
            calc_gsc_matrix_omp(gt, sd.dat);
        else
            calc_gsc_matrix(gt, sd.dat);

        sd.ind = gt.ind;

        write_square(sd, par.out + ".mat");
    }
    else {
        std::cerr << "INFO: reading genetic relationship matrix file...\n";
        if (read_square(par.grm, sd) != 0)
            return 1;
        std::cerr << "INFO: " << sd.ind.size() << " individuals\n";
    }

    auto m = static_cast<size_t>(par.top);

    Covariate ct;
    std::vector<double> eval;

    if (eigen(sd.dat, eval, ct.dat) != 0)
        return 1;

    ct.ind = sd.ind;

    for (size_t i = 0; i < m; ++i)
        ct.phe.push_back("EV" + std::to_string(i+1));

    ct.dat.resize(m);

    write_covar(ct, par.out + ".evec");

    std::ofstream ofs(par.out + ".eval");
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << par.out << ".eval\n";
        return 1;
    }

    double sum = 0;
    for (auto e : eval) {
        if (e > 0)
            sum += e;
    }

    for (size_t i = 0; i < m; ++i)
        ofs << eval[i] << "\t" << 100 * eval[i] / sum << "\n";

    return 0;
}
