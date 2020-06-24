#include "rtmgwasgsc.h"

#include "types.h"
#include "print.h"
#include "cmdline.h"
#include "vcf.h"
#include "pheno.h"
#include "openmp.h"
#include "matrix.h"
#include "cfile.h"

#ifndef RTM_GWAS_VERSION
#define RTM_GWAS_VERSION "unknown"
#endif // RTM_GWAS_VERSION

namespace {

    struct Parameter
    {
        std::string vcf;
        std::string grm;
        std::string out;
        int top = 10;
        int thread = 0;
    } par;

    void show_help()
    {
        const char *help =
            "Usage: rtm-gwas-gsc [options...]\n"
            "Options:\n"
            "  --vcf <>         SNP genotype file in VCF format\n"
            "  --out <gsc.out>  output file prefix\n"
            "  --top <10>       number of eigenvectors\n"
            "  --grm <>         genetic relationship matrix file\n"
            "  --thread <0>     set the number of threads\n"
            "\n";
        eprint(help);
    }

    void parse_cmdline(int argc, char *argv[])
    {
        CmdLine cmd;

        cmd.add("--vcf", "");
        cmd.add("--out", "gsc.out");
        cmd.add("--top", "10");
        cmd.add("--grm", "");
        cmd.add("--thread", "0");

        cmd.parse(argc, argv);

        par.vcf = cmd.get("--vcf");
        par.grm = cmd.get("--grm");
        par.out = cmd.get("--out");
        par.top = std::stoi(cmd.get("--top"));
        par.thread = std::stoi(cmd.get("--thread"));
    }

    std::pair<isize_t, isize_t> count_shared_hap(isize_t n, const int *x, const int *y)
    {
        isize_t a = 0, b = 0;
        for (isize_t i = 0; i < n; ++i) {
            if (x[i] && y[i]) {
                ++a;
                if (x[i] == y[i])
                    ++b;
            }
        }
        return { a, b };
    }

    int count_shared_dip_kernel(int a, int b, int c, int d)
    {
        if (a == b) return (a == c) + (a == d);
        if (c == d) return (c == a) + (c == b);
        return (a == c) + (a == d) + (b == c) + (b == d);
    }

    std::pair<isize_t, isize_t> count_shared_dip(isize_t n, const int *x, const int *y)
    {
        isize_t a = 0, b = 0;
        for (isize_t i = 0; i < n; ++i) {
            isize_t j = i * 2;
            if (x[j] && x[j+1] && y[j] && y[j+1]) {
                a += 2;
                b += count_shared_dip_kernel(x[j], x[j+1], y[j], y[j+1]);
            }
        }
        return { a, b };
    }

    void count_shared_matrix_hap(const Genotype &gt, std::vector< std::vector<isize_t> > &z)
    {
        isize_t n = length(gt.ind);

        z.assign(n, std::vector<isize_t>(n, 0));

        const isize_t chunk_size = 10000;
        std::vector< std::vector<int> > chunk(n, std::vector<int>(chunk_size));

        isize_t k = 0;
        for (auto &v : gt.dat) {
            for (isize_t i = 0; i < n; ++i)
                chunk[i][k] = v[i];

            if (++k < chunk_size)
                continue;

            for (isize_t i = 0; i < n; ++i) {
                for (isize_t j = i + 1; j < n; ++j) {
                    auto p = count_shared_hap(chunk_size, chunk[i].data(), chunk[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }

            k = 0;
        }

        if (k > 0) {
            for (isize_t i = 0; i < n; ++i) {
                for (isize_t j = i + 1; j < n; ++j) {
                    auto p = count_shared_hap(k, chunk[i].data(), chunk[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }
        }
    }

    void count_shared_matrix_dip(const Genotype &gt, std::vector< std::vector<isize_t> > &z)
    {
        isize_t n = length(gt.ind);

        z.assign(n, std::vector<isize_t>(n, 0));

        const isize_t chunk_size = 10000;
        std::vector< std::vector<int> > chunk(n, std::vector<int>(chunk_size*2));

        isize_t k = 0;
        for (auto &v : gt.dat) {
            for (isize_t i = 0; i < n; ++i) {
                chunk[i][k*2] = v[i*2];
                chunk[i][k*2+1] = v[i*2+1];
            }

            if (++k < chunk_size)
                continue;

            for (isize_t i = 0; i < n; ++i) {
                for (isize_t j = i + 1; j < n; ++j) {
                    auto p = count_shared_dip(chunk_size, chunk[i].data(), chunk[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }

            k = 0;
        }

        if (k > 0) {
            for (isize_t i = 0; i < n; ++i) {
                for (isize_t j = i + 1; j < n; ++j) {
                    auto p = count_shared_dip(k, chunk[i].data(), chunk[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }
        }
    }

    void count_shared_matrix_hap_omp(const Genotype &gt, std::vector< std::vector<isize_t> > &z)
    {
        isize_t n = length(gt.ind);

        z.assign(n, std::vector<isize_t>(n, 0));

        const isize_t chunk_size = 10000;
        std::vector< std::vector<int> > chunk(n, std::vector<int>(chunk_size));

        std::vector< std::pair<isize_t, isize_t> > pidx;

        isize_t npair = n*(n-1)/2;
        pidx.reserve(npair);

        for (isize_t i = 0; i < n; ++i) {
            for (isize_t j = i + 1; j < n; ++j)
                pidx.emplace_back(i, j);
        }

        isize_t k = 0;
        for (auto &v : gt.dat) {
            for (isize_t i = 0; i < n; ++i)
                chunk[i][k] = v[i];

            if (++k < chunk_size)
                continue;

            #pragma omp parallel for
            for (isize_t l = 0; l < npair; ++l) {
                isize_t i = pidx[l].first;
                isize_t j = pidx[l].second;
                auto p = count_shared_hap(chunk_size, chunk[i].data(), chunk[j].data());
                z[i][j] += p.first;
                z[j][i] += p.second;
            }

            k = 0;
        }

        if (k > 0) {
            #pragma omp parallel for
            for (isize_t l = 0; l < npair; ++l) {
                isize_t i = pidx[l].first;
                isize_t j = pidx[l].second;
                auto p = count_shared_hap(k, chunk[i].data(), chunk[j].data());
                z[i][j] += p.first;
                z[j][i] += p.second;
            }
        }
    }

    void count_shared_matrix_dip_omp(const Genotype &gt, std::vector< std::vector<isize_t> > &z)
    {
        isize_t n = length(gt.ind);

        z.assign(n, std::vector<isize_t>(n, 0));

        const isize_t chunk_size = 10000;
        std::vector< std::vector<int> > chunk(n, std::vector<int>(chunk_size*2));

        isize_t npair = n*(n-1)/2;

        std::vector< std::pair<isize_t, isize_t> > pidx;
        pidx.reserve(npair);
        for (isize_t i = 0; i < n; ++i) {
            for (isize_t j = i + 1; j < n; ++j)
                pidx.emplace_back(i, j);
        }

        isize_t k = 0;
        for (auto &v : gt.dat) {
            for (isize_t i = 0; i < n; ++i) {
                chunk[i][k*2] = v[i*2];
                chunk[i][k*2+1] = v[i*2+1];
            }

            if (++k < chunk_size)
                continue;

            #pragma omp parallel for
            for (isize_t l = 0; l < npair; ++l) {
                isize_t i = pidx[l].first;
                isize_t j = pidx[l].second;
                auto p = count_shared_dip(chunk_size, chunk[i].data(), chunk[j].data());
                z[i][j] += p.first;
                z[j][i] += p.second;
            }

            k = 0;
        }

        if (k > 0) {
            #pragma omp parallel for
            for (isize_t l = 0; l < npair; ++l) {
                isize_t i = pidx[l].first;
                isize_t j = pidx[l].second;
                auto p = count_shared_dip(k, chunk[i].data(), chunk[j].data());
                z[i][j] += p.first;
                z[j][i] += p.second;
            }
        }
    }

    int eigen_gsc(const std::vector< std::vector<double> > &a,
                  std::vector<double> &el,
                  std::vector< std::vector<double> > &ec)
    {
        isize_t n = length(a);

        mat b(n, n);
        for (isize_t i = 0; i < n; ++i)
            b.col(i) = mat::Map(a[i].data(), n, 1);

        mat eval, evec;
        int info = syevr(b, eval, evec);

        if (info != 0) {
            eprint("ERROR: syevr failed in eigen_gsc (%s:%d): %d\n", __FILE__, __LINE__, info);
            return 1;
        }

        el.resize(n);
        ec.assign(n, std::vector<double>(n));

        for (isize_t i = 0; i < n; ++i) {
            isize_t j = n - i - 1;
            el[i] = eval(j);
            mat::Map(ec[i].data(), n, 1) = evec.col(j);
        }

        return 0;
    }

} // namespace

int rtm_gwas_gsc(int argc, char *argv[])
{
    eprint("RTM-GWAS %s GSC (Built on %s %s)\n", RTM_GWAS_VERSION, __DATE__, __TIME__);

    if (argc < 2) {
        show_help();
        return 1;
    }

    parse_cmdline(argc, argv);

    if (par.thread > 0)
        call_omp_set_num_threads(par.thread);

    Genotype gt;
    SquareData sd;

    if (par.grm.empty()) {
        eprint("INFO: reading genotype file...\n");

        if (read_vcf(par.vcf, gt) != 0)
            return 1;

        isize_t n = length(gt.ind);
        isize_t m = length(gt.loc);

        eprint("INFO: %td individuals, %td loci\n", n, m);

        std::vector< std::vector<isize_t> > z;

        if (par.thread > 0) {
            if (gt.ploidy == 2)
                count_shared_matrix_dip_omp(gt, z);
            else
                count_shared_matrix_hap_omp(gt, z);
        }
        else {
            if (gt.ploidy == 2)
                count_shared_matrix_dip(gt, z);
            else
                count_shared_matrix_hap(gt, z);
        }

        sd.dat.assign(n, std::vector<double>(n, 1.0));
        for (isize_t i = 0; i < n; ++i) {
            for (isize_t j = i + 1; j < n; ++j) {
                double a = (double) z[j][i] / z[i][j];
                sd.dat[i][j] = sd.dat[j][i] = a;
            }
        }

        sd.ind = gt.ind;
        write_square(sd, par.out + ".mat");
    }
    else {
        eprint("INFO: reading genetic relationship matrix file...\n");

        if (read_square(par.grm, sd) != 0)
            return 1;

        eprint("INFO: %td individuals\n", length(sd.ind));
    }

    Covariate ct;
    std::vector<double> eval;

    if (eigen_gsc(sd.dat, eval, ct.dat) != 0)
        return 1;

    ct.ind = sd.ind;

    for (int i = 0; i < par.top; ++i)
        ct.phe.push_back("EV" + std::to_string(i+1));

    ct.dat.resize(par.top);

    write_covar(ct, par.out + ".evec");

    CFile file(par.out + ".eval", "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s.eval\n", par.out);
        return 1;
    }

    double sum = 0;
    for (auto e : eval) {
        if (e > 0)
            sum += e;
    }

    for (int i = 0; i < par.top; ++i)
        fprint(file, "EV%d\t%g\t%g\n", i+1, eval[i], 100 * eval[i] / sum);

    eprint("INFO: GSC has finished successfully\n");

    return 0;
}
