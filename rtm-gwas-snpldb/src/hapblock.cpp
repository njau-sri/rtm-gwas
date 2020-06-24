#include "hapblock.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <utility>
#include <algorithm>

#include "print.h"
#include "popgen.h"

namespace {

    struct GabrielBlock
    {
        int8_t llim;
        int8_t ulim;
        int16_t frac;
        int first;  // index of first snp
        int last;   // index of last snp
        int length;
    };

    //     BB        Bb      bb
    // AA  g[0][0]  g[0][1]  g[0][2]
    // Aa  g[1][0]  g[1][1]  g[1][2]
    // aa  g[2][0]  g[2][1]  g[2][2]
    // x,y  0:AA, 1:Aa, 2:aa, -1:missing
    void count_gtable(const std::vector<int8_t> &x, const std::vector<int8_t> &y, int g[3][3])
    {
        isize_t n = length(x);
        for (isize_t i = 0; i < n; ++i)
            if (x[i] != -1 && y[i] != -1)
                ++g[x[i]][y[i]];
    }

    // D' 95% confidence interval
    int est_dprime_ci(int maxit, double tol, int g[3][3], int &lower, int &upper)
    {
        using std::log;

        int n11 = 2 * g[0][0] + g[0][1] + g[1][0];
        int n12 = 2 * g[0][2] + g[0][1] + g[1][2];
        int n21 = 2 * g[2][0] + g[1][0] + g[2][1];
        int n22 = 2 * g[2][2] + g[2][1] + g[1][2];
        int ndh = g[1][1];

        int nn = n11 + n12 + n21 + n22 + ndh * 2;
        if (nn < 4)
            return 1;

        double p11 = 0.0, p12 = 0.0, p21 = 0.0, p22 = 0.0;

        if (ndh > 0) {
            double x = est_het(maxit, tol, n11, n12, n21, n22, ndh);
            double y1 = ndh * x;
            double y2 = ndh - y1;
            p11 = (n11 + y1) / nn;
            p12 = (n12 + y2) / nn;
            p21 = (n21 + y2) / nn;
            p22 = (n22 + y1) / nn;
        }
        else {
            p11 = (double) n11 / nn;
            p12 = (double) n12 / nn;
            p21 = (double) n21 / nn;
            p22 = (double) n22 / nn;
        }

        if (p11*p22 - p12*p21 < 0.0) {
            std::swap(p11, p12);
            std::swap(p21, p22);
            std::swap(n11, n12);
            std::swap(n21, n22);
        }

        double pA = p11 + p12;
        double pa = 1.0 - pA;
        double pB = p11 + p21;
        double pb = 1.0 - pB;

        double Dmax = std::min(pA*pb, pa*pB);

        if (p11 < 1e-10) p11 = 1e-10;
        if (p12 < 1e-10) p12 = 1e-10;
        if (p21 < 1e-10) p21 = 1e-10;
        if (p22 < 1e-10) p22 = 1e-10;
        double LL1 = n11*log(p11) + n12*log(p12) + n21*log(p21) + n22*log(p22) + ndh*log(p11*p22 + p12*p21);

        double tp = 0.0;
        double ls[101];
        double Dstep = Dmax / 100;

        for (int i = 0; i <= 100; ++i) {
            double q11 = pA*pB + i*Dstep;
            double q12 = pA - q11;
            double q21 = pB - q11;
            double q22 = pa - q21;
            if (i == 100) {
                if (q11 < 1e-10) q11 = 1e-10;
                if (q12 < 1e-10) q12 = 1e-10;
                if (q21 < 1e-10) q21 = 1e-10;
                if (q22 < 1e-10) q22 = 1e-10;
            }
            double LL2 = n11*log(q11) + n12*log(q12) + n21*log(q21) + n22*log(q22) + ndh*log(q11*q22 + q12*q21);
            double prob = std::exp(LL2 - LL1);
            ls[i] = prob;
            tp += prob;
        }

        double sp = 0.0;
        double tp5 = tp * 0.05;
        for (int i = 0; i <= 100; ++i) {
            sp += ls[i];
            if (sp > tp5 && sp - ls[i] < tp5) {
                lower = i - 1;
                break;
            }
        }

        sp = 0.0;
        for (int i = 100; i >= 0; --i) {
            sp += ls[i];
            if (sp > tp5 && sp - ls[i] < tp5) {
                upper = i + 1;
                break;
            }
        }

        return 0;
    }

    double count_informative_pair(int first, int last,
                                  const std::vector< std::vector<char> > &sp,
                                  const std::vector< std::vector<char> > &rp,
                                  const std::vector< std::vector<char> > &ss)
    {
        int s = 0;
        int r = 0;
        int n = last - first + 1;

        if (n > 4) {
            for (int i = first; i <= last; ++i) {
                for (int j = i + 1; j <= last; ++j) {
                    if (sp[i][j]) ++s;
                    if (rp[i][j]) ++r;
                }
            }
        }
        else if (n == 2) {
            if (ss[last][first]) ++s;
            if (rp[first][last]) ++r;
        }
        else { // 3 or 4
            for (int i = first; i <= last; ++i) {
                for (int j = i + 1; j <= last; ++j) {
                    if (ss[i][j]) ++s;
                    if (rp[i][j]) ++r;
                }
            }
        }

        int t = s + r;

        if (n == 2) {
            if (t < 1)
                return -1;
        }
        else if (n == 3) {
            if (t < 3)
                return -2;
        }
        else {
            if (t < 6)
                return -3;
        }

        return (double) s / t;
    }

    // NOTE: block length = stop - start + 1

    int find_gabriel_full(const HapBlockGabriel::Parameter &par,
                          const std::vector<int> &pos,
                          const std::vector< std::vector<int8_t> > &geno,
                          std::vector<GabrielBlock> &blk)
    {
        int n = length(pos);

        // D' CI lower
        std::vector< std::vector<int8_t> > lc(n, std::vector<int8_t>(n, -1));

        // D' CI upper
        std::vector< std::vector<int8_t> > uc(n, std::vector<int8_t>(n, -1));

        // strong pair
        std::vector< std::vector<char> > sp(n, std::vector<char>(n, 0));

        // recombination pair
        std::vector< std::vector<char> > rp(n, std::vector<char>(n, 0));

        // small strong pair, lower - 2 SNPs, upper - 3 or 4 SNPs
        std::vector< std::vector<char> > ss(n, std::vector<char>(n, 0));

        for (int i = 0; i < n; ++i) {
            int start = pos[i];
            for (int j = i + 1; j < n; ++j) {
                int stop = pos[j];

                if (stop <= start) {
                    eprint("ERROR: position must be in ascending order in find_gabriel_full (%s:%d): %d %d\n",
                           __FILE__, __LINE__, start, stop);
                    return 1;
                }

                if (stop - start + 1 > par.maxlen)
                    break;

                int g[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
                count_gtable(geno[i], geno[j], g);

                int llim = -1, ulim = -1;
                est_dprime_ci(par.maxit, par.tol, g, llim, ulim);

                bool strong = false, strong2 = false, strong34 = false;
                if (ulim >= par.ulim) {
                    if (llim > par.llim)
                        strong = true;
                    if (llim > 80)
                        strong2 = strong34 = true;
                    else if (llim > 50)
                        strong34 = true;
                }

                bool recomb = ulim < par.recomb;

                lc[i][j] = llim;
                uc[i][j] = ulim;

                sp[i][j] = strong;
                rp[i][j] = recomb;

                ss[i][j] = strong34;
                ss[j][i] = strong2;
            }
        }

        for (int i = 0; i < n; ++i) {
            int start = pos[i];
            for (int j = i + 1; j < n; ++j) {
                int stop = pos[j];

                int len = stop - start + 1;
                if (len > par.maxlen)
                    break;

                if (lc[i][j] < par.llim || uc[i][j] < par.ulim)
                    continue;

                int q = (int) (j - i + 1);
                if ((q == 2 && len > 20000) || (q == 3 && len > 30000))
                    continue;

                double frac = count_informative_pair(i, j, sp, rp, ss);

                if (frac - par.frac > std::numeric_limits<double>::epsilon()) {
                    GabrielBlock b;
                    b.llim = lc[i][j];
                    b.ulim = uc[i][j];
                    b.frac = (int16_t) (frac * 10000);
                    b.first = i;
                    b.last = j;
                    b.length = len;
                    blk.push_back(std::move(b));
                }
            }
        }

        return 0;
    }

    int find_gabriel_split(const HapBlockGabriel::Parameter &par,
                           const std::vector<int> &pos,
                           const std::vector< std::vector<int8_t> > &geno,
                           std::vector<GabrielBlock> &blk)
    {
        int n = length(pos);
        int w = par.batch;

        for (int i = 0; i < n; ++i) {
            int start = pos[i];
            for (int j = i + 1; j < n; ++j) {
                int stop = pos[j];

                if (stop <= start) {
                    eprint("ERROR: position must be in ascending order in find_gabriel_window (%s:%d): %d %d\n",
                           __FILE__, __LINE__, start, stop);
                    return 1;
                }

                if (stop - start + 1 > par.maxlen)
                    break;

                int wij = j - i + 1;
                if (w < wij)
                    w = wij;
            }
        }

        // D' CI lower
        std::vector< std::vector<int8_t> > lc(w, std::vector<int8_t>(w, -1));

        // D' CI upper
        std::vector< std::vector<int8_t> > uc(w, std::vector<int8_t>(w, -1));

        // strong pair
        std::vector< std::vector<char> > sp(w, std::vector<char>(w, 0));

        // recombination pair
        std::vector< std::vector<char> > rp(w, std::vector<char>(w, 0));

        // small strong pair, lower - 2 SNPs, upper - 3 or 4 SNPs
        std::vector< std::vector<char> > ss(w, std::vector<char>(w, 0));

        int curr = 0;

        for (;;) {
            for (int i = 0; i < w; ++i) {
                int ii = curr + i;
                int start = pos[ii];
                for (int j = i + 1; j < w; ++j) {
                    int jj = curr + j;
                    int stop = pos[jj];

                    if (stop - start + 1 > par.maxlen)
                        break;

                    int g[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
                    count_gtable(geno[ii], geno[jj], g);

                    int llim = -1, ulim = -1;
                    est_dprime_ci(par.maxit, par.tol, g, llim, ulim);

                    bool strong = false, strong2 = false, strong34 = false;
                    if (ulim >= par.ulim) {
                        if (llim > par.llim)
                            strong = true;
                        if (llim > 80)
                            strong2 = strong34 = true;
                        else if (llim > 50)
                            strong34 = true;
                    }

                    bool recomb = ulim < par.recomb;

                    lc[i][j] = llim;
                    uc[i][j] = ulim;

                    sp[i][j] = strong;
                    rp[i][j] = recomb;

                    ss[i][j] = strong34;
                    ss[j][i] = strong2;
                }
            }

            for (int i = 0; i < w; ++i) {
                int ii = curr + i;
                int start = pos[ii];
                for (int j = i + 1; j < w; ++j) {
                    int jj = curr + j;
                    int stop = pos[jj];

                    int len = stop - start + 1;
                    if (len > par.maxlen)
                        break;

                    if (lc[i][j] < par.llim || uc[i][j] < par.ulim)
                        continue;

                    int q = j - i + 1;
                    if ((q == 2 && len > 20000) || (q == 3 && len > 30000))
                        continue;

                    double frac = count_informative_pair(i, j, sp, rp, ss);

                    if (frac - par.frac > std::numeric_limits<double>::epsilon()) {
                        GabrielBlock b;
                        b.llim = lc[i][j];
                        b.ulim = uc[i][j];
                        b.frac = (int16_t) (frac * 10000);
                        b.first = ii;
                        b.last = jj;
                        b.length = len;
                        blk.push_back(std::move(b));
                    }
                }
            }

            if (curr + w >= n) {
                eprint("INFO: 100% complete\n");
                break;
            }

            int stop = pos[curr + w];
            for (int i = 1; i < w; ++i) {
                int start = pos[++curr];
                if (stop - start + 1 <= par.maxlen)
                    break;
            }

            if (curr + w >= n)
                w = n - curr;

            eprint("INFO: %g%% complete\n", 100.0 * curr / n);
        }

        return 0;
    }

    int find_gabriel_full_omp(const HapBlockGabriel::Parameter &par,
                              const std::vector<int> &pos,
                              const std::vector< std::vector<int8_t> > &geno,
                              std::vector<GabrielBlock> &blk)
    {
        int n = length(pos);

        std::vector< std::pair<int, int> > pidx;

        for (int i = 0; i < n; ++i) {
            int start = pos[i];
            for (int j = i + 1; j < n; ++j) {
                int stop = pos[j];

                if (stop <= start) {
                    eprint("ERROR: position must be in ascending order in find_gabriel_full_omp (%s:%d): %d %d\n",
                           __FILE__, __LINE__, start, stop);
                    return 1;
                }

                if (stop - start + 1 > par.maxlen)
                    break;

                pidx.emplace_back(i, j);
            }
        }

        // D' CI lower
        std::vector< std::vector<int8_t> > lc(n, std::vector<int8_t>(n, -1));

        // D' CI upper
        std::vector< std::vector<int8_t> > uc(n, std::vector<int8_t>(n, -1));

        // strong pair
        std::vector< std::vector<char> > sp(n, std::vector<char>(n, 0));

        // recombination pair
        std::vector< std::vector<char> > rp(n, std::vector<char>(n, 0));

        // small strong pair, lower - 2 SNPs, upper - 3 or 4 SNPs
        std::vector< std::vector<char> > ss(n, std::vector<char>(n, 0));

        isize_t m = length(pidx);

        #pragma omp parallel for
        for (isize_t k = 0; k < m; ++k) {
            int i = pidx[k].first;
            int j = pidx[k].second;

            int g[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
            count_gtable(geno[i], geno[j], g);

            int llim = -1, ulim = -1;
            est_dprime_ci(par.maxit, par.tol, g, llim, ulim);

            bool strong = false, strong2 = false, strong34 = false;
            if (ulim >= par.ulim) {
                if (llim > par.llim)
                    strong = true;
                if (llim > 80)
                    strong2 = strong34 = true;
                else if (llim > 50)
                    strong34 = true;
            }

            bool recomb = ulim < par.recomb;

            lc[i][j] = llim;
            uc[i][j] = ulim;

            sp[i][j] = strong;
            rp[i][j] = recomb;

            ss[i][j] = strong34;
            ss[j][i] = strong2;
        }

        std::vector<int16_t> vf(m, -1);

        #pragma omp parallel for schedule(dynamic)
        for (isize_t k = 0; k < m; ++k) {
            int i = pidx[k].first;
            int j = pidx[k].second;

            if (lc[i][j] < par.llim || uc[i][j] < par.ulim)
                continue;

            int len = pos[j] - pos[i] + 1;

            int q = j - i + 1;
            if ((q == 2 && len > 20000) || (q == 3 && len > 30000))
                continue;

            double frac = count_informative_pair(i, j, sp, rp, ss);

            if (frac - par.frac > std::numeric_limits<double>::epsilon())
                vf[k] = (int16_t) (frac * 10000);
        }

        usize_t reserve = 0;
        for (auto e : vf) if (e > 0) ++reserve;
        blk.reserve(reserve);

        for (isize_t k = 0; k < m; ++k) {
            if (vf[k] > 0) {
                int i = pidx[k].first;
                int j = pidx[k].second;
                GabrielBlock b;
                b.llim = lc[i][j];
                b.ulim = uc[i][j];
                b.frac = vf[k];
                b.first = i;
                b.last = j;
                b.length = pos[j] - pos[i] + 1;
                blk.push_back(std::move(b));
            }
        }

        return 0;
    }

    int find_gabriel_split_omp(const HapBlockGabriel::Parameter &par,
                               const std::vector<int> &pos,
                               const std::vector< std::vector<int8_t> > &geno,
                               std::vector<GabrielBlock> &blk)
    {
        int n = length(pos);
        int w = par.batch;

        for (int i = 0; i < n; ++i) {
            int start = pos[i];
            for (int j = i + 1; j < n; ++j) {
                int stop = pos[j];

                if (stop <= start) {
                    eprint("ERROR: position must be in ascending order in find_gabriel_split_omp (%s:%d): %d %d\n",
                           __FILE__, __LINE__, start, stop);
                    return 1;
                }

                if (stop - start + 1 > par.maxlen)
                    break;

                int wij = j - i + 1;
                if (w < wij)
                    w = wij;
            }
        }

        std::vector<int16_t> vf;
        std::vector< std::pair<int, int> > pidx;

        // D' CI lower
        std::vector< std::vector<int8_t> > lc(w, std::vector<int8_t>(w, -1));

        // D' CI upper
        std::vector< std::vector<int8_t> > uc(w, std::vector<int8_t>(w, -1));

        // strong pair
        std::vector< std::vector<char> > sp(w, std::vector<char>(w, 0));

        // recombination pair
        std::vector< std::vector<char> > rp(w, std::vector<char>(w, 0));

        // small strong pair, lower - 2 SNPs, upper - 3 or 4 SNPs
        std::vector< std::vector<char> > ss(w, std::vector<char>(w, 0));

        int curr = 0;

        for (;;) {
            pidx.clear();

            for (int i = 0; i < w; ++i) {
                int ii = curr + i;
                int start = pos[ii];
                for (int j = i + 1; j < w; ++j) {
                    int jj = curr + j;
                    int stop = pos[jj];

                    if (stop - start + 1 > par.maxlen)
                        break;

                    pidx.emplace_back(i, j);
                }
            }

            isize_t m = length(pidx);

            #pragma omp parallel for
            for (isize_t k = 0; k < m; ++k) {
                int i = pidx[k].first;
                int j = pidx[k].second;

                int ii = curr + i;
                int jj = curr + j;

                int g[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
                count_gtable(geno[ii], geno[jj], g);

                int llim = -1, ulim = -1;
                est_dprime_ci(par.maxit, par.tol, g, llim, ulim);

                bool strong = false, strong2 = false, strong34 = false;
                if (ulim >= par.ulim) {
                    if (llim > par.llim)
                        strong = true;
                    if (llim > 80)
                        strong2 = strong34 = true;
                    else if (llim > 50)
                        strong34 = true;
                }

                bool recomb = ulim < par.recomb;

                lc[i][j] = llim;
                uc[i][j] = ulim;

                sp[i][j] = strong;
                rp[i][j] = recomb;

                ss[i][j] = strong34;
                ss[j][i] = strong2;
            }

            vf.assign(m, -1);

            #pragma omp parallel for schedule(dynamic)
            for (isize_t k = 0; k < m; ++k) {
                int i = pidx[k].first;
                int j = pidx[k].second;

                if (lc[i][j] < par.llim || uc[i][j] < par.ulim)
                    continue;

                int ii = curr + i;
                int jj = curr + j;
                int len = pos[jj] - pos[ii] + 1;

                int q = j - i + 1;
                if ((q == 2 && len > 20000) || (q == 3 && len > 30000))
                    continue;

                double frac = count_informative_pair(i, j, sp, rp, ss);

                if (frac - par.frac > std::numeric_limits<double>::epsilon())
                    vf[k] = (int16_t) (frac * 10000);
            }

            auto reserve = blk.size();
            for (auto e : vf) if (e > 0) ++reserve;
            if (blk.capacity() < reserve)
                blk.reserve(reserve);

            for (isize_t k = 0; k < m; ++k) {
                if (vf[k] > 0) {
                    int i = pidx[k].first;
                    int j = pidx[k].second;
                    int ii = curr + i;
                    int jj = curr + j;
                    GabrielBlock b;
                    b.llim = lc[i][j];
                    b.ulim = uc[i][j];
                    b.frac = vf[k];
                    b.first = ii;
                    b.last = jj;
                    b.length = pos[jj] - pos[ii] + 1;
                    blk.push_back(std::move(b));
                }
            }

            if (curr + w >= n) {
                eprint("INFO: 100% complete\n");
                break;
            }

            int stop = pos[curr + w];
            for (int i = 1; i < w; ++i) {
                int start = pos[++curr];
                if (stop - start + 1 <= par.maxlen)
                    break;
            }

            if (curr + w >= n)
                w = n - curr;

            eprint("INFO: %g%% complete\n", 100.0 * curr / n);
        }

        return 0;
    }

} // namespace

int find_hapblock_gabriel(const std::vector<int> &pos, const std::vector< std::vector<int8_t> > &geno, HapBlockGabriel *out)
{
    isize_t n = length(pos);

    out->block.clear();

    if (n < 2)
        return 0;

    if (n > std::numeric_limits<int>::max()) {
        eprint("ERROR: exceed maximum number (%d) of SNPs in find_hapblock_gabriel: %td\n",
               std::numeric_limits<int>::max(), n);
        return 1;
    }

    int info = 0;
    std::vector<GabrielBlock> blk;

    if (out->par.batch < n) {
        if (out->par.thread > 0)
            info = find_gabriel_split_omp(out->par, pos, geno, blk);
        else
            info = find_gabriel_split(out->par, pos, geno, blk);
    }
    else {
        if (out->par.thread > 0)
            info = find_gabriel_full_omp(out->par, pos, geno, blk);
        else
            info = find_gabriel_full(out->par, pos, geno, blk);
    }

    if (info != 0)
        return 1;

    if (blk.empty())
        return 0;

    auto cmp = [](const GabrielBlock &a, const GabrielBlock &b) {
        if (a.length != b.length)
            return a.length > b.length;
        if (a.frac != b.frac)
            return a.frac > b.frac;
        if (a.llim != b.llim)
            return a.llim > b.llim;
        if (a.ulim != b.ulim)
            return a.ulim > b.ulim;
        return a.first < b.first;
    };

    std::sort(blk.begin(), blk.end(), cmp);

    std::vector<char> in(n + 1, 0);

    for (auto &e : blk) {
        if (in[e.first] || in[e.last])
            continue;
        out->block.emplace_back(pos[e.first], pos[e.last]);
        for (int i = e.first; i <= e.last; ++i)
            in[i] = 1;
    }

    std::sort(out->block.begin(), out->block.end());

    return 0;
}
