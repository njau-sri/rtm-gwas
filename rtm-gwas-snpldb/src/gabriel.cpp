#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>
#include "gabriel.h"


using std::size_t;


namespace {


static const double g_EPS = std::numeric_limits<double>::epsilon();


struct BlockInfo
{
    char llim;
    char ulim;
    short frac;
    int first;
    int last;
    int length;
};


// Gabriel et al., Science, 2002, 296:2225-2229. DOI: 10.1126/science.1069424
// Barrett et al., Bioinformatics, 2005, 21:263-265. DOI: 10.1093/bioinformatics/bth457


// Estimates haplotype frequencies via the EM algorithm
// AB, Ab, aB, ab, AaBb
void calc_hap_prob_EM(int n11, int n12, int n21, int n22, int ndh, double &p11, double &p12, double &p21, double &p22)
{
    static const int maxit = 1000;
    static const double tol = 1e-10;

    double n = n11 + n12 + n21 + n22 + ndh * 2;
    p11 = n11 / n;
    p12 = n12 / n;
    p21 = n21 / n;
    p22 = n22 / n;

    if (ndh == 0)
        return;

    auto cp11 = p11;
    auto cp12 = p12;
    auto cp21 = p21;
    auto cp22 = p22;

    auto h = ndh / n;
    auto x = h / 2;
    auto y = h - x;

    for (int i = 0; i < maxit; ++i) {
        p11 = cp11 + x;
        p12 = cp12 + y;
        p21 = cp21 + y;
        p22 = cp22 + x;
        auto z = h * p11 * p22 / (p11 * p22 + p12 * p21);
        if (std::fabs(x - z) < tol)
            break;
        x = z;
        y = h - x;
    }
}

// 0,1,2,3 -> AA,Aa,aa,NN
void count_freq_3x3(const std::vector<char> &x, const std::vector<char> &y, int freq[3][3])
{
    auto n = x.size();
    for (size_t i = 0; i < n; ++i)
        if (x[i] < 3 && y[i] < 3)
            ++freq[x[i]][y[i]];
}

// D' 95% confidence interval estimate
int calc_dprime(int freq[3][3], int &lower, int &upper)
{
    using std::log;

    int n11 = freq[0][0] * 2 + freq[0][1] + freq[1][0];
    int n12 = freq[0][2] * 2 + freq[0][1] + freq[1][2];
    int n21 = freq[2][0] * 2 + freq[1][0] + freq[2][1];
    int n22 = freq[2][2] * 2 + freq[2][1] + freq[1][2];
    int ndh = freq[1][1];

    double nn = n11 + n12 + n21 + n22 + ndh * 2;
    if (nn < 4)
        return 1;

    double p11 = 0.0;
    double p12 = 0.0;
    double p21 = 0.0;
    double p22 = 0.0;

    calc_hap_prob_EM(n11, n12, n21, n22, ndh, p11, p12, p21, p22);

    if (n11 > n22) {
        std::swap(n11, n22);
        std::swap(n12, n21);
        std::swap(p11, p22);
        std::swap(p12, p21);
    }

    if (n11 > n12 || n11 > n21) {
        if (n12 < n21) {
            std::swap(n11, n12);
            std::swap(n21, n22);
            std::swap(p11, p12);
            std::swap(p21, p22);
        }
        else {
            std::swap(n11, n21);
            std::swap(n12, n22);
            std::swap(p11, p21);
            std::swap(p12, p22);
        }
    }

    auto p1x = (n11 + n12 + ndh) / nn;
    auto p2x = 1 - p1x;
    auto px1 = (n11 + n21 + ndh) / nn;
    auto px2 = 1 - px1;

    auto D = p11 - p1x*px1;
    auto Dmax = D < 0.0 ? std::min(p1x*px1, p2x*px2) : std::min(p1x*px2, p2x*px1);

    if (p11 < 1e-10) p11 = 1e-10;
    if (p12 < 1e-10) p12 = 1e-10;
    if (p21 < 1e-10) p21 = 1e-10;
    if (p22 < 1e-10) p22 = 1e-10;
    auto LL1 = n11*log(p11) + n12*log(p12) + n21*log(p21) + n22*log(p22) + ndh*log(p11*p22 + p12*p21);

    if (D < 0.0) {
        std::swap(p1x, p2x);
        std::swap(n11, n21);
        std::swap(n12, n22);
    }

    double tp = 0.0;
    double ls[101];
    auto Dstep = Dmax / 100;

    for (int i = 0; i <= 100; ++i) {
        auto q11 = i*Dstep + p1x*px1;
        auto q12 = p1x - q11;
        auto q21 = px1 - q11;
        auto q22 = p2x - q21;
        if (i == 100) {
            if (q11 < 1e-10) q11 = 1e-10;
            if (q12 < 1e-10) q12 = 1e-10;
            if (q21 < 1e-10) q21 = 1e-10;
            if (q22 < 1e-10) q22 = 1e-10;
        }
        auto LL2 = n11*log(q11) + n12*log(q12) + n21*log(q21) + n22*log(q22) + ndh*log(q11*q22 + q12*q21);
        auto prob = std::exp(LL2 - LL1);
        ls[i] = prob;
        tp += prob;
    }

    double sp = 0.0;
    auto tp5 = tp * 0.05;
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

double calc_frac_info_pair(size_t x, size_t y, const std::vector< std::vector<char> > &sr,
                           const std::vector< std::vector<char> > &ss)
{
    int s = 0;
    int r = 0;
    auto n = y - x + 1;

    for (size_t i = x; i <= y; ++i) {
        for (auto j = i + 1; j <= y; ++j) {
            r += sr[j][i];
            if (n > 4)
                s += sr[i][j];
            else if (n == 2)
                s += ss[i][j];
            else if (n == 3 || n == 4)
                s += ss[j][i];
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

    return static_cast<double>(s) / t;
}

int find_f(const BlockGabriel &par, const std::vector<int> &pos, const std::vector< std::vector<char> > &dat,
           std::vector<BlockInfo> &blk)
{
    auto m = pos.size();

    std::vector< std::vector<char> > ci(m, std::vector<char>(m, -1));
    std::vector< std::vector<char> > sr(m, std::vector<char>(m, 0));
    std::vector< std::vector<char> > ss(m, std::vector<char>(m, 0));

    for (size_t i = 0; i < m; ++i) {
        for (auto j = i + 1; j < m; ++j) {
            if (pos[j] <= pos[i]) {
                std::cerr << "ERROR: position must be in ascending order: " << pos[i] << " " << pos[j] << "\n";
                return 1;
            }
            if (pos[j] - pos[i] > par.maxlen)
                break;

            int freq[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
            count_freq_3x3(dat[i], dat[j], freq);

            int llim = -1, ulim = -1;
            calc_dprime(freq, llim, ulim);

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

            ci[j][i] = static_cast<char>(llim);
            ci[i][j] = static_cast<char>(ulim);

            sr[i][j] = strong;
            sr[j][i] = recomb;

            ss[i][j] = strong2;
            ss[j][i] = strong34;
        }
    }

    for (size_t i = 0; i < m; ++i) {
        for (auto j = i + 1; j < m; ++j) {
            auto dist = pos[j] - pos[i];
            if (dist > par.maxlen)
                break;
            if (ci[j][i] < par.llim || ci[i][j] < par.ulim)
                continue;
            auto q = j - i + 1;
            if ((q == 2 && dist > 20000) || (q == 3 && dist > 30000))
                continue;
            auto frac = calc_frac_info_pair(i, j, sr, ss);
            if (frac - par.frac > g_EPS) {
                BlockInfo bi;
                bi.llim = ci[j][i];
                bi.ulim = ci[i][j];
                bi.frac = static_cast<short>(frac*10000);
                bi.first = static_cast<int>(i);
                bi.last = static_cast<int>(j);
                bi.length = dist;
                blk.push_back(bi);
            }
        }
    }

    return 0;
}

int find_w(const BlockGabriel &par, const std::vector<int> &pos, const std::vector< std::vector<char> > &dat,
           std::vector<BlockInfo> &blk)
{
    auto m = pos.size();
    auto w = static_cast<size_t>(par.batch);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            if (pos[j] <= pos[i]) {
                std::cerr << "ERROR: position must be in ascending order: " << pos[i] << " " << pos[j] << "\n";
                return 1;
            }
            if (pos[j] - pos[i] > par.maxlen)
                break;
            auto wij = j - i + 1;
            if (w < wij)
                w = wij;
        }
    }

    std::vector< std::vector<char> > ci(w, std::vector<char>(w, -1));
    std::vector< std::vector<char> > sr(w, std::vector<char>(w, 0));
    std::vector< std::vector<char> > ss(w, std::vector<char>(w, 0));

    size_t curr = 0;

    for (;;) {
        for (size_t i = 0; i < w; ++i) {
            auto ii = curr + i;
            for (auto j = i + 1; j < w; ++j) {
                auto jj = curr + j;
                if (pos[jj] - pos[ii] > par.maxlen)
                    break;

                int freq[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
                count_freq_3x3(dat[ii], dat[jj], freq);

                int llim = -1, ulim = -1;
                calc_dprime(freq, llim, ulim);

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

                ci[j][i] = static_cast<char>(llim);
                ci[i][j] = static_cast<char>(ulim);

                sr[i][j] = strong;
                sr[j][i] = recomb;

                ss[i][j] = strong2;
                ss[j][i] = strong34;
            }
        }

        for (size_t i = 0; i < w; ++i) {
            auto ii = curr + i;
            for (auto j = i + 1; j < w; ++j) {
                auto jj = curr + j;
                auto dist = pos[jj] - pos[ii];
                if (dist > par.maxlen)
                    break;
                if (ci[j][i] < par.llim || ci[i][j] < par.ulim)
                    continue;
                auto q = j - i + 1;
                if ((q == 2 && dist > 20000) || (q == 3 && dist > 30000))
                    continue;
                auto frac = calc_frac_info_pair(i, j, sr, ss);
                if (frac - par.frac > g_EPS) {
                    BlockInfo bi;
                    bi.llim = ci[j][i];
                    bi.ulim = ci[i][j];
                    bi.frac = static_cast<short>(frac*10000);
                    bi.first = static_cast<int>(ii);
                    bi.last = static_cast<int>(jj);
                    bi.length = dist;
                    blk.push_back(bi);
                }
            }
        }

        if (curr + w >= m) {
            std::cerr << "INFO: 100% complete\n";
            break;
        }

        auto pos2 = pos[curr + w];
        for (size_t i = 1; i < w; ++i) {
            auto pos1 = pos[++curr];
            auto dist = pos2 - pos1;
            if (dist <= par.maxlen)
                break;
        }

        if (curr + w >= m)
            w = m - curr;

        std::cerr << "INFO: " << 100.0 * curr / m << "% complete\n";
    }

    return 0;
}

int find_f_omp(const BlockGabriel &par, const std::vector<int> &pos, const std::vector< std::vector<char> > &dat,
               std::vector<BlockInfo> &blk)
{
    auto m = pos.size();

    std::vector< std::pair<size_t,size_t> > pidx;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            if (pos[j] <= pos[i]) {
                std::cerr << "ERROR: position must be in ascending order: " << pos[i] << " " << pos[j] << "\n";
                return 1;
            }
            if (pos[j] - pos[i] > par.maxlen)
                break;
            pidx.emplace_back(i,j);
        }
    }

    std::vector< std::vector<char> > ci(m, std::vector<char>(m, -1));
    std::vector< std::vector<char> > sr(m, std::vector<char>(m, 0));
    std::vector< std::vector<char> > ss(m, std::vector<char>(m, 0));

    auto n = pidx.size();

    #pragma omp parallel for
    for (size_t k = 0; k < n; ++k) {
        auto i = pidx[k].first;
        auto j = pidx[k].second;

        int freq[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
        count_freq_3x3(dat[i], dat[j], freq);

        int llim = -1, ulim = -1;
        calc_dprime(freq, llim, ulim);

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

        ci[j][i] = static_cast<char>(llim);
        ci[i][j] = static_cast<char>(ulim);

        sr[i][j] = strong;
        sr[j][i] = recomb;

        ss[i][j] = strong2;
        ss[j][i] = strong34;
    }

    std::vector<short> vfrac(n,-1);

    #pragma omp parallel for schedule(dynamic)
    for (size_t k = 0; k < n; ++k) {
        auto i = pidx[k].first;
        auto j = pidx[k].second;
        auto dist = pos[j] - pos[i];

        if (ci[j][i] < par.llim || ci[i][j] < par.ulim)
            continue;

        auto q = j - i + 1;
        if ((q == 2 && dist > 20000) || (q == 3 && dist > 30000))
            continue;

        auto frac = calc_frac_info_pair(i, j, sr, ss);
        if (frac - par.frac > g_EPS)
            vfrac[k] = static_cast<short>(frac*10000);
    }

    size_t nb = 0;
    for (auto e : vfrac) {
        if (e > 0)
            ++nb;
    }

    blk.reserve(nb);

    for (size_t k = 0; k < n; ++k) {
        auto i = pidx[k].first;
        auto j = pidx[k].second;
        auto dist = pos[j] - pos[i];
        if (vfrac[k] > 0) {
            BlockInfo bi;
            bi.llim = ci[j][i];
            bi.ulim = ci[i][j];
            bi.frac = vfrac[k];
            bi.first = static_cast<int>(i);
            bi.last = static_cast<int>(j);
            bi.length = dist;
            blk.push_back(bi);
        }
    }

    return 0;
}

int find_w_omp(const BlockGabriel &par, const std::vector<int> &pos, const std::vector< std::vector<char> > &dat,
               std::vector<BlockInfo> &blk)
{
    auto m = pos.size();
    auto w = static_cast<size_t>(par.batch);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            if (pos[j] <= pos[i]) {
                std::cerr << "ERROR: position must be in ascending order: " << pos[i] << " " << pos[j] << "\n";
                return 1;
            }
            if (pos[j] - pos[i] > par.maxlen)
                break;
            auto wij = j - i + 1;
            if (w < wij)
                w = wij;
        }
    }

    std::vector<short> vfrac;
    std::vector< std::pair<size_t,size_t> > pidx;
    std::vector< std::vector<char> > ci(w, std::vector<char>(w, 0));
    std::vector< std::vector<char> > sr(w, std::vector<char>(w, 0));
    std::vector< std::vector<char> > ss(w, std::vector<char>(w, 0));

    size_t curr = 0;

    for (;;) {
        pidx.clear();
        for (size_t i = 0; i < w; ++i) {
            auto ii = curr + i;
            for (size_t j = i + 1; j < w; ++j) {
                auto jj = curr + j;
                if (pos[jj] - pos[ii] > par.maxlen)
                    break;
                pidx.emplace_back(i,j);
            }
        }

        auto n = pidx.size();

        #pragma omp parallel for
        for (size_t k = 0; k < n; ++k) {
            auto i = pidx[k].first;
            auto j = pidx[k].second;

            auto ii = curr + i;
            auto jj = curr + j;

            int freq[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
            count_freq_3x3(dat[ii], dat[jj], freq);

            int llim = -1, ulim = -1;
            calc_dprime(freq, llim, ulim);

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

            ci[j][i] = static_cast<char>(llim);
            ci[i][j] = static_cast<char>(ulim);

            sr[i][j] = strong;
            sr[j][i] = recomb;

            ss[i][j] = strong2;
            ss[j][i] = strong34;
        }

        vfrac.assign(n, -1);

        #pragma omp parallel for schedule(dynamic)
        for (size_t k = 0; k < n; ++k) {
            auto i = pidx[k].first;
            auto j = pidx[k].second;

            if (ci[j][i] < par.llim || ci[i][j] < par.ulim)
                continue;

            auto ii = curr + i;
            auto jj = curr + j;
            auto dist = pos[jj] - pos[ii];

            auto q = j - i + 1;
            if ((q == 2 && dist > 20000) || (q == 3 && dist > 30000))
                continue;

            auto frac = calc_frac_info_pair(i, j, sr, ss);
            if (frac - par.frac > g_EPS)
                vfrac[k] = static_cast<short>(frac*10000);
        }

        size_t nb = 0;
        for (auto e : vfrac) {
            if (e > 0)
                ++nb;
        }

        if (blk.capacity() < (blk.size() + nb))
            blk.reserve(blk.size() + nb);

        for (size_t k = 0; k < n; ++k) {
            auto i = pidx[k].first;
            auto j = pidx[k].second;
            auto ii = curr + i;
            auto jj = curr + j;
            auto dist = pos[jj] - pos[ii];
            if (vfrac[k] > 0) {
                BlockInfo bi;
                bi.llim = ci[j][i];
                bi.ulim = ci[i][j];
                bi.frac = vfrac[k];
                bi.first = static_cast<int>(ii);
                bi.last = static_cast<int>(jj);
                bi.length = dist;
                blk.push_back(bi);
            }
        }

        if (curr + w >= m) {
            std::cerr << "INFO: 100% complete\n";
            break;
        }

        auto pos2 = pos[curr + w];
        for (size_t i = 1; i < w; ++i) {
            auto pos1 = pos[++curr];
            auto dist = pos2 - pos1;
            if (dist <= par.maxlen)
                break;
        }

        if (curr + w >= m)
            w = m - curr;

        std::cerr << "INFO: " << 100.0 * curr / m << "% complete\n";
    }

    return 0;
}


} // namespace


int find_block(const BlockGabriel &par, const std::vector<int> &pos, const std::vector< std::vector<char> > &dat,
               std::vector< std::pair<int,int> > &bpos)
{
    auto m = pos.size();
    if (m < 2)
        return 0;

    std::vector<BlockInfo> blk;

    int ret = static_cast<size_t>(par.batch) < m ?
                find_w(par, pos, dat, blk) : find_f(par, pos, dat, blk);
    if (ret != 0)
        return 1;

    if ( blk.empty() )
        return 0;

    auto cmp = [](const BlockInfo &a, const BlockInfo &b) {
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

    std::vector<char> in(m+1,0);

    for (auto &e : blk) {
        auto i = static_cast<size_t>(e.first);
        auto j = static_cast<size_t>(e.last);
        if (in[i] || in[j])
            continue;

        bpos.emplace_back(pos[i], pos[j]);

        for (auto k = i; k <= j; ++k)
            in[k] = 1;
    }

    std::sort(bpos.begin(), bpos.end());

    return 0;
}


int find_block_omp(const BlockGabriel &par, const std::vector<int> &pos, const std::vector< std::vector<char> > &dat,
                   std::vector< std::pair<int,int> > &bpos)
{
    auto m = pos.size();
    if (m < 2)
        return 0;

    std::vector<BlockInfo> blk;

    int ret = static_cast<size_t>(par.batch) < m ?
                find_w_omp(par, pos, dat, blk) : find_f_omp(par, pos, dat, blk);
    if (ret != 0)
        return 1;

    if ( blk.empty() )
        return 0;

    auto cmp = [](const BlockInfo &a, const BlockInfo &b) {
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

    std::vector<char> in(m+1,0);

    for (auto &e : blk) {
        auto i = static_cast<size_t>(e.first);
        auto j = static_cast<size_t>(e.last);
        if (in[i] || in[j])
            continue;

        bpos.emplace_back(pos[i], pos[j]);

        for (auto k = i; k <= j; ++k)
            in[k] = 1;
    }

    std::sort(bpos.begin(), bpos.end());

    return 0;
}
