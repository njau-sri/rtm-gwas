#include "haplotype.h"

#include <cmath>

#include "print.h"
#include "vectorutil.h"

namespace {

    isize_t calc_identity(const std::vector<allele_t> &x, const std::vector<allele_t> &y)
    {
        isize_t s = 0;

        isize_t n = length(x);
        for (isize_t i = 0; i < n; ++i) {
            if (x[i] && y[i] && x[i] == y[i])
                ++s;
        }

        return s;
    }

    std::string hapstr(const std::vector<std::vector<std::string> > &allele,
                       const std::vector<isize_t> &sidx,
                       const std::vector<allele_t> &hap)
    {
        isize_t n = length(hap);

        std::string s(n, 'N');

        for (isize_t i = 0; i < n; ++i) {
            allele_t a = hap[i];
            if (a) {
                isize_t j = sidx[i];
                s[i] = allele[j][a-1][0];
            }
        }

        return s;
    }

} // namespace

int form_haplotype(const Genotype &gt, const std::vector<isize_t> &sidx, Haplotype *out)
{
    isize_t n = length(gt.ind);
    isize_t m = length(sidx);

    if (m < 2) {
        eprint("ERROR: requires at least two SNPs to form haplotype\n");
        return 1;
    }

    bool diploid = gt.ploidy == 2;
    isize_t nn = diploid ? n*2 : n;

    out->chr = gt.chr[sidx.front()];
    out->start = gt.pos[sidx.front()];
    out->stop = gt.pos[sidx.back()];
    out->size = m;
    out->filter = 0;
    out->dat.assign(nn, 0);
    out->hap.clear();

    // NOTE: assuming phase is known

    std::vector< std::vector<allele_t> > dat(nn, std::vector<allele_t>(m));

    for (isize_t i = 0; i < nn; ++i) {
        for (isize_t j = 0; j < m; ++j) {
            isize_t jj = sidx[j];
            dat[i][j] = gt.dat[jj][i];
        }
    }

    auto hap = unique(dat);

    std::vector<allele_t> mis(m, 0);
    auto itr = std::find(hap.begin(), hap.end(), mis);
    if (itr != hap.end())
        hap.erase(itr);

    if (hap.empty())
        return 0;

    isize_t nh = length(hap);

    if (nh == 1) {
        out->dat.assign(nn, 1);
        out->hap.push_back( hapstr(gt.allele, sidx, hap[0]) );
        return 0;
    }

    const auto& par = out->par;

    isize_t identity = 0;
    if (par.identity > 0.0 && par.identity < 1.0)
        identity = (isize_t) std::ceil(m * par.identity);

    // TODO: improve filtering

    std::vector<isize_t> freq;
    for (auto &e : hap)
        freq.push_back(count(dat, e));

    auto ord = order_desc(freq);
    subset(hap, ord).swap(hap);
    subset(freq, ord).swap(freq);

    isize_t na = nh;
    if (par.maf > 0.0 && par.maf < 1.0) {
        isize_t maf = (isize_t) std::ceil(nn * par.maf);
        na = std::count_if(freq.begin(), freq.end(), [maf](isize_t e) { return e >= maf; });
    }
    if (na == 1)
        na = 2;

    std::vector<int> codec(nh);
    std::iota(codec.begin(), codec.end(), 1);

    for (isize_t i = na; i < nh; ++i) {
        std::vector<isize_t> sc(na);
        for (isize_t j = 0; j < na; ++j)
            sc[j] = calc_identity(hap[j], hap[i]);
        isize_t k = index_max(sc);
        codec[i] = codec[k];
        if (sc[k] < identity)
            codec[i] = 0;
        out->filter += (int) freq[i];
    }

    for (isize_t i = 0; i < nn; ++i) {
        isize_t k = index(hap, dat[i]);
        if (k != -1)
            out->dat[i] = codec[k];
    }

    for (isize_t i = 0; i < na; ++i)
        out->hap.push_back( hapstr(gt.allele, sidx, hap[i]) );

    return 0;
}

int form_haplotype_ril(const Genotype &gt, const std::vector<isize_t> &sidx, Haplotype *out)
{
    isize_t n = length(gt.ind);
    isize_t m = length(sidx);

    if (m < 2) {
        eprint("ERROR: requires at least two SNPs to form RIL haplotype\n");
        return 1;
    }

    if (n < 3) {
        eprint("ERROR: requires at least three individuals to form RIL haplotype\n");
        return 1;
    }

    bool diploid = gt.ploidy == 2;
    isize_t nn = diploid ? (n-2)*2 : n-2;

    out->chr = gt.chr[sidx.front()];
    out->start = gt.pos[sidx.front()];
    out->stop = gt.pos[sidx.back()];
    out->size = m;
    out->filter = 0;
    out->dat.assign(nn, 0);
    out->pdat.assign(diploid ? 4 : 2, 0);
    out->hap.clear();

    // NOTE: assuming parents are homozygous

    std::vector<allele_t> p1, p2;
    p1.reserve(m); p2.reserve(m);

    for (auto j : sidx) {
        p1.push_back(gt.dat[j][0]);
        isize_t i = diploid ? 2 : 1;
        p2.push_back(gt.dat[j][i]);
    }

    if (count(p1, (allele_t) 0) == m || count(p2, (allele_t) 0) == m) {
        auto hs1 = hapstr(gt.allele, sidx, p1);
        auto hs2 = hapstr(gt.allele, sidx, p2);
        eprint("ERROR: invalid parental haplotype in RIL\n");
        eprint("       %s\n", hs1);
        eprint("       %s\n", hs2);
        return 1;
    }

    const auto& par = out->par;

    isize_t identity = 0;
    if (par.identity > 0.0 && par.identity < 1.0)
        identity = (isize_t) std::ceil(m * par.identity);

    std::vector<allele_t> g;
    g.reserve(m);

    isize_t curr = diploid ? 4 : 2;

    for (isize_t i = 0; i < nn; ++i) {
        g.clear();
        for (auto j : sidx)
            g.push_back(gt.dat[j][i+curr]);

        if (g == p1)
            out->dat[i] = 1;
        else if (g == p2)
            out->dat[i] = 2;
        else {
            isize_t s1 = calc_identity(g, p1);
            isize_t s2 = calc_identity(g, p2);
            int a = s1 < identity || s2 < identity ? 0 : (s1 >= s2 ? 1 : 2);
            out->dat[i] = a;
            out->filter += 1;
        }
    }

    if (diploid) {
        out->pdat[0] = out->pdat[1] = 1;
        out->pdat[2] = out->pdat[3] = 2;
    }
    else {
        out->pdat[0] = 1;
        out->pdat[1] = 2;
    }

    out->hap.push_back( hapstr(gt.allele, sidx, p1) );
    out->hap.push_back( hapstr(gt.allele, sidx, p2) );

    return 0;
}

int form_haplotype_nam(const Genotype &gt, const std::vector<isize_t> &sidx, const std::vector<isize_t> &fam,
                       Haplotype *out)
{
    isize_t n = length(gt.ind);
    isize_t f = length(fam);
    isize_t m = length(sidx);

    if (m < 2) {
        eprint("ERROR: requires at least two SNPs to form NAM haplotype\n");
        return 1;
    }

    if (f < 2) {
        eprint("ERROR: requires at least two families to form NAM haplotype\n");
        return 1;
    }

    if (n < f*2 + 1) {
        eprint("ERROR: requires at least %td individuals to form NAM haplotype\n", f*2 + 1);
        return 1;
    }

    isize_t nf = std::accumulate(fam.begin(), fam.end(), (isize_t) 0);

    if (n != nf + f + 1) {
        eprint("ERROR: inconsistent number of lines to form NAM haplotype: %td != %td\n", n, nf + f + 1);
        return 1;
    }

    bool diploid = gt.ploidy == 2;
    isize_t nn = diploid ? nf*2 : nf;

    out->chr = gt.chr[sidx.front()];
    out->start = gt.pos[sidx.front()];
    out->stop = gt.pos[sidx.back()];
    out->size = m;
    out->filter = 0;
    out->dat.assign(nn, 0);
    out->pdat.assign(diploid ? (f+1)*2 : f+1, 0);
    out->hap.clear();

    // NOTE: assuming parents are homozygous

    std::vector< std::vector<allele_t> > pdat(f + 1, std::vector<allele_t>(m));

    for (isize_t i = 0; i <= f; ++i) {
        isize_t ii = diploid ? i*2 : i;
        for (isize_t j = 0; j < m; ++j) {
            isize_t jj = sidx[j];
            pdat[i][j] = gt.dat[jj][ii];
        }
        if (count(pdat[i], (allele_t) 0) == m) {
            auto hs = hapstr(gt.allele, sidx, pdat[i]);
            eprint("ERROR: invalid parental haplotype in NAM\n");
            eprint("       %td: %s\n", i+1, hs);
            return 1;
        }
    }

    const auto& par = out->par;

    isize_t identity = 0;
    if (par.identity > 0.0 && par.identity < 1.0)
        identity = (isize_t) std::ceil(m * par.identity);

    auto hap = stable_unique(pdat);

    std::vector<int> codec;
    for (auto &e : pdat)
        codec.push_back((int) index(hap, e) + 1);

    if (diploid) {
        isize_t lc = length(codec);
        for (isize_t i = 0; i < lc; ++i)
            out->pdat[i*2] = out->pdat[i*2+1] = codec[i];
    }
    else
        out->pdat = codec;

    std::vector<allele_t> g;
    g.reserve(m);

    isize_t curr1 = 0;
    isize_t curr2 = diploid ? (f+1)*2 : f+1;

    for (isize_t k = 0; k < f; ++k) {
        const auto& p1 = pdat[0];
        const auto& p2 = pdat[k+1];

        isize_t nk = diploid ? fam[k] * 2 : fam[k];
        for (isize_t i = 0; i < nk; ++i) {
            isize_t i1 = curr1 + i;
            isize_t i2 = curr2 + i;

            g.clear();
            for (auto j : sidx)
                g.push_back(gt.dat[j][i2]);

            if (g == p1)
                out->dat[i1] = 1;
            else if (g == p2)
                out->dat[i1] = codec[k+1];
            else {
                isize_t s1 = calc_identity(g, p1);
                isize_t s2 = calc_identity(g, p2);
                int a = s1 < identity || s2 < identity ? 0 : (s1 >= s2 ? 1 : 2);
                out->dat[i1] = a != 2 ? a : codec[k+1];
                out->filter += 1;
            }
        }
        curr1 += nk;
        curr2 += nk;
    }

    for (auto &e : hap)
        out->hap.push_back(hapstr(gt.allele, sidx, e));

    return 0;
}
