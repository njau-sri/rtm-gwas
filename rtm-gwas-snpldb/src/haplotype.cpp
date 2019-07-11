#include <cmath>
#include "haplotype.h"
#include "util.h"


using std::size_t;


namespace {

size_t calc_similarity(const std::vector<allele_t> &x, const std::vector<allele_t> &y)
{
    size_t s = 0;

    auto n = x.size();
    for (size_t i = 0; i < n; ++i) {
        if (x[i] == y[i])
            ++s;
    }

    return s;
}

} // namespace


size_t infer_haplotype(double maf, const Genotype &gt, const std::vector<std::size_t> &sidx, Genotype &ht)
{
    bool diploid = gt.ploidy == 2;
    size_t p = diploid ? 2 : 1;

    auto n = gt.ind.size();
    std::vector< std::vector<allele_t> > dat;

    for (size_t i = 0; i < n; ++i) {
        std::vector<allele_t> v1, v2;
        for (auto j : sidx) {
            if ( diploid ) {
                v1.push_back(gt.dat[j][i*2]);
                v2.push_back(gt.dat[j][i*2+1]);
            }
            else
                v1.push_back(gt.dat[j][i]);
        }
        dat.push_back(v1);
        if ( ! v2.empty() )
            dat.push_back(v2);
    }

    auto haps = unique(dat);
    auto nhap = haps.size();

    std::vector<size_t> freq;
    for (auto &e : haps)
        freq.push_back( count(dat,e) );

    auto ord = order(freq, true);
    subset(haps,ord).swap(haps);
    subset(freq,ord).swap(freq);

    auto na = haps.size();
    if (maf > 0.0 && maf < 1.0) {
        auto imaf = static_cast<size_t>( std::ceil(p * n * maf) );
        na = count_if(freq, [imaf](size_t a) { return a >= imaf; });
    }

    // TODO: use cluster analysis if na = 1
    // TODO: maf_new may less than maf_threshold
    if (na == 1)
        na = 2;

    if (na > 255)
        na = 255;

    size_t rec = 0;
    std::vector<allele_t> codec(nhap);
    std::iota(codec.begin(), codec.end(), allele_t(1));

    for (auto i = na; i < nhap; ++i) {
        rec += freq[i];
        std::vector<size_t> v(na);
        for (size_t j = 0; j < na; ++j)
            v[j] = calc_similarity(haps[i], haps[j]);
        auto k = index(v, * std::max_element(v.begin(),v.end()));
        codec[i] = static_cast<allele_t>(k+1);
    }

    if ( diploid )
        rec /= 2;

    std::vector<allele_t> v;

    if ( diploid ) {
        for (size_t i = 0; i < n; ++i) {
            auto k = index(haps, dat[i*2]);
            v.push_back(codec[k]);
            k = index(haps, dat[i*2+1]);
            v.push_back(codec[k]);
        }
    }
    else {
        for (size_t i = 0; i < n; ++i) {
            auto k = index(haps, dat[i]);
            v.push_back(codec[k]);
        }
    }

    if (na == 0)
        std::fill(v.begin(), v.end(), 0);

    ht.dat.push_back(v);

    // first snp position as block position
    int start = std::numeric_limits<int>::max();
    for (auto j : sidx) {
        if (gt.pos[j] < start)
            start = gt.pos[j];
    }

    ht.pos.push_back(start);

    std::vector<std::string> allele;
    for (size_t i = 0; i < na; ++i) {
        std::string si;
        auto p = sidx.size();
        for (size_t j = 0; j < p; ++j) {
            auto jj = sidx[j];
            auto a = haps[i][j];
            if ( a )
                si.append(gt.allele[jj][a-1]);
            else
                si.push_back('N');
        }
        allele.push_back(si);
    }

    ht.allele.push_back(allele);

    return rec;
}

std::size_t infer_haplotype_ril(const Genotype &igt, const Genotype &pgt,
                                const std::vector<std::size_t> &sidx,
                                const std::vector< std::vector<std::size_t> > &fam,
                                Genotype &iht, Genotype &pht)
{
    bool diploid = igt.ploidy == 2;

    // number of parents
    auto np = fam.size() + 1;

    // genotype of parents
    std::vector< std::vector<allele_t> > pdat;
    for (size_t i = 0; i < np; ++i) {
        std::vector<allele_t> v;
        for (auto j : sidx)
            v.push_back(pgt.dat[j][i]);
        pdat.push_back(v);
    }

    // unique parental haplotype
    auto phap = stable_unique(pdat);

    // parental haplotype code
    std::vector<allele_t> codec;
    for (auto &e : pdat)
        codec.push_back( static_cast<allele_t>( index(phap, e) + 1 ) );

    pht.dat.push_back(codec);

    // number of lines
    auto n = igt.ind.size();

    // genotype of lines
    std::vector< std::vector<allele_t> > dat;

    for (size_t i = 0; i < n; ++i) {
        std::vector<allele_t> v1, v2;
        for (auto j : sidx) {
            if ( diploid ) {
                v1.push_back(igt.dat[j][i*2]);
                v2.push_back(igt.dat[j][i*2+1]);
            }
            else
                v1.push_back(igt.dat[j][i]);
        }
        dat.push_back(v1);
        if ( ! v2.empty() )
            dat.push_back(v2);
    }

    size_t rec = 0;
    std::vector<allele_t> v;

    auto& p1hap = pdat[0];
    for (size_t j = 1; j < np; ++j) {
        auto& p2hap = pdat[j];
        for (auto i : fam[j-1]) {
            auto i1 = diploid ? i*2 : i;
            allele_t a1 = 0;
            if (dat[i1] == p1hap)
                a1 = codec[0];
            else if (dat[i1] == p2hap)
                a1 = codec[j];
            else {
                ++rec;
                auto s1 = calc_similarity(dat[i1], p1hap);
                auto s2 = calc_similarity(dat[i1], p2hap);
                a1 = s1 >= s2 ? codec[0] : codec[j];
            }
            v.push_back(a1);

            if ( diploid ) {
                auto i2 = i1 + 1;
                allele_t a2 = 0;
                if (dat[i2] == p1hap)
                    a2 = codec[0];
                else if (dat[i2] == p2hap)
                    a2 = codec[j];
                else {
                    ++rec;
                    auto s1 = calc_similarity(dat[i2], p1hap);
                    auto s2 = calc_similarity(dat[i2], p2hap);
                    a2 = s1 >= s2 ? codec[0] : codec[j];
                }
                v.push_back(a2);
            }
        }
    }

    if ( diploid )
        rec /= 2;

    iht.dat.push_back(v);

    // first snp position as block position
    int start = std::numeric_limits<int>::max();
    for (auto j : sidx) {
        if (igt.pos[j] < start)
            start = igt.pos[j];
    }

    iht.pos.push_back(start);

    auto na = * std::max_element(codec.begin(), codec.end());

    std::vector<std::string> allele;
    for (size_t k = 0; k < na; ++k) {
        auto i = index(codec, k+1);
        std::string si;
        auto p = sidx.size();
        for (size_t j = 0; j < p; ++j) {
            auto a = pdat[i][j];
            auto jj = sidx[j];
            if ( a )
                si.append(igt.allele[jj][a-1]);
            else
                si.push_back('N');
        }
        allele.push_back(si);
    }

    iht.allele.push_back(allele);

    return rec;
}
