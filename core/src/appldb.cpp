#include <cmath>
#include <map>
#include <limits>
#include <memory>
#include <fstream>
#include <iostream>
#include "appldb.h"
#include "cmdline.h"
#include "util.h"
#include "strsplit.h"
#include "haploprob.h"
#include "vcfio.h"
#include "plinkio.h"
#include "hapmapio.h"
#include "rtmio.h"

namespace {

struct Block
{
    int first;
    int last;
    int length;
    float frac;
};

size_t match(const vector<allele_t> &x, const vector<allele_t> &y)
{
    auto n = x.size();
    size_t c = 0;
    for (size_t i = 0; i < n; ++i) {
        if (x[i] == y[i])
            ++c;
    }
    return c;
}

int read_block(const string &filename, vector<string> &chr, vector<int> &pos)
{
    std::ifstream ifs(filename);

    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file: " << filename << "\n";
        return 1;
    }

    auto delim = [](char c) { return c == ' ' || c == '\t' || c == '\r'; };

    size_t ln = 0;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<string> vs;
        strsplit(delim, line.begin(), line.end(), vs);
        if (vs.empty())
            continue;

        if (vs.size() < 3) {
            std::cerr << "ERROR: expected at least 3 columns at line " << ln << ": " << filename << "\n";
            return 1;
        }

        auto start = number<int>(vs[1]);
        auto stop = number<int>(vs[2]);

        if (stop <= start) {
            std::cerr << "ERROR: incorrect block definition at line " << ln << ":"
                      << vs[0] << " " << vs[1] << " " << vs[2] << "\n";
            return 1;
        }

        chr.push_back(vs[0]);
        pos.push_back(start);
        pos.push_back(stop);
    }

    return 0;
}

} // namespace

int AppLDB::run(int argc, char *argv[])
{
    auto cmd = std::make_shared<CmdLine>("ldb [options]");

    cmd->add("--vcf", "VCF file", "");
    cmd->add("--ped", "PLINK PED/MAP file prefix", "");
    cmd->add("--hmp", "HapMap file", "");
    cmd->add("--geno", "genotype file", "");
    cmd->add("--block", "predefined block file", "");
    cmd->add("--out", "output file", "appldb");
    cmd->add("--mhf", "minimum minor haplotype frequency", "0.01");
    cmd->add("--maxlen", "maximum length of blocks", "200000");

    cmd->add("--low", "lower CI for strong LD", "0.7");
    cmd->add("--high", "upper CI for string LD", "0.98");
    cmd->add("--rec", "upper CI for strong recombination", "0.9");
    cmd->add("--frac", "minimum fraction of informative strong LD", "0.95");

    if (argc < 2) {
        cmd->help();
        return 1;
    }

    cmd->parse(argc, argv);

    m_par.vcf = cmd->get("--vcf");
    m_par.ped = cmd->get("--ped");
    m_par.hmp = cmd->get("--hmp");
    m_par.geno = cmd->get("--geno");
    m_par.block = cmd->get("--block");
    m_par.out = cmd->get("--out");

    m_par.mhf = number<double>(cmd->get("--mhf"));
    m_par.maxlen = number<int>(cmd->get("--maxlen"));

    m_par.low = std::round(100*number<double>(cmd->get("--low")));
    m_par.high = std::round(100*number<double>(cmd->get("--high")));
    m_par.rec = std::round(100*number<double>(cmd->get("--rec")));
    m_par.frac = number<double>(cmd->get("--frac"));

    int info = perform();

    return info;
}

int AppLDB::perform()
{
    load_genotype();

    if ( ! m_par.block.empty() )
        load_block();

    auto chrid = stable_unique(m_gt.chr);
    auto nchr = chrid.size();
    auto indx = index_locus(chrid);

    if ( m_chrid.empty() ) {
        for (size_t i = 0; i < nchr; ++i) {
            vector<int> start, stop;
            find_block(indx[i], start, stop);
            m_chrid.insert(m_chrid.end(), start.size(), chrid[i]);
            m_start.insert(m_start.end(), start.begin(), start.end());
            m_stop.insert(m_stop.end(), stop.begin(), stop.end());
        }
    }

    auto m = m_gt.loc.size();
    auto nb = m_chrid.size();

    auto gt = std::make_shared<Genotype>();
    vector<bool> inblock(m, false);

    for (size_t i = 0; i < nchr; ++i) {
        auto ldb = std::make_shared<Genotype>();
        for (size_t k = 0; k < nb; ++k) {
            if (m_chrid[k] != chrid[i])
                continue;
            vector<size_t> snps;
            for (auto j : indx[i]) {
                if (m_gt.pos[j] < m_start[k] || m_gt.pos[j] > m_stop[k])
                    continue;
                inblock[j] = true;
                snps.push_back(j);
            }
            if ( snps.empty() ) {
                std::cerr << "WARNING: there is no SNPs within block: "
                          << chrid[i] << " " << m_start[k] << " " << m_stop[k] << "\n";
                continue;
            }
            make_snpldb(snps, *ldb);
        }

        for (auto j : indx[i]) {
            if ( inblock[j] ) {
                auto found = std::find(ldb->pos.begin(), ldb->pos.end(), m_gt.pos[j]);
                if (found != ldb->pos.end()) {
                    auto wh = found - ldb->pos.begin();
                    gt->loc.push_back(ldb->loc[wh]);
                    gt->chr.push_back(ldb->chr[wh]);
                    gt->pos.push_back(ldb->pos[wh]);
                    gt->dist.push_back(ldb->dist[wh]);
                    gt->dat.push_back(ldb->dat[wh]);
                    gt->allele.push_back(ldb->allele[wh]);
                }
            }
            else {
                gt->loc.push_back(m_gt.loc[j]);
                gt->chr.push_back(m_gt.chr[j]);
                gt->pos.push_back(m_gt.pos[j]);
                gt->dist.push_back(m_gt.dist[j]);
                gt->dat.push_back(m_gt.dat[j]);
                gt->allele.push_back(m_gt.allele[j]);
            }
        }
    }

    gt->ind = m_gt.ind;
    gt->ploidy = m_gt.ploidy;

    int info = write_vcf(*gt, m_par.out + ".vcf");
    if (info != 0)
        return 1;

    return 0;
}

void AppLDB::load_genotype()
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

void AppLDB::load_block()
{
    std::cerr << "INFO: reading block file...\n";

    vector<int> pos;
    vector<string> chr;
    int ret = read_block(m_par.block, chr, pos);

    if (ret == 0) {
        auto n = chr.size();
        m_chrid.swap(chr);
        for (size_t i = 0; i < n; ++i) {
            m_start.push_back(pos[i*2]);
            m_stop.push_back(pos[i*2+1]);
        }
    }

    std::cerr << "INFO: " << m_chrid.size() << " blocks were observed\n";
}

void AppLDB::make_snpldb(const vector<size_t> &snps, Genotype &gt) const
{
    int n = m_gt.ind.size();
    bool haploid = m_gt.ploidy == 1;

    vector< vector<allele_t> > dat;
    for (int i = 0; i < n; ++i) {
        vector<allele_t> v1, v2;
        for (auto j : snps) {
            if ( haploid ) {
                v1.push_back(m_gt.dat[j][i]);
            }
            else {
                v1.push_back(m_gt.dat[j][i*2]);
                v2.push_back(m_gt.dat[j][i*2+1]);
            }
        }
        dat.push_back(v1);
        if ( ! v2.empty() )
            dat.push_back(v2);
    }

    auto haps = unique(dat);
    int m = haps.size();

    vector<int> freq;
    for (auto &e : haps) {
        auto c = std::count(dat.begin(), dat.end(), e);
        freq.push_back(c);
    }

    auto ord = order(freq);
    std::reverse(ord.begin(), ord.end());
    subset(haps,ord).swap(haps);
    subset(freq,ord).swap(freq);

    auto mhc = static_cast<int>(std::ceil(m_par.mhf*n));
    int na = std::count_if(freq.begin(), freq.end(), [mhc](int a) { return a >= mhc; });

    vector<int> hc(m);
    std::iota(hc.begin(), hc.end(), 0);
    for (int i = na; i < m; ++i) {
        vector<int> ns;
        for (int j = 0; j < na; ++j) {
            int c = match(haps[i], haps[j]);
            ns.push_back(c);
        }
        int wh = index(ns, *std::max_element(ns.begin(), ns.end()));
        hc[i] = wh;
    }

    auto start = m_gt.pos[snps[0]];
    auto stop = start;
    for (auto j : snps) {
        auto pos = m_gt.pos[j];
        if (pos < start)
            start = pos;
        else if (pos > stop)
            stop = pos;
    }

    auto loc = m_gt.chr[snps[0]] + "_LDB_" + std::to_string(start) + "_" + std::to_string(stop);
    gt.loc.push_back(loc);
    gt.chr.push_back(m_gt.chr[snps[0]]);
    gt.pos.push_back(start);
    gt.dist.push_back(0.0);

    vector<allele_t> v;
    if ( haploid ) {
        for (int i = 0; i < n; ++i) {
            auto wh = index(haps, dat[i]);
            auto code = static_cast<allele_t>(hc[wh] + 1);
            v.push_back(code);
        }
    }
    else {
        for (int i = 0; i < n; ++i) {
            auto wh = index(haps, dat[i*2]);
            auto code = static_cast<allele_t>(hc[wh] + 1);
            v.push_back(code);
            wh = index(haps, dat[i*2+1]);
            code = static_cast<allele_t>(hc[wh] + 1);
            v.push_back(code);
        }
    }
    gt.dat.push_back(v);

    vector<string> allele;
    for (int i = 0; i < na; ++i) {
        string s;
        for (size_t j = 0; j < snps.size(); ++j)
            s.append(m_gt.allele[snps[j]][haps[i][j]-1]);
        allele.push_back(s);
    }
    gt.allele.push_back(allele);
}

vector< vector<size_t> > AppLDB::index_locus(const vector<string> &chrid) const
{
    int nchr = chrid.size();

    std::map<string,int> mapchr;
    for (int i = 0; i < nchr; ++i)
        mapchr[chrid[i]] = i;

    auto m = m_gt.loc.size();

    vector< vector<size_t> > indx(nchr);
    for (size_t i = 0; i < m; ++i) {
        auto j = mapchr[ m_gt.chr[i] ];
        indx[j].push_back(i);
    }

    return indx;
}

void AppLDB::find_block(const vector<size_t> &loc, vector<int> &start, vector<int> &stop) const
{
    int n = loc.size();

    int w = 0;
    for (int i = 0; i < n; ++i) {
        auto pi = m_gt.pos[loc[i]];
        for (int j = i + 1; j < n; ++j) {
            auto pj = m_gt.pos[loc[j]];
            if (std::abs(pj-pi) <= m_par.maxlen && w < (j-i+1))
                w = j-i+1;
        }
    }

    vector<Block> vb;
    vector< vector<int> > cild(w,vector<int>(w,-1));

    for (int i = 1; i < w; ++i) {
        auto li = loc[i-1];
        auto pi = m_gt.pos[li];
        for (int j = i + 1; j < w; ++j) {
            auto lj = loc[j-1];
            auto pj = m_gt.pos[lj];
            if (std::abs(pj-pi) <= m_par.maxlen)
                calc_cild(li, lj, cild[j][i], cild[i][j]);
        }
    }

    for (int i = 1; i < w; ++i) {
        auto li = loc[i-1];
        auto pi = m_gt.pos[li];
        for (int j = i + 1; j < w; ++j) {
            auto lj = loc[j-1];
            auto pj = m_gt.pos[lj];
            auto len = std::abs(pj-pi);
            float frac = test_block(i, j, len, cild);
            if (frac > 0.0)
                vb.emplace_back( Block{i-1, j-1, len, frac} );
        }
    }

    for (int k = 0; k <= n-w; ++k) {
        for (int i = 1; i < w; ++i)
            for (int j = 1; j < w; ++j)
                cild[i-1][j-1] = cild[i][j];
        for (int i = 0, j = w-1; i < w-1; ++i) {
            auto li = loc[k+i], lj = loc[k+j];
            auto pi = m_gt.pos[li], pj = m_gt.pos[lj];
            cild[i][j] = cild[j][i] = -1;
            if (std::abs(pj-pi) <= m_par.maxlen)
                calc_cild(li, lj, cild[j][i], cild[i][j]);
        }
        for (int i = 0, j = w-1; i < w-1; ++i) {
            auto li = loc[k+i], lj = loc[k+j];
            auto pi = m_gt.pos[li], pj = m_gt.pos[lj];
            auto len = std::abs(pj-pi);
            float frac = test_block(i, j, len, cild);
            if (frac > 0)
                vb.emplace_back( Block{k+i, k+j, len, frac} );
        }
    }

    auto comp = [&](const Block &a, const Block &b) {
        if (a.length > b.length) return true;
        if (a.length < b.length) return false;
        if ((b.first > a.first && b.first < a.last) || (b.last > a.first && b.last < a.last)) {
            if (a.frac > b.frac) return true;
            if (a.frac < b.frac) return false;
            int a1 = -1, a2 = -1, b1 = -1, b2 = -1;
            calc_cild(loc[a.first], loc[a.last], a1, a2);
            calc_cild(loc[b.first], loc[b.last], b1, b2);
            if (a1 > b1) return true;
            if (a1 < b1) return false;
        }
        return a.first < b.first;
    };

    std::sort(vb.begin(), vb.end(), comp);

    vector<bool> inblock(n + 1, false);

    for (auto &e : vb) {
        auto first = e.first, last = e.last;
        if (inblock[first] || inblock[last])
            continue;
        start.push_back(m_gt.pos[loc[first]]);
        stop.push_back(m_gt.pos[loc[last]]);
        for (auto i = first; i <= last; ++i)
            inblock[i] = true;
    }

    auto ord = order(start);
    subset(start,ord).swap(start);
    subset(stop,ord).swap(stop);
}

// Modified from Haploview-v4.2: edu.mit.wi.haploview/HaploData.java
// https://github.com/jazzywhit/Haploview
int AppLDB::calc_cild(size_t l1, size_t l2, int &low, int &high) const
{
    using std::log;

    if (m_gt.allele[l1].size() != 2 || m_gt.allele[l2].size() != 2) {
        std::cerr << "WARNING: not polymorphic or bi-allelic marker: " << l1 << "/" << l2 << "\n";
        return -1;
    }

    static const double eps = std::numeric_limits<double>::epsilon();

    auto n = m_gt.ind.size();
    size_t nv = 0;

    auto& g1 = m_gt.dat[l1];
    auto& g2 = m_gt.dat[l2];
    int freq[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };

    if (g1.size() == n*2) {
        for (size_t i = 0; i < n; ++i) {
            auto i1 = i*2, i2 = i*2 + 1;
            if ( g1[i1] && g1[i2] && g2[i1] && g2[i2] ) {
                auto a = g1[i1] + g1[i2] - 2;
                auto b = g2[i1] + g2[i2] - 2;
                freq[a][b] += 1;
                ++nv;
            }
        }
    }
    else {
        for (size_t i = 0; i < n; ++i) {
            if ( g1[i] && g2[i] ) {
                auto a = g1[i]*2 - 2;
                auto b = g2[i]*2 - 2;
                freq[a][b] += 1;
                ++nv;
            }
        }
    }

    if (nv == 0)
        return -1;

    HaploProb hp;
    hp.n_AABB = freq[0][0]; hp.n_AABb = freq[0][1]; hp.n_AAbb = freq[0][2];
    hp.n_AaBB = freq[1][0]; hp.n_AaBb = freq[1][1]; hp.n_Aabb = freq[1][2];
    hp.n_aaBB = freq[2][0]; hp.n_aaBb = freq[2][1]; hp.n_aabb = freq[2][2];
    hp.solve();

    auto nAB = hp.n_AABB*2 + hp.n_AABb + hp.n_AaBB;
    auto nAb = hp.n_AAbb*2 + hp.n_AABb + hp.n_Aabb;
    auto naB = hp.n_aaBB*2 + hp.n_AaBB + hp.n_aaBb;
    auto nab = hp.n_aabb*2 + hp.n_aaBb + hp.n_Aabb;

    auto pA = (double) (nAB + nAb + hp.n_AaBb) / (nv*2);
    auto pB = (double) (nAB + naB + hp.n_AaBb) / (nv*2);
    auto pa = 1.0 - pA;
    auto pb = 1.0 - pB;

    auto D = hp.p_AB - pA*pB;
    auto Dmax = D > 0.0 ? std::min(pA*pb, pa*pB) : std::min(pA*pB, pa*pb);

    if (hp.p_AB < eps) hp.p_AB = eps;
    if (hp.p_Ab < eps) hp.p_Ab = eps;
    if (hp.p_aB < eps) hp.p_aB = eps;
    if (hp.p_ab < eps) hp.p_ab = eps;

    auto LL1 = nAB*log(hp.p_AB) + nAb*log(hp.p_Ab) + naB*log(hp.p_aB) + nab*log(hp.p_ab) +
            hp.n_AaBb*log(hp.p_AB*hp.p_ab + hp.p_Ab*hp.p_aB);

    if (D < 0.0) {
        std::swap(pA, pa);
        std::swap(nAB, naB);
        std::swap(nAb, nab);
    }

    double tp = 0.0;
    vector<double> ls(101);

    for (int i = 0; i < 101; ++i) {
        auto p_AB = 0.01*i*Dmax + pA*pB;
        auto p_Ab = pA - p_AB;
        auto p_aB = pB - p_AB;
        auto p_ab = pa - p_aB;
        if (i == 100) {
            if (p_AB < eps) p_AB = eps;
            if (p_Ab < eps) p_Ab = eps;
            if (p_aB < eps) p_aB = eps;
            if (p_ab < eps) p_ab = eps;
        }
        auto LL2 = nAB*log(p_AB) + nAb*log(p_Ab) + naB*log(p_aB) + nab*log(p_ab) +
                hp.n_AaBb*log(p_AB*p_ab + p_Ab*p_aB);
        ls[i] = std::exp(LL2 - LL1);
        tp += ls[i];
    }

    auto tp5 = tp * 0.05;
    double sp = 0.0;

    for (int i = 0; i < 101; ++i) {
        sp += ls[i];
        if (sp > tp5 && sp - ls[i] < tp5) {
            low = i - 1;
            break;
        }
    }

    sp = 0.0;
    for (int i = 100; i >= 0; --i) {
        sp += ls[i];
        if (sp > tp5 && sp - ls[i] < tp5) {
            high = i + 1;
            break;
        }
    }

    if (low < 0) low = 0;
    if (high > 100) high = 100;

    return 0;
}

// Modified from Haploview-v4.2: edu/mit/wi/haploview/FindBlocks.java
// https://github.com/jazzywhit/Haploview
double AppLDB::test_block(int x, int y, int len, const vector<vector<int> > &cild) const
{
    if (len > m_par.maxlen)
        return -1;

    if (cild[y][x] < m_par.low || cild[x][y] < m_par.high)
        return -2;

    auto n = y - x + 1;
    if (n < 4 && len > m_par.max_dist_var[n])
        return -3;

    int strong = 0, recomb = 0;

    for (int i = x; i <= y; ++i) {
        for (int j = i + 1; j <= y; ++j) {
            auto low = cild[j][i], high = cild[i][j];
            if (low < 0 || high < 0)
                continue;

            if (n < 5) {
                if (low > m_par.cut_low_var[n] && high >= m_par.high)
                    ++strong;
            }
            else {
                if (low > m_par.low && high >= m_par.high)
                    ++strong;
            }

            if (high < m_par.rec)
                ++recomb;
        }
    }

    if (n > 3) {
        if (strong + recomb < 6)
            return -4;
    }
    else if (n > 2) {
        if (strong + recomb < 3)
            return -5;
    }
    else {
        if (strong + recomb < 1)
            return -6;
    }

    static const double eps = std::numeric_limits<double>::epsilon();

    auto frac = (double) strong / (strong + recomb);
    if (frac - m_par.frac > eps)
        return frac;

    return -7;
}
