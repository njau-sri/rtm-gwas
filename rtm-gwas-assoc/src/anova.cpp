#include <limits>
#include <sstream>
#include <algorithm>
#include <functional>
#include "anova.h"
#include "statsutil.h"
#include "lsfit.h"


using std::size_t;


namespace {


template<typename T1, typename T2>
std::vector<T1> subset(const std::vector<T1> &vec, const std::vector<T2> &idx)
{
    std::vector<T1> out;

    out.reserve(idx.size());

    for (auto i : idx)
        out.push_back(vec[i]);

    return out;
}


} // namespace


void ANOVA::add_reg(const std::string &name, const std::vector<double> &x)
{
    auto tm = std::make_shared<Term>();
    tm->name = name;
    tm->par.push_back(name);
    tm->dat.push_back(x);
    tms_.push_back(tm);
}

void ANOVA::add_main(const std::string &name, const std::vector<std::string> &a)
{
    std::vector<int> gi;
    std::vector<std::string> gn;
    factor(a, gn, gi);

    std::vector< std::vector<double> > x;
    idummy3(gi, x);

    auto tm = std::make_shared<Term>();
    tm->name = name;
    for (auto &e : gn)
        tm->par.push_back(name + " " + e);
    tm->contr.emplace_back(gn.size(), 1.0);
    tm->dat.swap(x);
    tms_.push_back(tm);
}

void ANOVA::add_crossed(const std::string &name, const std::vector<std::string> &a, const std::vector<std::string> &b)
{
    std::vector<int> gia;
    std::vector<std::string> gna;
    factor(a, gna, gia);

    std::vector<int> gib;
    std::vector<std::string> gnb;
    factor(b, gnb, gib);

    std::vector< std::vector<double> > xa;
    idummy3(gia, xa);

    std::vector< std::vector<double> > xb;
    idummy3(gib, xb);

    auto tm = std::make_shared<Term>();
    tm->name = name;

    auto na = gna.size();
    auto nb = gnb.size();

    // A1*B1 A1*B2 A1*B3, A2*B1 A2*B2 A2*B3
    for (size_t i = 0; i < na; ++i) {
        for (size_t j = 0; j < nb; ++j) {
            tm->par.push_back(name + " " + gna[i] + "*" + gnb[j]);
            auto v = xb[j];
            std::transform(v.begin(), v.end(), xa[i].begin(), v.begin(), std::multiplies<double>());
            tm->dat.push_back(v);
        }
    }

    // zero-sum constraints for A
    for (size_t j = 0; j < nb; ++j) {
        std::vector<double> v(na*nb, 0);
        for (size_t i = 0; i < na; ++i)
            v[i*nb+j] = 1.0;
        tm->contr.push_back(v);
    }

    // zero-sum constraints for B
    for (size_t i = 0; i < na; ++i) {
        std::vector<double> v(na*nb, 0);
        std::fill_n(&v[i*nb], nb, 1.0);
        tm->contr.push_back(v);
    }

    tms_.push_back(tm);
}

void ANOVA::add_nested(const std::string &name, const std::vector<std::string> &a, const std::vector<std::string> &b)
{
    std::vector<int> gia;
    std::vector<std::string> gna;
    factor(a, gna, gia);

    std::vector<int> gib;
    std::vector<std::string> gnb;
    factor(b, gnb, gib);

    std::vector< std::vector<double> > xa;
    idummy3(gia, xa);

    std::vector< std::vector<double> > xb;
    idummy3(gib, xb);

    auto tm = std::make_shared<Term>();
    tm->name = name;

    auto na = gna.size();
    auto nb = gnb.size();

    // A1(B1) A2(B1) A3(B1), A1(B2) A2(B2) A3(B2)
    for (size_t i = 0; i < nb; ++i) {
        for (size_t j = 0; j < na; ++j) {
            tm->par.push_back(name + " " + gna[j] + "(" + gnb[i] + ")");
            auto v = xa[j];
            std::transform(v.begin(), v.end(), xb[i].begin(), v.begin(), std::multiplies<double>());
            tm->dat.push_back(v);
        }
        // zero-sum constraints for A
        std::vector<double> v(na*nb, 0);
        std::fill_n(&v[i*na], na, 1.0);
        tm->contr.push_back(v);
    }

    auto m = tm->dat.size();
    std::vector<size_t> idx;

    for (size_t i = 0; i < m; ++i) {
        for (auto e : tm->dat[i]) {
            if (static_cast<int>(e) != 0) {
                idx.push_back(i);
                break;
            }
        }
    }

    subset(tm->par,idx).swap(tm->par);
    subset(tm->dat,idx).swap(tm->dat);
    for (auto &v : tm->contr)
        subset(v,idx).swap(v);

    tms_.push_back(tm);
}

ANOVA::Table ANOVA::solve1(const std::vector<double> &y) const
{
    auto n = y.size();

    Table tbl;
    std::vector<double> x;

    size_t q0 = 1;
    x.assign(n, 1);

    double dfe0 = n - 1;
    double sse0 = calc_css(y);

    tbl.total.push_back(dfe0);
    tbl.total.push_back(sse0);

    for (auto &tm : tms_) {
        auto q1 = q0 + tm->dat.size();
        for (auto &e : tm->dat)
            x.insert(x.end(), e.begin(), e.end());

        double dfe1, sse1;
        std::vector<double> b1;
        lsfit(y, x, b1, dfe1, sse1);

        tbl.src.push_back(tm->name);
        tbl.df.push_back(dfe0 - dfe1);
        tbl.ss.push_back(sse0 - sse1);

        q0 = q1;
        dfe0 = dfe1;
        sse0 = sse1;
    }

    tbl.error.push_back(dfe0);
    tbl.error.push_back(sse0);
    tbl.error.push_back(sse0 / dfe0);

    auto m = tms_.size();
    tbl.ms.resize(m);
    tbl.f.resize(m);
    tbl.p.assign(m, std::numeric_limits<double>::quiet_NaN());
    for (size_t j = 0; j < m; ++j) {
        tbl.ms[j] = tbl.ss[j] / tbl.df[j];
        tbl.f[j] = tbl.ms[j] / tbl.error[2];
        if (tbl.df[j] > 0.0 && tbl.ss[j] > 0.0)
            tbl.p[j] = fpval(tbl.f[j], tbl.df[j], tbl.error[0]);
    }

    return tbl;
}

ANOVA::Table ANOVA::solve3(const std::vector<double> &y) const
{
    auto n = y.size();

    Table tbl;
    std::vector<double> x;

    size_t q1 = 1;
    x.assign(n, 1);

    double dft = n - 1;
    double sst = calc_css(y);

    tbl.total.push_back(dft);
    tbl.total.push_back(sst);

    for (auto &tm : tms_) {
        q1 += tm->dat.size();
        for (auto &e : tm->dat)
            x.insert(x.end(), e.begin(), e.end());
    }

    std::vector<double> zt;
    int pos = 1, nz = 0;
    for (auto &tm : tms_) {
        for (auto &e : tm->contr) {
            std::vector<double> v(q1, 0);
            std::copy(e.begin(), e.end(), v.begin() + pos);
            zt.insert(zt.end(), v.begin(), v.end());
            nz += 1;
        }
        pos += tm->dat.size();
    }

    double dfe, sse;
    std::vector<double> b;
    lsfitc(y, x, zt, b, dfe, sse);

    tbl.error.push_back(dfe);
    tbl.error.push_back(sse);
    tbl.error.push_back(sse / dfe);

    for (auto &curr : tms_) {
        size_t q0 = 1;
        x.assign(n, 1);

        for (auto &tm : tms_) {
            if (&tm == &curr)
                continue;
            q0 += tm->dat.size();
            for (auto &e : tm->dat)
                x.insert(x.end(), e.begin(), e.end());
        }

        zt.clear();
        pos = 1;
        nz = 0;
        for (auto &tm : tms_) {
            if (&tm == &curr)
                continue;
            for (auto &e : tm->contr) {
                std::vector<double> v(q0, 0);
                std::copy(e.begin(), e.end(), v.begin() + pos);
                zt.insert(zt.end(), v.begin(), v.end());
                nz += 1;
            }
            pos += tm->dat.size();
        }

        double dfe0, sse0;
        std::vector<double> b0;
        lsfitc(y, x, zt, b0, dfe0, sse0);

        tbl.src.push_back(curr->name);
        tbl.df.push_back(dfe0 - dfe);
        tbl.ss.push_back(sse0 - sse);
    }

    auto m = tms_.size();
    tbl.ms.resize(m);
    tbl.f.resize(m);
    tbl.p.assign(m, std::numeric_limits<double>::quiet_NaN());
    for (size_t j = 0; j < m; ++j) {
        tbl.ms[j] = tbl.ss[j] / tbl.df[j];
        tbl.f[j] = tbl.ms[j] / tbl.error[2];
        if (tbl.df[j] > 0.0 && tbl.ss[j] > 0.0)
            tbl.p[j] = fpval(tbl.f[j], tbl.df[j], tbl.error[0]);
    }

    return tbl;
}

ANOVA::Solution ANOVA::solution(const std::vector<double> &y) const
{
    auto n = y.size();

    std::vector<double> x;
    std::vector<std::string> par;

    size_t p = 1;
    x.assign(n, 1);
    par.push_back("Constant");

    for (auto &tm : tms_) {
        par.insert(par.end(), tm->par.begin(), tm->par.end());
        p += tm->dat.size();
        for (auto &e : tm->dat)
            x.insert(x.end(), e.begin(), e.end());
    }

    std::vector<double> zt;
    int pos = 1, q = 0;
    for (auto &tm : tms_) {
        for (auto &e : tm->contr) {
            std::vector<double> v(p, 0);
            std::copy(e.begin(), e.end(), v.begin() + pos);
            zt.insert(zt.end(), v.begin(), v.end());
            q += 1;
        }
        pos += tm->dat.size();
    }

    double dfe, sse;
    std::vector<double> b;
    lsfitc(y, x, zt, b, dfe, sse);

    Solution sol;
    sol.par.swap(par);
    sol.est.swap(b);

    return sol;
}

std::string ANOVA::Table::to_string() const
{
    std::ostringstream oss;
    oss << "Source\tDF\tSS\tMS\tF\tp\n";

    auto n = src.size();
    for (size_t i = 0; i < n; ++i)
        oss << src[i] << "\t" << df[i] << "\t" << ss[i] << "\t" << ms[i] << "\t" << f[i] << "\t" << p[i] << "\n";

    oss << "Error";
    for (auto &e : error)
        oss << "\t" << e;
    oss << "\n";

    oss << "Total";
    for (auto &e : total)
        oss << "\t" << e;
    oss << "\n";

    return oss.str();
}

std::string ANOVA::Solution::to_string() const
{
    std::ostringstream oss;
    oss << "Parameter\tEstimate\n";

    auto n = par.size();
    for (size_t i = 0; i < n; ++i)
        oss << par[i] << "\t" << est[i] << "\n";

    return oss.str();
}
