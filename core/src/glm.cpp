#include <memory>
#include <algorithm>
#include <functional>
#include "glm.h"
#include "util.h"
#include "stat.h"
#include "lstsqr.h"

void GLM::add_reg(const string &name, const vector<double> &x)
{
    auto tm = std::make_shared<ModelTerm>();
    tm->name = name;
    tm->coeffnames.push_back(name);
    tm->data.push_back(x);
    m_terms.push_back(tm);
}

void GLM::add_main(const string &name, const vector<string> &a)
{
    vector<int> g;
    vector<string> gn;
    factor(a, gn, g);

    vector< vector<double> > x;
    design(3, g, x);

    auto tm = std::make_shared<ModelTerm>();
    tm->name = name;
    for (auto &e : gn)
        tm->coeffnames.push_back(name + "\t" + e);
    tm->constr.emplace_back(gn.size(), 1.0);
    tm->data.swap(x);
    m_terms.push_back(tm);
}

void GLM::add_crossed(const string &name, const vector<string> &a, const vector<string> &b)
{
    vector<int> ga;
    vector<string> gna;
    factor(a, gna, ga);

    vector<int> gb;
    vector<string> gnb;
    factor(b, gnb, gb);

    vector< vector<double> > xa;
    design(3, ga, xa);

    vector< vector<double> > xb;
    design(3, gb, xb);

    auto tm = std::make_shared<ModelTerm>();
    tm->name = name;

    auto na = gna.size();
    auto nb = gnb.size();

    // A1*B1 A1*B2 A1*B3, A2*B1 A2*B2 A2*B3
    for (size_t i = 0; i < na; ++i) {
        for (size_t j = 0; j < nb; ++j) {
            tm->coeffnames.push_back(name + "=" + gna[i] + "*" + gnb[j]);
            auto v = xb[j];
            std::transform(v.begin(), v.end(), xa[i].begin(), v.begin(), std::multiplies<double>());
            tm->data.push_back(v);
        }
    }

    // zero-sum constraints for A
    for (size_t j = 0; j < nb; ++j) {
        vector<double> v(na*nb, 0);
        for (size_t i = 0; i < na; ++i)
            v[i*nb+j] = 1.0;
        tm->constr.push_back(v);
    }

    // zero-sum constraints for B
    for (size_t i = 0; i < na; ++i) {
        vector<double> v(na*nb, 0);
        std::fill_n(v.begin() + i*nb, nb, 1.0);
        tm->constr.push_back(v);
    }

    m_terms.push_back(tm);
}

void GLM::add_nested(const string &name, const vector<string> &a, const vector<string> &b)
{
    vector<int> ga;
    vector<string> gna;
    factor(a, gna, ga);

    vector<int> gb;
    vector<string> gnb;
    factor(b, gnb, gb);

    vector< vector<double> > xa;
    design(3, ga, xa);

    vector< vector<double> > xb;
    design(3, gb, xb);

    auto tm = std::make_shared<ModelTerm>();
    tm->name = name;

    auto na = gna.size();
    auto nb = gnb.size();

    // A1(B1) A2(B1) A3(B1), A1(B2) A2(B2) A3(B2)
    for (size_t i = 0; i < nb; ++i) {
        for (size_t j = 0; j < na; ++j) {
            tm->coeffnames.push_back(name + "\t" + gna[j] + "(" + gnb[i] + ")");
            auto v = xa[j];
            std::transform(v.begin(), v.end(), xb[i].begin(), v.begin(), std::multiplies<double>());
            tm->data.push_back(v);
        }
        // zero-sum constraints for A
        vector<double> v(na*nb, 0);
        std::fill_n(v.begin() + i*na, na, 1.0);
        tm->constr.push_back(v);
    }

    m_terms.push_back(tm);
}

GLM::AnovaTable GLM::solve1(const vector<double> &y) const
{
    auto n = y.size();

    auto lstsqr = std::make_shared<LstSqr>();

    AnovaTable atab;
    vector<double> x;

    size_t q0 = 1;
    x.assign(n, 1);

    double df0 = n - 1;
    double rss0 = sumsqc(y);

    atab.total.push_back(df0);
    atab.total.push_back(rss0);

    for (auto &term : m_terms) {
        auto q1 = q0 + term->data.size();
        for (auto &e : term->data)
            x.insert(x.end(), e.begin(), e.end());

        lstsqr->solve(q1, y, x);
        auto df1 = lstsqr->stats().dfe;
        auto rss1 = lstsqr->stats().sse;

        atab.names.push_back(term->name);
        atab.df.push_back(df0 - df1);
        atab.ss.push_back(rss0 - rss1);

        q0 = q1;
        df0 = df1;
        rss0 = rss1;
    }

    atab.error.push_back(df0);
    atab.error.push_back(rss0);
    atab.error.push_back(rss0 / df0);

    auto m = m_terms.size();
    atab.ms.assign(m, 0);
    atab.f.assign(m, 0);
    atab.p.assign(m, 1);
    for (size_t j = 0; j < m; ++j) {
        if (atab.df[j] > 0.0 && atab.ss[j] > 0.0) {
            atab.ms[j] = atab.ss[j] / atab.df[j];
            atab.f[j] = atab.ms[j] / atab.error[2];
            atab.p[j] = fcdf(atab.f[j], atab.df[j], atab.error[0], false);
        }
    }

    return atab;
}

GLM::AnovaTable GLM::solve3(const vector<double> &y) const
{
    auto n = y.size();

    auto lstsqr = std::make_shared<LstSqr>();

    AnovaTable atab;
    vector<double> x;

    size_t q1 = 1;
    x.assign(n, 1);

    double dft = n - 1;
    double sst = sumsqc(y);

    atab.total.push_back(dft);
    atab.total.push_back(sst);

    for (auto &term : m_terms) {
        q1 += term->data.size();
        for (auto &e : term->data)
            x.insert(x.end(), e.begin(), e.end());
    }

    lstsqr->solve(q1, y, x);
    auto dfe = lstsqr->stats().dfe;
    auto sse = lstsqr->stats().sse;

    atab.error.push_back(dfe);
    atab.error.push_back(sse);
    atab.error.push_back(sse / dfe);

    for (auto &curr : m_terms) {
        size_t q0 = 1;
        x.assign(n, 1);

        for (auto &term : m_terms) {
            if (&term == &curr)
                continue;
            q0 += term->data.size();
            for (auto &e : term->data)
                x.insert(x.end(), e.begin(), e.end());
        }

        lstsqr->solve(q0, y, x);
        auto df0 = lstsqr->stats().dfe;
        auto rss0 = lstsqr->stats().sse;

        atab.names.push_back(curr->name);
        atab.df.push_back(df0 - dfe);
        atab.ss.push_back(rss0 - sse);
    }

    auto m = m_terms.size();
    atab.ms.assign(m, 0.0);
    atab.f.assign(m, 0.0);
    atab.p.assign(m, 1.0);
    for (size_t j = 0; j < m; ++j) {
        if (atab.df[j] > 0.0 && atab.ss[j] > 0.0) {
            atab.ms[j] = atab.ss[j] / atab.df[j];
            atab.f[j] = atab.ms[j] / atab.error[2];
            atab.p[j] = fcdf(atab.f[j], atab.df[j], atab.error[0], false);
        }
    }

    return atab;
}

GLM::Coefficient GLM::coeff(const vector<double> &y) const
{
    auto n = y.size();

    vector<string> params;
    vector<double> x;

    size_t p = 1;
    x.assign(n, 1);
    params.push_back("Constant");

    for (auto &term : m_terms) {
        params.insert(params.end(),
            term->coeffnames.begin(),
            term->coeffnames.end());
        p += term->data.size();
        for (auto &e : term->data)
            x.insert(x.end(), e.begin(), e.end());
    }

    vector<double> ct;
    size_t pos = 1, q = 0;
    for (auto &term : m_terms) {
        for (auto &e : term->constr) {
            vector<double> v(p, 0);
            std::copy(e.begin(), e.end(), v.begin() + pos);
            ct.insert(ct.end(), v.begin(), v.end());
            q += 1;
        }
        pos += term->data.size();
    }

    auto lstsqr = std::make_shared<LstSqr>();
    lstsqr->solve(p, q, y, x, ct);

    Coefficient est;
    est.params.swap(params);
    est.coeffs = lstsqr->stats().b;

    return est;
}

void GLM::reset()
{
    m_terms.clear();
    m_terms.shrink_to_fit();
}
