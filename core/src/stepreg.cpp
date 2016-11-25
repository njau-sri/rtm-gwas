#include <memory>
#include <iostream>
#include <algorithm>
#include "stepreg.h"
#include "lstsqr.h"
#include "stat.h"

void StepReg::forward(const vector<double> &y)
{
    auto n = y.size();
    auto m = m_cols.size();

    auto sst = sumsqc(y);

    vector<double*> xs;
    xs.push_back(m_effs.data());
    for (size_t i = 1; i < m; ++i)
        xs.push_back(xs[i-1] + n*m_cols[i-1]);

    auto lstsqr = std::make_shared<LstSqr>();
    vector<double> x;

    auto q0 = 1 + m_cofs.size() / n;
    x.assign(n, 1.0);
    x.insert(x.end(), m_cofs.begin(), m_cofs.end());

    vector<bool> ignore(m, false);
    for (auto i : m_model) {
        x.insert(x.end(), xs[i], xs[i] + n*m_cols[i]);
        ignore[i] = true;
        q0 += m_cols[i];
    }

    lstsqr->solve(q0, y, x);
    auto dfe0 = lstsqr->stats().dfe;
    auto sse0 = lstsqr->stats().sse;

    size_t nx = std::count(ignore.begin(), ignore.end(), false);

    size_t maxstep = m_par.maxstep;
    if (maxstep <= 0 || maxstep > nx)
        maxstep = nx;

    std::cerr << "INFO: model selection forward step 0: " << 1-sse0/sst << "\n";

    for (size_t step = 0; step < maxstep; ++step) {
        size_t idx = 0;
        double dfe = 0.0, sse = 0.0, pval = 1.0;

        for (size_t j = 0; j < m; ++j) {
            if ( ignore[j] )
                continue;

            size_t q1 = q0 + m_cols[j];
            x.resize(n*q0);
            x.insert(x.end(), xs[j], xs[j] + n*m_cols[j]);

            lstsqr->solve(q1, y, x);
            auto dfe1 = lstsqr->stats().dfe;
            auto sse1 = lstsqr->stats().sse;

            if (dfe1 > 0.0 && sse1 > 0.0 && dfe0 > dfe1 && sse0 > sse1) {
                auto f = ((sse0-sse1)/(dfe0-dfe1)) / (sse1/dfe1);
                auto p = fcdf(f, dfe0-dfe1, dfe1, false);
                if (p < pval) {
                    idx = j;
                    pval = p;
                    dfe = dfe1;
                    sse = sse1;
                }
            }
        }

        double alpha = m_par.sle;
        if (m_par.multtest == 1)
            alpha /= nx;
        else if (m_par.multtest == 2)
            alpha = alpha * (step + 1) / nx;
        else if (m_par.multtest == 3)
            alpha /= (nx - step);

        if (pval > alpha)
            break;

        if (1-sse/sst > m_par.maxrsq)
            break;

        m_model.push_back(idx);
        ignore[idx] = true;

        x.resize(n*q0);
        x.insert(x.end(), xs[idx], xs[idx] + n*m_cols[idx]);
        q0 += m_cols[idx];
        dfe0 = dfe;
        sse0 = sse;

        std::cerr << "INFO: model selection forward step " << step+1 << ": " << pval << " " << 1-sse0/sst << "\n";
    }
}

void StepReg::backward(const vector<double> &y)
{
    auto n = y.size();
    auto m = m_cols.size();

    vector<double*> xs;
    xs.push_back(m_effs.data());
    for (size_t i = 1; i < m; ++i)
        xs.push_back(xs[i-1] + n*m_cols[i-1]);

    auto lstsqr = std::make_shared<LstSqr>();
    vector<double> x, ps;

    x.assign(n, 1.0);
    x.insert(x.end(), m_cofs.begin(), m_cofs.end());

    int step = 0;
    while ( ! m_model.empty() ) {
        x.resize(n + m_cofs.size());
        size_t q1 = 1 + m_cofs.size() / n;
        for (auto i : m_model) {
            x.insert(x.end(), xs[i], xs[i] + n*m_cols[i]);
            q1 += m_cols[i];
        }

        lstsqr->solve(q1, y, x);

        auto dfe1 = lstsqr->stats().dfe;
        auto sse1 = lstsqr->stats().sse;

        if (dfe1 == 0.0 || sse1 == 0.0)
            return;

        ps.clear();

        for (auto i : m_model) {
            x.resize(n + m_cofs.size());
            size_t q0 = 1 + m_cofs.size() / n;
            for (auto j : m_model) {
                if (j != i) {
                    x.insert(x.end(), xs[j], xs[j] + n*m_cols[j]);
                    q0 += m_cols[j];
                }
            }

            lstsqr->solve(q0, y, x);

            auto dfe0 = lstsqr->stats().dfe;
            auto sse0 = lstsqr->stats().sse;

            if (sse0 > sse1 && dfe0 > dfe1) {
                auto fval = ((sse0-sse1)/(dfe0-dfe1)) / (sse1/dfe1);
                auto pval = fcdf(fval, dfe0-dfe1, dfe1, false);
                ps.push_back(pval);
            }
            else {
                ps.push_back(1.0);
            }
        }

        auto itr = std::max_element(ps.begin(), ps.end());
        if (*itr <= m_par.sls)
            break;

        m_model.erase(m_model.begin() + (itr - ps.begin()));
        std::cerr << "INFO: model selection backward step " << ++step << ": " << *itr << "\n";
    }
}

void StepReg::stepwise(const vector<double> &y)
{
    forward(y);
    backward(y);
}
