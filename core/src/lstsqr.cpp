#include <cmath>
#include <algorithm>
#include "lstsqr.h"
#include "lapack.h"

int LstSqr::solve(size_t p, const vector<double> &y, const vector<double> &x)
{
    static const double tol = 1e-8;

    int info = 0;
    auto n = y.size();
    auto ldy  = std::max(n, p);
    auto rank = std::min(n, p);

    m_st.sse = m_st.dfe = 0.0;

    m_x = x;
    m_y.assign(ldy, 0.0);
    std::copy(y.begin(), y.end(), m_y.begin());

    bool singular = true;

    if (n >= p) {
        singular = false;

        info = la_dgels('N', n, p, 1, m_x.data(), n, m_y.data(), ldy);

        if (info == 0) {
            for (size_t i = 0; i < p; ++i) {
                if (std::fabs(m_x[i*n+i]) < tol) {
                    singular = true;
                    break;
                }
            }
        }
        else
            singular = true;
    }

    if ( singular ) {
        m_x = x;
        m_y.assign(ldy, 0.0);
        std::copy(y.begin(), y.end(), m_y.begin());

        bint lrank = 0;
        vector<bint> jpvt(p,0);
        double rcond = tol;

        info = la_dgelsy(n, p, 1, m_x.data(), n, m_y.data(), ldy, jpvt.data(), rcond, &lrank);

        rank = lrank;
    }

    m_st.b.assign(m_y.begin(), m_y.begin() + p);

    if (n > rank) {
        m_y = y;
        la_dgemv('N', n, p, -1.0, x.data(), n, m_st.b.data(), 1, 1.0, m_y.data(), 1);
        m_st.sse = la_dnrm2(n, m_y.data(), 1);
        m_st.sse *= m_st.sse;
        m_st.dfe = n - rank;
    }

    return 0;
}

int LstSqr::solve(size_t p, size_t q, const vector<double> &y, const vector<double> &x, vector<double> ct)
{
    static const double tol = 1e-8;

    auto n = y.size();
    auto minpq = std::min(p,q);

    vector<bint> jpvt(q,0);
    vector<double> tau(minpq);

    la_dgeqp3(p, q, ct.data(), p, jpvt.data(), tau.data());

    size_t k = 0;
    for (size_t i = 0; i < minpq; ++i) {
        if (std::fabs(ct[i*p+i]) > tol)
            ++k;
    }

    vector<double> Q(p*p);
    std::copy_n(ct.begin(), p*minpq, Q.begin());

    la_dorgqr(p, p, minpq, Q.data(), p, tau.data());

    auto Q0 = Q.data() + p*k;

    vector<double> X(n*(p-k));
    la_dgemm('N', 'N', n, p-k, p, 1.0, x.data(), n, Q0, p, 0.0, X.data(), n);

    solve(p-k, y, X);

    vector<double> b(p);
    la_dgemv('N', p, p-k, 1.0, Q0, p, m_st.b.data(), 1, 0.0, b.data(), 1);

    m_st.b = b;

    return 0;
}
