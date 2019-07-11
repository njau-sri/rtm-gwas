#include <cmath>
#include <memory>
#include <algorithm>
#include "lsfit.h"
#include "lapack.h"


using std::size_t;


namespace {

static const double RCOND_DGELSY_LSFIT = 1e-8;

static const double TOL_DGEQP3_LSFITC = 1e-8;

} // namespace


int lsfit(const std::vector<double> &y, const std::vector<double> &x, std::vector<double> &b, double &dfe, double &sse)
{
    auto n = y.size();
    auto q = x.size() / n;
    auto ldy = std::max(n, q);

    b.clear();
    dfe = sse = 0.0;

    auto cx = x;
    std::vector<double> cy(ldy, 0.0);
    std::copy(y.begin(), y.end(), cy.begin());

    bint irank = 0;
    std::vector<bint> jpvt(q, 0);
    double rcond = RCOND_DGELSY_LSFIT;

    int info = C_dgelsy(n, q, 1, cx.data(), n, cy.data(), ldy, jpvt.data(), rcond, &irank);

    if (info != 0)
        return 1;

    auto rank = static_cast<size_t>(irank);

    b = cy;
    b.resize(q);

    if (rank < n) {
        dfe = n - rank;
        cy = y;
        C_dgemv('N', n, q, -1.0, x.data(), n, b.data(), 1, 1.0, cy.data(), 1);
        sse = C_dnrm2(n, cy.data(), 1);
        sse *= sse;
    }

    return 0;
}


int lsfitc(const std::vector<double> &y, const std::vector<double> &x, std::vector<double> z,
           std::vector<double> &b, double &dfe, double &sse)
{
    auto n = y.size();
    auto q = x.size() / n;
    auto p = z.size() / q;
    auto m = std::min(q, p);

    std::vector<bint> jpvt(p, 0);
    std::vector<double> tau(m);

    C_dgeqp3(q, p, z.data(), q, jpvt.data(), tau.data());

    size_t k = 0;
    for (size_t i = 0; i < m; ++i) {
        if (std::fabs(z[i*q+i]) > TOL_DGEQP3_LSFITC)
            ++k;
    }

    std::vector<double> Q(q*q);
    std::copy_n(z.begin(), q*m, Q.begin());
    C_dorgqr(q, q, m, Q.data(), q, tau.data());

    auto Q0 = Q.data() + q*k;

    std::vector<double> xq(n*(q-k));
    C_dgemm('N', 'N', n, q-k, q, 1.0, x.data(), n, Q0, q, 0.0, xq.data(), n);

    std::vector<double> cb;
    lsfit(y, xq, cb, dfe, sse);

    b.resize(q);
    C_dgemv('N', q, q-k, 1.0, Q0, q, cb.data(), 1, 0.0, b.data(), 1);

    return 0;
}


int glsfit(const std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &v,
           std::vector<double> &b)
{
    int info = 0;

    auto n = y.size();
    auto q = x.size() / n;

    // LL' = V
    auto L = v;
    info = C_dpotrf('L', n, L.data(), n);
    if (info != 0)
        return 1;

    // L^-1 * y
    auto ly = y;
    C_dtrsv('L', 'N', 'N', n, L.data(), n, ly.data(), 1);

    // L^-1 * y
    auto lx = x;
    C_dtrsm('L', 'L', 'N', 'N', n, q, 1.0, L.data(), n, lx.data(), n);

    // b = (X'*X)^-1 * X' * y
    double dfe = 0.0, sse = 0.0;
    info = lsfit(ly, lx, b, dfe, sse);
    if (info != 0)
        return 2;

    return 0;
}


double glsfstat(std::size_t p, const std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &v)
{
    int info = 0;

    auto n = y.size();
    auto q = x.size() / n;

    if (p > q)
        return -1;

    auto nn = n*n;
    auto qq = q*q;
    auto qn = q*n;

    auto lwork = nn + qn + qq + q*2;
    std::unique_ptr<double[]> work(new double[lwork]);

    auto iV = work.get();
    auto XtV = iV + nn;
    auto iXVX = XtV + qn;
    auto XVy = iXVX + qq;
    auto b = XVy + q;
    auto iMiM = work.get();
    auto ib = iMiM + p*p;

    // V^-1  ->  iV [n,n]
    std::copy_n(v.begin(), nn, iV);
    info = M_dsyinv(n, iV);
    if (info != 0)
        return -2;

    // X' V^-1  ->  XtV [q,n]
    C_dgemm('T', 'N', q, n, n, 1.0, x.data(), n, iV, n, 0.0, XtV, q);

    // X' V^-1 X  ->  iXX [q,q]
    C_dgemm('N', 'N', q, q, n, 1.0, XtV, q, x.data(), n, 0.0, iXVX, q);

    // (X' V^-1 X)^-1  ->  iXX [q,q]
    info = M_dsyinv(q, iXVX);
    if (info != 0)
        return -3;

    // X' * V^-1 * y  ->  Xy [q,1]
    C_dgemv('N', q, n, 1.0, XtV, q, y.data(), 1, 0.0, XVy, 1);

    // (X' * V^-1 * X)^-1 * X' * V^-1 * y  ->  b [q,1]
    C_dgemv('N', q, q, 1.0, iXVX, q, XVy, 1, 0.0, b, 1);

    // M (X' V^-1 X)^-1 M'  ->  iMiM [p,p]
    //   M = [0 0 0 1 0; 0 0 0 0 1]
    //   bottom right corner p*p submatrix of X' V^-1 X
    if (p == 1) {
        iMiM[0] = iXVX[qq-1];
    }
    else {
        auto k = q - p;
        for (size_t i = 0; i < p; ++i) {
            auto s1 = i * p;
            auto s2 = (k + i) * q;
            for (size_t j = 0; j < p; ++j)
                iMiM[s1+j] = iXVX[s2+k+j];
        }
    }

    // [M (X' V^-1 X)^-1 M']^-1  ->  iMiM
    info = M_dsyinv(p, iMiM);
    if (info != 0)
        return -4;

    // Mb [p,1]
    auto Mb = b + q - p;

    // [M (X' V^-1 X)^-1 M']^-1 (Mb)  ->  ib [p,1]
    C_dgemv('N', p, p, 1.0, iMiM, p, Mb, 1, 0.0, ib, 1);

    // (Mb)' [M (X' V^-1 X)^-1 M']^-1 (Mb)
    auto f = C_ddot(p, Mb, 1, ib, 1);
    f /= p;

    return f;
}
