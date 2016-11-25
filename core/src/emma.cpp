#include <cmath>
#include <limits>
#include <memory>
#include <algorithm>
#include "emma.h"
#include "lapack.h"
#include "stat.h"

#ifndef M_2PI
#define M_2PI   6.283185307179586476925286766559    // 2*pi
#endif

// Kang, H.M. et al. Efficient control of population structure in model organism association mapping.
//   Genetics 178, 1709-23 (2008). doi: 10.1534/genetics.107.080101

namespace
{

double LL_REML_wo_Z(size_t nq, double ldelta, double *lambda, double *etas)
{
    using std::log;
    auto delta = std::exp(ldelta);
    double a = 0.0, b = 0.0;
    for (size_t i = 0; i < nq; ++i) {
        a += etas[i]*etas[i]/(lambda[i]+delta);
        b += log(lambda[i]+delta);
    }
    return 0.5*(nq*(log(nq/M_2PI)-1.0-log(a))-b);
}

double dLL_REML_wo_Z(size_t nq, double ldelta, double *lambda, double *etas)
{
    auto delta = std::exp(ldelta);
    double a = 0.0, b = 0.0, c = 0.0;
    for (size_t i = 0; i < nq; ++i) {
        auto sqr = etas[i] * etas[i];
        auto sum = lambda[i] + delta;
        a += sqr/(sum*sum);
        b += sqr/sum;
        c += 1.0/sum;
    }
    return 0.5*(nq*a/b-c);
}

int invsym(size_t n, double *a)
{
    int info = 0;

    std::unique_ptr<bint[]> ipiv(new bint[n]);

    info = la_dsytrf('U', n, a, n, ipiv.get());
    if (info > 0)
        return 1;

    info = la_dsytri('U', n, a, n, ipiv.get());
    if (info > 0)
        return 2;

    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < i; ++j)
            a[j*n+i] = a[i*n+j];
    }

    return 0;
}

int eigen_R_wo_Z(size_t n, size_t q, double *X, double *K, double *eval, double *evec)
{
    auto lwork = q*q + n*q + 3*n*n;
    std::unique_ptr<double[]> work(new double[lwork]);

    auto iXtX   = work.get();
    auto XiXtX  = iXtX + q*q;
    auto S      = XiXtX + n*q;
    auto SK     = S + n*n;
    auto SKS    = SK + n*n;
    auto w      = work.get();
    auto z      = w + n;

    la_dsyrk('U', 'T', q, n, 1.0, X, n, 0.0, iXtX, q);
    if ( invsym(q, iXtX) )
        return 1;

    la_dgemm('N', 'N', n, q, q, 1.0, X, n, iXtX, q, 0.0, XiXtX, n);

    std::fill_n(S, n*n, 0.0);
    for (size_t i = 0; i < n; ++i)
        S[i*n+i] = 1.0;

    la_dgemm('N', 'T', n, n, q, -1.0, XiXtX, n, X, n, 1.0, S, n);

    for (size_t i = 0; i < n; ++i)
        K[i*n+i] += 1.0;

    la_dgemm('N', 'N', n, n, n, 1.0, S, n, K, n, 0.0, SK, n);
    la_dgemm('N', 'N', n, n, n, 1.0, SK, n, S, n, 0.0, SKS, n);

    for (size_t i = 0; i < n; ++i)
        K[i*n+i] -= 1.0;

    bint m = 0;
    std::unique_ptr<bint[]> sup(new bint[2*n]);

    la_dsyevr('V', 'A', 'U', n, SKS, n, 0.0, 0.0, 0, 0, 0.0, &m, w, z, n, sup.get());

    auto nq = n - q;
    for (size_t j = 0; j < nq; ++j) {
        size_t k = n - j - 1;
        eval[j] = w[k] - 1.0;
        for (size_t i = 0; i < n; ++i)
            evec[j*n+i] = z[k*n+i];
    }

    return 0;
}

struct dLL_REML_wo_Z_wrapper_info
{
    double *lambda;
    double *etas;
    size_t nq;
};

double dLL_REML_wo_Z_wrapper(double ldelta, void *info)
{
    auto p = static_cast<dLL_REML_wo_Z_wrapper_info*>(info);
    return dLL_REML_wo_Z(p->nq, ldelta, p->lambda, p->etas);
}

} // namespace

void EMMA::solve(const vector<double> &y, const vector< vector<double> > &X, const vector< vector<double> > &K)
{
    m_n = y.size();
    m_q = X.size();

    m_y = y;

    m_X.clear();
    for (auto &e : X)
        m_X.insert(m_X.end(), e.begin(), e.end());

    m_K.clear();
    for (auto &e : K)
        m_K.insert(m_K.end(), e.begin(), e.end());

    m_vc.REML = m_vc.delta = m_vc.ve = m_vc.vg = std::numeric_limits<double>::quiet_NaN();

    m_L.clear();
    m_Ly.clear();

    compute();
}

void EMMA::compute()
{
    auto n = m_n;
    auto q = m_q;
    auto nq = n - q;
    auto m = m_par.ngrids + 1;

    auto lwork = nq + n*nq + nq + 3*m;
    std::unique_ptr<double[]> work(new double[lwork]);

    auto eval   = work.get();
    auto evec   = eval + nq;
    auto etas   = evec + n*nq;
    auto ldelta = etas + nq;
    auto LL     = ldelta + m;
    auto dLL    = LL + m;

    if ( eigen_R_wo_Z(n, q, m_X.data(), m_K.data(), eval, evec) )
        return;

    la_dgemv('T', n, nq, 1.0, evec, n, m_y.data(), 1, 0.0, etas, 1);

    for (int i = 0; i < m; ++i) {
        ldelta[i] = m_par.llim + i*(m_par.ulim-m_par.llim)/m_par.ngrids;
        LL[i] =  LL_REML_wo_Z(nq, ldelta[i], eval, etas);
        dLL[i] = dLL_REML_wo_Z(nq, ldelta[i], eval, etas);
        if (std::isnan(m_vc.REML) || LL[i] > m_vc.REML) {
            m_vc.REML = LL[i];
            m_vc.delta = ldelta[i];
        }
    }

    dLL_REML_wo_Z_wrapper_info info = {eval, etas, nq};

    for (int i = 0; i < m_par.ngrids; ++i) {
        if (dLL[i] > 0.0 && dLL[i+1] < 0.0 && !std::isnan(LL[i])) {
            double tol = m_par.tol;
            int maxit = m_par.maxit;
            auto root = zeroin(ldelta[i], ldelta[i+1], dLL[i], dLL[i+1], dLL_REML_wo_Z_wrapper, &info, tol, maxit);
            auto optLL = LL_REML_wo_Z(nq, root, eval, etas);
            if ( std::isnan(m_vc.REML) || optLL > m_vc.REML ) {
                m_vc.REML = optLL;
                m_vc.delta = root;
            }
        }
    }

    if ( ! std::isnan(m_vc.REML) ) {
        m_vc.delta = std::exp(m_vc.delta);
        m_vc.vg = 0.0;
        for (size_t i = 0; i < nq; ++i)
            m_vc.vg += etas[i]*etas[i]/(eval[i]+m_vc.delta);
        m_vc.vg /= nq;
        m_vc.ve = m_vc.vg*m_vc.delta;
    }
}

int EMMA::cholfact()
{
    int info = 0;
    auto n = m_n;

    m_L = m_K;

    info = la_dscal(m_L.size(), m_vc.vg, m_L.data(), 1);

    for (size_t i = 0; i < n; ++i)
        m_L[i*n+i] += m_vc.ve;

    info = la_dpotrf('L', n, m_L.data(), n);
    if (info > 0)
        return 1;

    m_Ly = m_y;

    info = la_dtrtrs('L', 'N', 'N', n, 1, m_L.data(), n, m_Ly.data(), n);
    if (info > 0)
        return 2;

    return 0;
}

void EMMA::ftest(size_t p, const vector<vector<double> > &X, double &fval, double &pval)
{
    int info = 0;

    fval = pval = std::numeric_limits<double>::quiet_NaN();

    if ( std::isnan(m_vc.delta) )
        return;

    if ( m_L.empty() && cholfact() )
        return;

    auto n = m_n;
    auto q = X.size();

    m_LX.clear();
    for (auto &e : X)
        m_LX.insert(m_LX.end(), e.begin(), e.end());

    info = la_dtrtrs('L', 'N', 'N', n, q, m_L.data(), n, m_LX.data(), n);
    if (info > 0)
        return;

    auto lwork = q*q + 2*q + p*p;
    std::unique_ptr<double[]> work(new double[lwork]);

    auto iXtX = work.get();
    auto Xty = iXtX + q*q;
    auto b = Xty + q;
    auto iMiM = b + q;

    la_dsyrk('U', 'T', q, n, 1.0, m_LX.data(), n, 0.0, iXtX, q);

    if ( invsym(q, iXtX) )
        return;

    la_dgemv('T', n, q, 1.0, m_LX.data(), n, m_Ly.data(), 1, 0.0, Xty, 1);

    la_dgemv('N', q, q, 1.0, iXtX, q, Xty, 1, 0.0, b, 1);

    if (p == 1) {
        fval = b[q-1]*b[q-1] / iXtX[q*q-1];
    }
    else {
        for (size_t i = 0; i < p; ++i)
            std::copy_n(iXtX + (q-p+i)*q + (q-p), p, iMiM + i*p);
        if ( invsym(p, iMiM) )
            return;
        auto Mb = b + q-p;
        auto iMb = Xty;
        la_dgemv('N', p, p, 1.0, iMiM, p, Mb, 1, 0.0, iMb, 1);
        fval = la_ddot(p, Mb, 1, iMb, 1) / p;
    }

    pval = fcdf(fval, p, n-q, false);
}

void EMMA::ftest(size_t p, const vector<double> &y, const vector<vector<double> > &X,
                 const vector<vector<double> > &K, double &fval, double &pval)
{
    int info = 0;

    fval = pval = std::numeric_limits<double>::quiet_NaN();

    if ( std::isnan(m_vc.delta) )
        return;

    auto n = y.size();
    auto q = X.size();

    vector<double> L;
    for (auto &e : K)
        L.insert(L.end(), e.begin(), e.end());

    la_dscal(L.size(), m_vc.vg, L.data(), 1);

    for (size_t i = 0; i < n; ++i)
        L[i*n+i] += m_vc.ve;

    info = la_dpotrf('L', n, L.data(), n);
    if (info > 0)
        return;

    auto Ly = y;

    info = la_dtrtrs('L', 'N', 'N', n, 1, L.data(), n, Ly.data(), n);
    if (info > 0)
        return;

    m_LX.clear();
    for (auto &e : X)
        m_LX.insert(m_LX.end(), e.begin(), e.end());

    info = la_dtrtrs('L', 'N', 'N', n, q, L.data(), n, m_LX.data(), n);
    if (info > 0)
        return;

    auto lwork = q*q + 2*q + p*p;
    std::unique_ptr<double[]> work(new double[lwork]);

    auto iXtX = work.get();
    auto Xty = iXtX + q*q;
    auto b = Xty + q;
    auto iMiM = b + q;

    la_dsyrk('U', 'T', q, n, 1.0, m_LX.data(), n, 0.0, iXtX, q);

    if ( invsym(q, iXtX) )
        return;

    la_dgemv('T', n, q, 1.0, m_LX.data(), n, Ly.data(), 1, 0.0, Xty, 1);

    la_dgemv('N', q, q, 1.0, iXtX, q, Xty, 1, 0.0, b, 1);

    if (p == 1) {
        fval = b[q-1]*b[q-1] / iXtX[q*q-1];
    }
    else {
        for (size_t i = 0; i < p; ++i)
            std::copy_n(iXtX + (q-p+i)*q + (q-p), p, iMiM + i*p);
        if ( invsym(p, iMiM) )
            return;
        auto Mb = b + q-p;
        auto iMb = Xty;
        la_dgemv('N', p, p, 1.0, iMiM, p, Mb, 1, 0.0, iMb, 1);
        fval = la_ddot(p, Mb, 1, iMb, 1) / p;
    }

    pval = fcdf(fval, p, n-q, false);
}
