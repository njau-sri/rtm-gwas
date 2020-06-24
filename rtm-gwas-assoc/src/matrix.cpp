#include "matrix.h"

#include <memory>
#include <limits>
#include <algorithm>

#include "print.h"
#include "lapack.h"

void design1(isize_t n, const isize_t *g, mat &x)
{
    auto p = std::minmax_element(g, g + n);

    if (*p.first != 0 || *p.second <= 0) {
        x.resize(0, 0);
        return;
    }

    isize_t m = *p.second;

    x.setZero(n, m);

    for (isize_t i = 0; i < n; ++i) {
        isize_t k = g[i];
        if (k != m)
            x(i, k) = 1.0;
        else
            x.row(i).setConstant(-1.0);
    }
}

void design2(isize_t n, const isize_t *g, mat &x)
{
    auto p = std::minmax_element(g, g + n);

    if (*p.first != 0 || *p.second <= 0) {
        x.resize(0, 0);
        return;
    }

    isize_t m = *p.second;

    x.setZero(n, m);

    for (isize_t i = 0; i < n; ++i) {
        isize_t k = g[i];
        if (k != m)
            x(i, k) = 1.0;
    }
}

void design3(isize_t n, const isize_t *g, mat &x)
{
    auto p = std::minmax_element(g, g + n);

    if (*p.first != 0 || *p.second <= 0) {
        x.resize(0, 0);
        return;
    }

    isize_t m = *p.second;

    x.setZero(n, m+1);

    for (isize_t i = 0; i < n; ++i)
        x(i, g[i]) = 1.0;
}

double nrm2(const Eigen::Ref<const mat> &x)
{
    return call_dnrm2(length(x), x.data(), 1);
}

double dot(const Eigen::Ref<const mat> &x, const Eigen::Ref<const mat> &y)
{
    return call_ddot(length(x), x.data(), 1, y.data(), 1);
}

int gemv(bool trans, double alpha, const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &x, double beta, mat &y)
{
    char at = trans ? 'T' : 'N';
    isize_t m = nrow(a);
    isize_t n = ncol(a);

    if (beta == 0.0)
        y.resize(trans ? n : m, 1);

    isize_t lda = nrow(a);

    call_dgemv(at, m, n, alpha, a.data(), lda, x.data(), 1, beta, y.data(), 1);

    return 0;
}

int gemm(bool transa, bool transb, double alpha, const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b, double beta, mat &c)
{
    char at = 'N';
    char bt = 'N';
    isize_t m = nrow(a);
    isize_t n = ncol(b);
    isize_t k = ncol(a);

    if (transa) {
        at = 'T';
        std::swap(m, k);
    }

    if (transb) {
        bt = 'T';
        n = nrow(b);
    }

    if (beta == 0.0)
        c.resize(m, n);

    isize_t lda = nrow(a);
    isize_t ldb = nrow(b);
    isize_t ldc = nrow(c);

    call_dgemm(at, bt, m, n, k, alpha, a.data(), lda, b.data(), ldb, beta, c.data(), ldc);

    return 0;
}

int syrk(bool trans, double alpha, const Eigen::Ref<const mat> &a, double beta, mat &c)
{
    char at = 'N';
    isize_t n = nrow(a);
    isize_t k = ncol(a);

    if (trans) {
        at = 'T';
        std::swap(n, k);
    }

    if (beta == 0.0)
        c.resize(n, n);

    isize_t lda = nrow(a);
    isize_t ldc = nrow(c);

    call_dsyrk('U', at, n, k, alpha, a.data(), lda, beta, c.data(), ldc);

    for (isize_t i = 0; i < n; ++i)
        for (isize_t j = i + 1; j < n; ++j)
            c(j, i) = c(i, j);

    return 0;
}

mat mult(const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b)
{
    mat c;
    gemm(false, false, 1.0, a, b, 0.0, c);
    return c;
}

mat crossprod(const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b)
{
    mat c;
    gemm(true, false, 1.0, a, b, 0.0, c);
    return c;
}

mat tcrossprod(const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b)
{
    mat c;
    gemm(false, true, 1.0, a, b, 0.0, c);
    return c;
}

mat crossprod(const Eigen::Ref<const mat> &a)
{
    mat c;
    syrk(true, 1.0, a, 0.0, c);
    return c;
}

mat tcrossprod(const Eigen::Ref<const mat> &a)
{
    mat c;
    syrk(false, 1.0, a, 0.0, c);
    return c;
}

double norm1(const Eigen::Ref<const mat> &a)
{
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t lda = nrow(a);
    return call_dlange('1', m, n, a.data(), lda);
}

double normi(const Eigen::Ref<const mat> &a)
{
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t lda = nrow(a);
    return call_dlange('I', m, n, a.data(), lda);
}

double normf(const Eigen::Ref<const mat> &a)
{
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t lda = nrow(a);
    return call_dlange('F', m, n, a.data(), lda);
}

double norm2(const Eigen::Ref<const mat> &a)
{
    mat s;
    int info = gesdd(a, s);

    if (info != 0)
        eprint("ERROR: gesdd failed in norm2 (%s:%d): %d\n", __FILE__, __LINE__, info);

    return s(0, 0);
}

double cond(const Eigen::Ref<const mat> &a)
{
    mat s;

    int info = gesdd(a, s);

    if (info != 0)
        eprint("ERROR: gesdd failed in cond (%s:%d): %d\n", __FILE__, __LINE__, info);

    auto p = std::minmax_element(s.data(), s.data() + length(s));
    return (*p.second) / (*p.first);
}

mat inv_sym(const Eigen::Ref<const mat> &a)
{
    mat b = a;

    int info = dsyinv(nrow(a), b.data());

    if (info != 0)
        eprint("ERROR: dsyinv failed in inv_sym (%s:%d): %d\n", __FILE__, __LINE__, info);

    return b;
}

mat pinv(double rcond, const Eigen::Ref<const mat> &a)
{
    mat x;
    isize_t m = nrow(a);

    int rank = gelsd(rcond, a, mat::Identity(m, m), x);

    if (rank < 0)
        eprint("ERROR: gelsd failed in pinv (%s:%d): %d\n", __FILE__, __LINE__, rank);

    return x;
}

int syev(const Eigen::Ref<const mat> &a, mat &eval)
{
    mat aa = a;
    isize_t n = nrow(a);
    eval.resize(n, 1);

    int info = call_dsyev('N', 'U', n, aa.data(), n, eval.data());

    if (info != 0)
        eprint("ERROR: call_dsyev failed in syev (%s:%d): %d", __FILE__, __LINE__, info);

    return info;
}

int syev(const Eigen::Ref<const mat> &a, mat &eval, mat &evec)
{
    isize_t n = nrow(a);
    eval.resize(n, 1);
    evec = a;

    int info = call_dsyev('V', 'U', n, evec.data(), n, eval.data());

    if (info != 0)
        eprint("ERROR: call_dsyev failed in syev (%s:%d): %d", __FILE__, __LINE__, info);

    return info;
}

int syevr(const Eigen::Ref<const mat> &a, mat &eval)
{
    mat aa = a;
    isize_t n = nrow(a);
    bint m = 0;
    std::unique_ptr<bint[]> sup(new bint[n*2]);
    eval.resize(n, 1);

    int info = call_dsyevr('N', 'A', 'U', n, aa.data(), n, 0, 0, 0, 0, 0, &m, eval.data(), nullptr, n, sup.get());

    if (info != 0)
        eprint("ERROR: call_dsyevr failed in syevr (%s:%d): %d", __FILE__, __LINE__, info);

    return info;
}

int syevr(const Eigen::Ref<const mat> &a, mat &eval, mat &evec)
{
    mat aa = a;
    isize_t n = nrow(a);
    bint m = 0;
    std::unique_ptr<bint[]> sup(new bint[n*2]);
    eval.resize(n, 1);
    evec.resize(n, n);

    int info = call_dsyevr('V', 'A', 'U', n, aa.data(), n, 0, 0, 0, 0, 0, &m, eval.data(), evec.data(), n, sup.get());

    if (info != 0)
        eprint("ERROR: call_dsyevr failed in syevr (%s:%d): %d", __FILE__, __LINE__, info);

    return info;
}

int qr(const Eigen::Ref<const mat> &a, mat &q, mat &r)
{
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t k = m < n ? m : n;

    r = a;
    std::unique_ptr<double[]> tau(new double[k]);

    int info = call_dgeqrf(m, n, r.data(), m, tau.get());

    q.resize(m, m);
    q.leftCols(k) = r.leftCols(k);

    if (info != 0) {
        eprint("ERROR: call_dgeqrf failed in qr (%s:%d): %d", __FILE__, __LINE__, info);
        return info;
    }

    info = call_dorgqr(m, m, k, q.data(), m, tau.get());

    if (info != 0) {
        eprint("ERROR: call_dorgqr failed in qr (%s:%d): %d", __FILE__, __LINE__, info);
        return info;
    }

    r.triangularView<Eigen::StrictlyLower>().setZero();

    return 0;
}

int qr_eco(const Eigen::Ref<const mat> &a, mat &q, mat &r)
{
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t k = m < n ? m : n;

    if (m <= n)
        return qr(a, q, r);

    q = a;
    std::unique_ptr<double[]> tau(new double[k]);

    int info = call_dgeqrf(m, n, q.data(), m, tau.get());

    r = q.topRows(n).triangularView<Eigen::Upper>();

    if (info != 0) {
        eprint("ERROR: call_dgeqrf failed in qr_eco (%s:%d): %d", __FILE__, __LINE__, info);
        return info;
    }

    info = call_dorgqr(m, n, k, q.data(), m, tau.get());

    if (info != 0)
        eprint("ERROR: call_dorgqr failed in qr_eco (%s:%d): %d", __FILE__, __LINE__, info);

    return info;
}

int qr(const Eigen::Ref<const mat> &a, mat &q, mat &r, imat &p)
{
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t k = m < n ? m : n;

    p.setZero(n, 1);

    r = a;
    std::unique_ptr<bint[]> jpvt(new bint[n]);
    std::fill_n(jpvt.get(), n, 0);
    std::unique_ptr<double[]> tau(new double[k]);

    int info = call_dgeqp3(m, n, r.data(), m, jpvt.get(), tau.get());

    q.resize(m, m);
    q.leftCols(k) = r.leftCols(k);

    if (info != 0) {
        eprint("ERROR: call_dgeqp3 failed in qr (%s:%d): %d", __FILE__, __LINE__, info);
        return info;
    }

    info = call_dorgqr(m, m, k, q.data(), m, tau.get());

    r.triangularView<Eigen::StrictlyLower>().setZero();

    if (info != 0) {
        eprint("ERROR: call_dorgqr failed in qr (%s:%d): %d", __FILE__, __LINE__, info);
        return info;
    }

    for (isize_t i = 0; i < n; ++i)
        p(i) = static_cast<int>(jpvt[i]);

    return 0;
}

int qr_eco(const Eigen::Ref<const mat> &a, mat &q, mat &r, imat &p)
{
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t k = m < n ? m : n;

    if (m <= n)
        return qr(a, q, r, p);

    p.setZero(n, 1);

    q = a;
    std::unique_ptr<bint[]> jpvt(new bint[n]);
    std::fill_n(jpvt.get(), n, 0);
    std::unique_ptr<double[]> tau(new double[k]);

    int info = call_dgeqp3(m, n, q.data(), m, jpvt.get(), tau.get());

    r = q.topRows(n).triangularView<Eigen::Upper>();

    if (info != 0) {
        eprint("ERROR: call_dgeqp3 failed in qr_eco (%s:%d): %d", __FILE__, __LINE__, info);
        return info;
    }

    for (isize_t i = 0; i < n; ++i)
        p(i) = static_cast<int>(jpvt[i]);

    info = call_dorgqr(m, n, k, q.data(), m, tau.get());

    if (info != 0)
        eprint("ERROR: call_dorgqr failed in qr_eco (%s:%d): %d", __FILE__, __LINE__, info);

    return info;
}

int gesvd(const Eigen::Ref<const mat> &a, mat &s)
{
    mat aa = a;
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    s.resize(m < n ? m : n, 1);

    int info = call_dgesvd('N', 'N', m, n, aa.data(), m, s.data(), nullptr, m, nullptr, n);

    if (info != 0)
        eprint("ERROR: call_dgesvd failed in gesvd (%s:%d): %d", __FILE__, __LINE__, info);

    return info;
}

int gesvd(const Eigen::Ref<const mat> &a, mat &u, mat &s, mat &vt)
{
    mat aa = a;
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t k = m < n ? m : n;

    u.resize(m, m);
    s.setZero(m, n);
    vt.resize(n, n);

    std::unique_ptr<double[]> ss(new double[k]);

    int info = call_dgesvd('A', 'A', m, n, aa.data(), m, ss.get(), u.data(), m, vt.data(), n);

    if (info != 0)
        eprint("ERROR: call_dgesvd failed in gesvd (%s:%d): %d", __FILE__, __LINE__, info);

    for (isize_t i = 0; i < k; ++i)
        s(i, i) = ss[i];

    return info;
}

int gesvd_eco(const Eigen::Ref<const mat> &a, mat &u, mat &s, mat &vt)
{
    mat aa = a;
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t k = m < n ? m : n;

    u.resize(m, k);
    s.setZero(k, k);
    vt.resize(k, n);

    int info = call_dgesvd('S', 'S', m, n, aa.data(), m, s.data(), u.data(), m, vt.data(), k);

    if (info != 0)
        eprint("ERROR: call_dgesvd failed in gesvd (%s:%d): %d", __FILE__, __LINE__, info);

    for (isize_t i = 1; i < k; ++i) {
        s(i, i) = s(i, 0);
        s(i, 0) = 0.0;
    }

    return info;
}

int gesdd(const Eigen::Ref<const mat> &a, mat &s)
{
    mat aa = a;
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    s.resize(m < n ? m : n, 1);

    int info = call_dgesdd('N', m, n, aa.data(), m, s.data(), nullptr, m, nullptr, n);

    if (info != 0)
        eprint("ERROR: call_dgesdd failed in gesdd (%s:%d): %d", __FILE__, __LINE__, info);

    return info;
}

int gesdd(const Eigen::Ref<const mat> &a, mat &u, mat &s, mat &vt)
{
    mat aa = a;
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t k = m < n ? m : n;

    u.resize(m, m);
    s.setZero(m, n);
    vt.resize(n, n);

    std::unique_ptr<double[]> ss(new double[k]);

    int info = call_dgesdd('A', m, n, aa.data(), m, ss.get(), u.data(), m, vt.data(), n);

    if (info != 0)
        eprint("ERROR: call_dgesdd failed in gesdd (%s:%d): %d", __FILE__, __LINE__, info);

    for (isize_t i = 0; i < k; ++i)
        s(i, i) = ss[i];

    return info;
}

int gesdd_eco(const Eigen::Ref<const mat> &a, mat &u, mat &s, mat &vt)
{
    mat aa = a;
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t k = m < n ? m : n;

    u.resize(m, k);
    s.setZero(k, k);
    vt.resize(k, n);

    int info = call_dgesdd('S', m, n, aa.data(), m, s.data(), u.data(), m, vt.data(), k);

    if (info != 0)
        eprint("ERROR: call_dgesdd failed in gesdd_eco (%s:%d): %d", __FILE__, __LINE__, info);

    for (isize_t i = 1; i < k; ++i) {
        s(i, i) = s(i, 0);
        s(i, 0) = 0.0;
    }

    return info;
}

int gels(const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b, mat &x)
{
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t nrhs = ncol(b);
    mat aa = a;

    mat bb;
    if (m < n) {
        bb.resize(n, nrhs);
        bb.topRows(m) = b;
    }
    else
        bb = b;

    isize_t lda = nrow(aa);
    isize_t ldb = nrow(bb);

    int info = call_dgels('N', m, n, nrhs, aa.data(), lda, bb.data(), ldb);

    if (info != 0)
        eprint("ERROR: call_dgels failed in gels (%s:%d): %d", __FILE__, __LINE__, info);

    if (info == 0 && m >= n) {
        mat d = aa.diagonal().cwiseAbs();
        double tol = m * d.maxCoeff() * std::numeric_limits<double>::epsilon();
        if ((d.array() < tol).any()) {
            eprint("ERROR: rank deficient in gels, tol = %g\n", tol);
            info = -999;
        }
    }

    if (m <= n)
        x = bb;
    else
        x = bb.topRows(n);

    return info;
}

int gelsy(double rcond, const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b, mat &x)
{
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t nrhs = ncol(b);
    mat aa = a;

    mat bb;
    if (m < n) {
        bb.resize(n, nrhs);
        bb.topRows(m) = b;
    }
    else
        bb = b;

    isize_t lda = nrow(aa);
    isize_t ldb = nrow(bb);

    bint rank = 0;
    std::unique_ptr<bint[]> jpvt(new bint[n]);
    std::fill_n(jpvt.get(), n, 0);

    int info = call_dgelsy(m, n, nrhs, aa.data(), lda, bb.data(), ldb, jpvt.get(), rcond, &rank);

    if (info != 0)
        eprint("ERROR: call_dgelsy failed in gelsy (%s:%d): %d", __FILE__, __LINE__, info);

    if (m <= n)
        x = bb;
    else
        x = bb.topRows(n);

    return info == 0 ? rank : -1;
}

int gelsd(double rcond, const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b, mat &x)
{
    isize_t m = nrow(a);
    isize_t n = ncol(a);
    isize_t nrhs = ncol(b);
    mat aa = a;

    mat bb;
    if (m < n) {
        bb.resize(n, nrhs);
        bb.topRows(m) = b;
    }
    else
        bb = b;

    isize_t lda = nrow(aa);
    isize_t ldb = nrow(bb);

    bint rank = 0;
    std::unique_ptr<double[]> s(new double[m > n ? m : n]);

    int info = call_dgelsd(m, n, nrhs, aa.data(), lda, bb.data(), ldb, s.get(), rcond, &rank);

    if (info != 0)
        eprint("ERROR: call_dgelsd failed in gelsd (%s:%d): %d", __FILE__, __LINE__, info);

    if (m <= n)
        x = bb;
    else
        x = bb.topRows(n);

    return info == 0 ? rank : -1;
}
