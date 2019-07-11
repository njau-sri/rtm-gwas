#include <memory>
#include "lapack.h"


double C_dnrm2(bint n, const double *x, bint incx)
{
    return dnrm2_(&n, x, &incx);
}

double C_ddot(bint n, const double *x, bint incx, const double *y, bint incy)
{
    return ddot_(&n, x, &incx, y, &incy);
}

int C_dscal(bint n, double alpha, double *x, bint incx)
{
    dscal_(&n, &alpha, x, &incx);
    return 0;
}

int C_daxpy(bint n, double alpha, const double *x, bint incx, double *y, bint incy)
{
    daxpy_(&n, &alpha, x, &incx, y, &incy);
    return 0;
}

int C_dgemv(char trans, bint m, bint n, double alpha, const double *a, bint lda,
            const double *x, bint incx, double beta, double *y, bint incy)
{
    dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
    return 0;
}

int C_dgemm(char transa, char transb, bint m, bint n, bint k, double alpha, const double *a, bint lda,
            const double *b, bint ldb, double beta, double *c, bint ldc)
{
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    return 0;
}

int C_dsyrk(char uplo, char trans, bint n, bint k, double alpha, const double *a, bint lda,
            double beta, double *c, bint ldc)
{
    dsyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
    return 0;
}

int C_dtrsv(char uplo, char trans, char diag, bint n, const double *a, bint lda, double *x, bint incx)
{
    dtrsv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
    return 0;
}

int C_dtrsm(char side, char uplo, char transa, char diag, bint m, bint n, double alpha,
            const double *a, bint lda, double *b, bint ldb)
{
    dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
    return 0;
}

int C_dgels(char trans, bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        std::unique_ptr<double[]> work(new double[lwork]);
        dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work.get(), &lwork, &info);
    }

    return info;
}

int C_dgelsy(bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb, bint *jpvt, double rcond, bint *rank)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        std::unique_ptr<double[]> work(new double[lwork]);
        dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, work.get(), &lwork, &info);
    }

    return info;
}

int C_dgeqp3(bint m, bint n, double *a, bint lda, bint *jpvt, double *tau)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dgeqp3_(&m, &n, a, &lda, jpvt, tau, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        std::unique_ptr<double[]> work(new double[lwork]);
        dgeqp3_(&m, &n, a, &lda, jpvt, tau, work.get(), &lwork, &info);
    }

    return info;
}

int C_dorgqr(bint m, bint n, bint k, double *a, bint lda, const double *tau)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dorgqr_(&m, &n, &k, a, &lda, tau, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        std::unique_ptr<double[]> work(new double[lwork]);
        dorgqr_(&m, &n, &k, a, &lda, tau, work.get(), &lwork, &info);
    }

    return info;
}

int C_dpotrf(char uplo, bint n, double *a, bint lda)
{
    bint info = 0;

    dpotrf_(&uplo, &n, a, &lda, &info);

    return info;
}

int C_dsyevr(char jobz, char range, char uplo, bint n, double *a, bint lda, double vl, double vu, bint il, bint iu,
             double abstol, bint *m, double *w, double *z, bint ldz, bint *isuppz)
{
    bint info = 0;

    double wkopt;
    bint iwkopt;
    bint lwork = -1, liwork = -1;

    dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz,
            &wkopt, &lwork, &iwkopt, &liwork, &info);

    if (info != 0)
        return info;

    lwork = static_cast<bint>(wkopt);
    liwork = iwkopt;

    std::unique_ptr<bint[]> iwork(new bint[liwork]);
    std::unique_ptr<double[]> work(new double[lwork]);

    dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz,
            work.get(), &lwork, iwork.get(), &liwork, &info);

    return info;
}

int C_dsytrf(char uplo, bint n, double *a, bint lda, bint *ipiv)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dsytrf_(&uplo, &n, a, &lda, ipiv, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        std::unique_ptr<double[]> work(new double[lwork]);
        dsytrf_(&uplo, &n, a, &lda, ipiv, work.get(), &lwork, &info);
    }

    return info;
}

int C_dsytri(char uplo, bint n, double *a, bint lda, const bint *ipiv)
{
    bint info = 0;

    std::unique_ptr<double[]> work(new double[n]);
    dsytri_(&uplo, &n, a, &lda, ipiv, work.get(), &info);

    return info;
}

int C_dsycon(char uplo, bint n, const double *a, bint lda, const bint *ipiv, double anorm, double *rcond)
{
    bint info = 0;

    std::unique_ptr<bint[]> iwork(new bint[n]);
    std::unique_ptr<double[]> work(new double[n*2]);

    dsycon_(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work.get(), iwork.get(), &info);

    return info;
}

double C_dlansy(char norm, char uplo, bint n, const double *a, bint lda)
{
    std::unique_ptr<double[]> work(new double[n]);
    return dlansy_(&norm, &uplo, &n, a, &lda, work.get());
}

int M_dsyinv(int n, double *a)
{
    if (n < 1)
        return 1;

    if (n == 1) {
        a[0] = 1.0 / a[0];
        return 0;
    }

    char uplo = 'U';
    char norm = 'I';

    bint m = n;
    bint lda = n;
    bint info = 0;

    std::unique_ptr<bint[]> iwork(new bint[m*2]);
    auto ipiv = iwork.get() + m;

    double wkopt;
    bint lwork = -1;

    dsytrf_(&uplo, &m, a, &lda, ipiv, &wkopt, &lwork, &info);

    if (info != 0)
        return 2;

    lwork = static_cast<bint>(wkopt);
    std::unique_ptr<double[]> work(new double[lwork + m * 2]);

    double anorm = dlansy_(&norm, &uplo, &m, a, &lda, work.get());

    dsytrf_(&uplo, &m, a, &lda, ipiv, work.get(), &lwork, &info);

    if (info != 0)
        return 3;

    double rcond = 0.0;

    dsycon_(&uplo, &m, a, &lda, ipiv, &anorm, &rcond, work.get(), iwork.get(), &info);

    if (info != 0 || rcond < 1e-8)
        return 4;

    dsytri_(&uplo, &m, a, &lda, ipiv, work.get(), &info);

    if (info != 0)
        return 5;

    for (bint i = 1; i < m; ++i) {
        for (bint j = 0; j < i; ++j)
            a[j*m+i] = a[i*m+j];
    }

    return 0;
}
