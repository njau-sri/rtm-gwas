#include "lapack.h"

#ifdef BLAS_INT64
    static_assert(sizeof(bint) == 8, "BLAS integer size is not 8");
#else
    static_assert(sizeof(bint) == 4, "BLAS integer size is not 4");
#endif

extern "C"
{

double dnrm2_(bint *n, const double *x, bint *incx);

double ddot_(bint *n, const double *x, bint *incx, const double *y, bint *incy);

void dscal_(bint *n, double *alpha, double *x, bint *incx);

void daxpy_(bint *n, double *alpha, const double *x, bint *incx, double *y, bint *incy);

void dgemv_(char *trans, bint *m, bint *n, double *alpha, const double *a, bint *lda,
            const double *x, bint *incx, double *beta, double *y, bint *incy);

void dgemm_(char *transa, char *transb, bint *m, bint *n, bint *k, double *alpha, const double *a, bint *lda,
            const double *b, bint *ldb, double *beta, double *c, bint *ldc);

void dsyrk_(char *uplo, char *trans, bint *n, bint *k, double *alpha, const double *a, bint *lda,
            double *beta, double *c, bint *ldc);

void dgels_(char *trans, bint *m, bint *n, bint *nrhs, double *a, bint *lda, double *b, bint *ldb,
            double *work, bint *lwork, bint *info);

void dgelsy_(bint *m, bint *n, bint *nrhs, double *a, bint *lda, double *b, bint *ldb, bint *jpvt,
             double *rcond, bint *rank, double *work, bint *lwork, bint *info);

void dsyevr_(char *jobz, char *range, char *uplo, bint *n, double *a, bint *lda, double *vl, double *vu,
             bint *il, bint *iu, double *abstol, bint *m, double *w, double *z, bint *ldz, bint *isuppz,
             double *work, bint *lwork, bint *iwork, bint *liwork, bint *info);

void dtrtrs_(char *uplo, char *trans, char *diag, bint *n, bint *nrhs, const double *a, bint *lda,
             double *b, bint *ldb, bint *info);

void dpotrf_(char *uplo, bint *n, double *a, bint *lda, bint *info);

void dgeqrf_(bint *m, bint *n, double *a, bint* lda, double *tau, double *work, bint *lwork, bint *info);

void dgeqp3_(bint *m, bint *n, double *a, bint *lda, bint *jpvt, double *tau, double *work, bint *lwork, bint *info);

void dormqr_(char *side, char *trans, bint *m, bint *n, bint *k, const double *a, bint *lda, const double *tau,
             double *c, bint *ldc, double *work, bint *lwork, bint *info);

void dorgqr_(bint *m, bint *n, bint *k, double *a, bint *lda, const double *tau, double *work, bint *lwork, bint *info);

void dsytrf_(char *uplo, bint *n, double *a, bint *lda, bint* ipiv, double *work, bint *lwork, bint *info);

void dsytri_(char *uplo, bint *n, double *a, bint *lda, const bint *ipiv, double *work, bint *info);

} // extern "C"

double la_dnrm2(bint n, const double *x, bint incx)
{
    return dnrm2_(&n, x, &incx);
}

double la_ddot(bint n, const double *x, bint incx, const double *y, bint incy)
{
    return ddot_(&n, x, &incx, y, &incy);
}

int la_dscal(bint n, double alpha, double *x, bint incx)
{
    dscal_(&n, &alpha, x, &incx);
    return 0;
}

int la_daxpy(bint n, double alpha, const double *x, bint incx, double *y, bint incy)
{
    daxpy_(&n, &alpha, x, &incx, y, &incy);
    return 0;
}

int la_dgemv(char trans, bint m, bint n, double alpha, const double *a, bint lda, const double *x, bint incx,
             double beta, double *y, bint incy)
{
    dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
    return 0;
}

int la_dgemm(char transa, char transb, bint m, bint n, bint k, double alpha, const double *a, bint lda,
             const double *b, bint ldb, double beta, double *c, bint ldc)
{
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    return 0;
}

int la_dsyrk(char uplo, char trans, bint n, bint k, double alpha, const double *a, bint lda,
             double beta, double *c, bint ldc)
{
    dsyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
    return 0;
}

int la_dgels(char trans, bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        double *work = new double[lwork];
        dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
        delete[] work;
    }

    return info;
}

int la_dgelsy(bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb,
              bint *jpvt, double rcond, bint *rank)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        double *work = new double[lwork];
        dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, work, &lwork, &info);
        delete[] work;
    }

    return info;
}

int la_dsyevr(char jobz, char range, char uplo, bint n, double *a, bint lda, double vl, double vu, bint il, bint iu,
              double abstol, bint *m, double *w, double *z, bint ldz, bint *isuppz)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1, iwkopt, liwork = -1;
    dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz,
            &wkopt, &lwork, &iwkopt, &liwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        liwork = iwkopt;
        bint *iwork = new bint[liwork];
        double *work = new double[lwork];
        dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz,
                work, &lwork, iwork, &liwork, &info);
        delete[] work;
        delete[] iwork;
    }

    return info;
}

int la_dtrtrs(char uplo, char trans, char diag, bint n, bint nrhs, const double *a, bint lda, double *b, bint ldb)
{
    bint info = 0;
    dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

int la_dpotrf(char uplo, bint n, double *a, bint lda)
{
    bint info = 0;
    dpotrf_(&uplo, &n, a, &lda, &info);
    return info;
}

int la_dgeqrf(bint m, bint n, double *a, bint lda, double *tau)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dgeqrf_(&m, &n, a, &lda, tau, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        double *work = new double[lwork];
        dgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
        delete[] work;
    }

    return info;
}

int la_dgeqp3(bint m, bint n, double *a, bint lda, bint *jpvt, double *tau)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dgeqp3_(&m, &n, a, &lda, jpvt, tau, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        double *work = new double[lwork];
        dgeqp3_(&m, &n, a, &lda, jpvt, tau, work, &lwork, &info);
        delete[] work;
    }

    return info;
}

int la_dormqr(char side, char trans, bint m, bint n, bint k, const double *a, bint lda,
              const double *tau, double *c, bint ldc)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dormqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        double *work = new double[lwork];
        dormqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
        delete[] work;
    }

    return info;
}

int la_dorgqr(bint m, bint n, bint k, double *a, bint lda, const double *tau)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dorgqr_(&m, &n, &k, a, &lda, tau, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        double *work = new double[lwork];
        dorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
        delete[] work;
    }

    return info;
}

int la_dsytrf(char uplo, bint n, double *a, bint lda, bint *ipiv)
{
    bint info = 0;

    double wkopt;
    bint lwork = -1;
    dsytrf_(&uplo, &n, a, &lda, ipiv, &wkopt, &lwork, &info);

    if (info == 0) {
        lwork = static_cast<bint>(wkopt);
        double *work = new double[lwork];
        dsytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
        delete[] work;
    }

    return info;
}

int la_dsytri(char uplo, bint n, double *a, bint lda, const bint *ipiv)
{
    bint info = 0;

    double *work = new double[n];
    dsytri_(&uplo, &n, a, &lda, ipiv, work, &info);
    delete[] work;

    return info;
}
