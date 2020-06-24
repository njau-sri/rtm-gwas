#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "lapack.h"


#ifdef BLASINT64
static_assert(sizeof(bint) == 8, "blas integer size is not 8");
#else
static_assert(sizeof(bint) == 4, "blas integer size is not 4");
#endif


double dnrm2_(bint *n, const double *x, bint *incx);

double ddot_(bint *n, const double *x, bint *incx, const double *y, bint *incy);

void dscal_(bint *n, double *alpha, double *x, bint *incx);

void daxpy_(bint *n, double *alpha, const double *x, bint *incx, double *y, bint *incy);

void dgemv_(char *trans, bint *m, bint *n, double *alpha, const double *a, bint *lda,
            const double *x, bint *incx, double *beta, double *y, bint *incy);

void dgemm_(char *transa, char *transb, bint *m, bint *n, bint *k, double *alpha, const double *a,
            bint *lda, const double *b, bint *ldb, double *beta, double *c, bint *ldc);

void dsyrk_(char *uplo, char *trans, bint *n, bint *k, double *alpha, const double *a, bint *lda,
            double *beta, double *c, bint *ldc);

void dtrsv_(char *uplo, char *trans, char *diag, bint *n, const double *a, bint *lda,
            double *x, bint *incx);

void dtrsm_(char *side, char *uplo, char *transa, char *diag, bint *m, bint *n, double *alpha,
            const double *a, bint *lda, double *b, bint *ldb);

void dgels_(char *trans, bint *m, bint *n, bint *nrhs, double *a, bint *lda,
            double *b, bint *ldb, double *work, bint *lwork, bint *info);

void dgelsy_(bint *m, bint *n, bint *nrhs, double *a, bint *lda, double *b, bint *ldb, bint *jpvt,
             double *rcond, bint *rank, double *work, bint *lwork, bint *info);

void dgelsd_(bint *m, bint *n, bint *nrhs, double *a, bint *lda, double *b, bint *ldb, double *s,
             double *rcond, bint *rank, double *work, bint *lwork, bint *iwork, bint *info);

void dgeqrf_(bint *m, bint *n, double *a, bint *lda, double *tau,
             double *work, bint *lwork, bint *info);

void dgeqp3_(bint *m, bint *n, double *a, bint *lda, bint *jpvt, double *tau,
             double *work, bint *lwork, bint *info);

void dorgqr_(bint *m, bint *n, bint *k, double *a, bint *lda, const double *tau,
             double *work, bint *lwork, bint *info);

void dpotrf_(char *uplo, bint *n, double *a, bint *lda, bint *info);

void dsyev_(char *jobz, char *uplo, bint *n, double *a, bint *lda, double *w, double *work, bint *lwork, bint *info);

void dsyevr_(char *jobz, char *range, char *uplo, bint *n, double *a, bint *lda, double *vl,
             double *vu, bint *il, bint *iu, double *abstol, bint *m, double *w, double *z,
             bint *ldz, bint *isuppz, double *work, bint *lwork, bint *iwork,
             bint *liwork, bint *info);

void dsytrf_(char *uplo, bint *n, double *a, bint *lda, bint *ipiv,
             double *work, bint *lwork, bint *info);

void dsytri_(char *uplo, bint *n, double *a, bint *lda, const bint *ipiv,
             double *work, bint *info);

void dsycon_(char *uplo, bint *n, const double *a, bint *lda, const bint *ipiv, double *anorm,
             double *rcond, double *work, bint *iwork, bint *info);

double dlansy_(char *norm, char *uplo, bint *n, const double *a, bint *lda, double *work);

double dlange_(char *norm, bint *m, bint *n, const double *a, bint *lda, double *work);

void dgesvd_(char *jobu, char *jobvt, bint *m, bint *n, double *a, bint *lda, double *s, double *u, bint *ldu,
             double *vt, bint *ldvt, double *work, bint *lwork, bint *info);

void dgesdd_(char *jobz, bint *m, bint *n, double *a, bint *lda, double *s, double *u, bint *ldu,
             double *vt, bint *ldvt, double *work, bint *lwork, bint *iwork, bint *info);


double call_dnrm2(bint n, const double *x, bint incx)
{
    return dnrm2_(&n, x, &incx);
}

double call_ddot(bint n, const double *x, bint incx, const double *y, bint incy)
{
    return ddot_(&n, x, &incx, y, &incy);
}

int call_dscal(bint n, double alpha, double *x, bint incx)
{
    dscal_(&n, &alpha, x, &incx);
    return 0;
}

int call_daxpy(bint n, double alpha, const double *x, bint incx, double *y, bint incy)
{
    daxpy_(&n, &alpha, x, &incx, y, &incy);
    return 0;
}

int call_dgemv(char trans, bint m, bint n, double alpha, const double *a, bint lda,
               const double *x, bint incx, double beta, double *y, bint incy)
{
    dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
    return 0;
}

int call_dgemm(char transa, char transb, bint m, bint n, bint k, double alpha, const double *a,
               bint lda, const double *b, bint ldb, double beta, double *c, bint ldc)
{
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    return 0;
}

int call_dsyrk(char uplo, char trans, bint n, bint k, double alpha, const double *a, bint lda,
               double beta, double *c, bint ldc)
{
    dsyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
    return 0;
}

int call_dtrsv(char uplo, char trans, char diag, bint n, const double *a, bint lda,
               double *x, bint incx)
{
    dtrsv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
    return 0;
}

int call_dtrsm(char side, char uplo, char transa, char diag, bint m, bint n, double alpha,
               const double *a, bint lda, double *b, bint ldb)
{
    dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
    return 0;
}

int call_dgels(char trans, bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb)
{
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    double *work = NULL;

    dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, &wkopt, &lwork, &info);
    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dgels\n", stderr);
        return -999;
    }

    dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    free(work);

    return info;
}

int call_dgelsy(bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb,
                bint *jpvt, double rcond, bint *rank)
{
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    double *work = NULL;

    dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, &wkopt, &lwork, &info);
    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dgelsy\n", stderr);
        return -999;
    }

    dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, work, &lwork, &info);
    free(work);

    return info;
}

int call_dgelsd(bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb,
                double *s, double rcond, bint *rank)
{
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    bint liwork;
    double *work = NULL;
    bint *iwork = NULL;

    dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, &wkopt, &lwork, &liwork, &info);
    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dgelsd\n", stderr);
        return -999;
    }

    iwork = (bint*) malloc(sizeof(bint) * liwork);
    if (iwork == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dgelsd\n", stderr);
        free(work);
        return -999;
    }

    dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, iwork, &info);
    free(iwork);
    free(work);

    return info;
}

int call_dgeqrf(bint m, bint n, double *a, bint lda, double *tau)
{
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    double *work = NULL;

    dgeqrf_(&m, &n, a, &lda, tau, &wkopt, &lwork, &info);
    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dgeqrf\n", stderr);
        return -999;
    }

    dgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
    free(work);

    return info;
}

int call_dgeqp3(bint m, bint n, double *a, bint lda, bint *jpvt, double *tau)
{
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    double *work = NULL;

    dgeqp3_(&m, &n, a, &lda, jpvt, tau, &wkopt, &lwork, &info);
    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dgeqp3\n", stderr);
        return -999;
    }

    dgeqp3_(&m, &n, a, &lda, jpvt, tau, work, &lwork, &info);
    free(work);

    return info;
}

int call_dorgqr(bint m, bint n, bint k, double *a, bint lda, const double *tau)
{
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    double *work = NULL;

    dorgqr_(&m, &n, &k, a, &lda, tau, &wkopt, &lwork, &info);
    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dorgqr\n", stderr);
        return -999;
    }

    dorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    free(work);

    return info;
}

int call_dpotrf(char uplo, bint n, double *a, bint lda)
{
    bint info = 0;

    dpotrf_(&uplo, &n, a, &lda, &info);

    return info;
}

int call_dsyev(char jobz, char uplo, bint n, double *a, bint lda, double *w)
{
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    double *work = NULL;

    dsyev_(&jobz, &uplo, &n, a, &lda, w, &wkopt, &lwork, &info);

    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dsyev\n", stderr);
        return -999;
    }

    dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

    free(work);

    return 0;
}

int call_dsyevr(char jobz, char range, char uplo, bint n, double *a, bint lda,
                double vl, double vu, bint il, bint iu, double abstol, bint *m,
                double *w, double *z, bint ldz, bint *isuppz)
{
    bint info = 0;
    double wkopt;
    bint iwkopt;
    bint lwork = -1;
    bint liwork = -1;
    bint *iwork = NULL;
    double *work = NULL;

    dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz,
            &wkopt, &lwork, &iwkopt, &liwork, &info);

    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    liwork = iwkopt;

    iwork = (bint*) malloc(sizeof(bint) * liwork);
    if (iwork == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dsyevr\n", stderr);
        return -999;
    }

    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        free(iwork);
        fputs("ERROR: memory allocation (malloc) failed in call_dsyevr\n", stderr);
        return -999;
    }

    dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz,
            work, &lwork, iwork, &liwork, &info);

    free(iwork);
    free(work);

    return info;
}

int call_dsytrf(char uplo, bint n, double *a, bint lda, bint *ipiv)
{
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    double *work = NULL;

    dsytrf_(&uplo, &n, a, &lda, ipiv, &wkopt, &lwork, &info);
    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dsytrf\n", stderr);
        return -999;
    }

    dsytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
    free(work);

    return info;
}

int call_dsytri(char uplo, bint n, double *a, bint lda, const bint *ipiv)
{
    bint info = 0;
    double *work = NULL;

    work = (double*) malloc(sizeof(double) * n);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dsytri\n", stderr);
        return -999;
    }

    dsytri_(&uplo, &n, a, &lda, ipiv, work, &info);
    free(work);

    return info;
}

int call_dsycon(char uplo, bint n, const double *a, bint lda, const bint *ipiv,
                double anorm, double *rcond)
{
    bint info = 0;
    bint *iwork = NULL;
    double *work = NULL;

    iwork = (bint*) malloc(sizeof(bint) * n);
    if (iwork == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dsycon\n", stderr);
        return -999;
    }

    work = (double*) malloc(sizeof(double) * n * 2);
    if (work == NULL) {
        free(iwork);
        fputs("ERROR: memory allocation (malloc) failed in call_dsycon\n", stderr);
        return -999;
    }

    dsycon_(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, iwork, &info);
    free(iwork);
    free(work);

    return info;
}

double call_dlansy(char norm, char uplo, bint n, const double *a, bint lda)
{
    double res = 0.0;
    double *work = NULL;

    work = (double*) malloc(sizeof(double) * n);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dlansy\n", stderr);
        return res;
    }

    res = dlansy_(&norm, &uplo, &n, a, &lda, work);
    free(work);

    return res;
}

double call_dlange(char norm, bint m, bint n, const double *a, bint lda)
{
    double res = 0.0;
    double *work = NULL;

    if (norm == 'I' || norm == 'i') {
        work = (double*) malloc(sizeof(double) * m);
        if (work == NULL) {
            fputs("ERROR: memory allocation (malloc) failed in call_dlange\n", stderr);
            return res;
        }
    }

    res = dlange_(&norm, &m, &n, a, &lda, work);

    if (work)
        free(work);

    return res;
}

int call_dgesvd(char jobu, char jobvt, bint m, bint n, double *a, bint lda, double *s,
                double *u, bint ldu, double *vt, bint ldvt)
{
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    double *work = NULL;

    dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info);
    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dgesvd\n", stderr);
        return -999;
    }

    dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
    free(work);

    return info;
}

int call_dgesdd(char jobz, bint m, bint n, double *a, bint lda, double *s,
                double *u, bint ldu, double *vt, bint ldvt)
{
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    double *work = NULL;
    bint *iwork = NULL;

    dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, iwork, &info);
    if (info != 0)
        return info;

    lwork = (bint) wkopt;
    work = (double*) malloc(sizeof(double) * lwork);
    if (work == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in call_dgesdd\n", stderr);
        return -999;
    }

    iwork = (bint*) malloc(sizeof(bint) * 8 * (m > n ? n : m));

    dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
    free(iwork);
    free(work);

    return info;
}

int dsyinv(bint n, double *a)
{
    char uplo = 'U';
    char norm = '1';

    bint lda = n;
    bint info = 0;
    double wkopt;
    bint lwork = -1;
    bint lwork2 = 0;
    bint *ipiv = NULL;
    bint *iwork = NULL;
    double *work = NULL;

    bint i = 0;
    bint j = 0;
    double rcond = 0.0;
    double anorm = 0.0;
    double tol = DBL_EPSILON * n * 100;

    if (n < 1)
        return -999;

    if (n == 1) {
        a[0] = 1.0 / a[0];
        return 0;
    }

    ipiv = (bint*) malloc(sizeof(bint) * n*2);
    if (ipiv == NULL) {
        fputs("ERROR: memory allocation (malloc) failed in dsyinv\n", stderr);
        return -999;
    }

    dsytrf_(&uplo, &n, a, &lda, ipiv, &wkopt, &lwork, &info);
    if (info != 0) {
        free(ipiv);
        return info;
    }

    lwork = lwork2 = (bint) wkopt;
    if (lwork2 < n*2)
        lwork2 = n*2;

    work = (double*) malloc(sizeof(double) * lwork2);
    if (work == NULL) {
        free(ipiv);
        fputs("ERROR: memory allocation (malloc) failed in dsyinv\n", stderr);
        return -999;
    }

    anorm = dlansy_(&norm, &uplo, &n, a, &lda, work);

    dsytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
    if (info != 0) {
        free(ipiv);
        free(work);
        return info;
    }

    iwork = ipiv + n;
    dsycon_(&uplo, &n, a, &lda, ipiv, &anorm, &rcond, work, iwork, &info);

    if (info == 0 && rcond < tol) {
        fprintf(stderr, "ERROR: matrix is close to singular in dsyinv, rcond = %g, tol = %g\n", rcond, tol);
        info = -999;
    }

    if (info == 0)
        dsytri_(&uplo, &n, a, &lda, ipiv, work, &info);
    
    free(ipiv);
    free(work);

    if (info != 0)
        return info;

    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j)
            a[j*n+i] = a[i*n+j];
    }

    return 0;
}
