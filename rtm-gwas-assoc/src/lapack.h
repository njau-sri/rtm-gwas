#ifndef LAPACK_H
#define LAPACK_H


#ifdef BLASINT64
#include <stddef.h>
#define  bint  ptrdiff_t
#else
#define  bint  int
#endif


#ifdef __cplusplus
extern "C" {
#endif // __cplusplus


    double call_dnrm2(bint n, const double *x, bint incx);

    double call_ddot(bint n, const double *x, bint incx, const double *y, bint incy);

    int call_dscal(bint n, double alpha, double *x, bint incx);

    int call_daxpy(bint n, double alpha, const double *x, bint incx, double *y, bint incy);

    int call_dgemv(char trans, bint m, bint n, double alpha, const double *a, bint lda,
                   const double *x, bint incx, double beta, double *y, bint incy);

    int call_dgemm(char transa, char transb, bint m, bint n, bint k, double alpha, const double *a,
                   bint lda, const double *b, bint ldb, double beta, double *c, bint ldc);

    int call_dsyrk(char uplo, char trans, bint n, bint k, double alpha, const double *a, bint lda,
                   double beta, double *c, bint ldc);

    int call_dtrsv(char uplo, char trans, char diag, bint n, const double *a, bint lda,
                   double *x, bint incx);

    int call_dtrsm(char side, char uplo, char transa, char diag, bint m, bint n, double alpha,
                   const double *a, bint lda, double *b, bint ldb);

    int call_dgels(char trans, bint m, bint n, bint nrhs, double *a, bint lda,
                   double *b, bint ldb);

    int call_dgelsy(bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb,
                    bint *jpvt, double rcond, bint *rank);

    int call_dgelsd(bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb,
                    double *s, double rcond, bint *rank);

    int call_dgeqrf(bint m, bint n, double *a, bint lda, double *tau);

    int call_dgeqp3(bint m, bint n, double *a, bint lda, bint *jpvt, double *tau);

    int call_dorgqr(bint m, bint n, bint k, double *a, bint lda, const double *tau);

    int call_dpotrf(char uplo, bint n, double *a, bint lda);

    int call_dsyev(char jobz, char uplo, bint n, double *a, bint lda, double *w);

    int call_dsyevr(char jobz, char range, char uplo, bint n, double *a, bint lda, double vl,
                    double vu, bint il, bint iu, double abstol, bint *m, double *w,
                    double *z, bint ldz, bint *isuppz);

    int call_dsytrf(char uplo, bint n, double *a, bint lda, bint *ipiv);

    int call_dsytri(char uplo, bint n, double *a, bint lda, const bint *ipiv);

    int call_dsycon(char uplo, bint n, const double *a, bint lda, const bint *ipiv,
                    double anorm, double *rcond);

    double call_dlansy(char norm, char uplo, bint n, const double *a, bint lda);

    double call_dlange(char norm, bint m, bint n, const double *a, bint lda);

    int call_dgesvd(char jobu, char jobvt, bint m, bint n, double *a, bint lda, double *s,
                    double *u, bint ldu, double *vt, bint ldvt);

    int call_dgesdd(char jobz, bint m, bint n, double *a, bint lda, double *s,
                    double *u, bint ldu, double *vt, bint ldvt);

    // Symmetric Matrix Inverse
    //
    //   Input
    //     n   the order of the symmetric matrix A
    //     a   the n*n symmetric matrix A
    //
    int dsyinv(bint n, double *a);


#ifdef __cplusplus
}
#endif // __cplusplus


#endif // LAPACK_H
