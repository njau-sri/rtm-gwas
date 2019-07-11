#ifndef LAPACK_H
#define LAPACK_H


#ifdef BLAS64
#include <cstddef>
using bint = std::ptrdiff_t;
static_assert(sizeof(bint) == 8, "BLAS integer size is not 8");
#else
using bint = int;
static_assert(sizeof(bint) == 4, "BLAS integer size is not 4");
#endif


#ifdef __cplusplus
extern "C" {
#endif // __cplusplus


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

void dtrsv_(char *uplo, char *trans, char *diag, bint *n, const double *a, bint *lda, double *x, bint *incx);

void dtrsm_(char *side, char *uplo, char *transa, char *diag, bint *m, bint *n, double *alpha,
            const double *a, bint *lda, double *b, bint *ldb);

void dgels_(char *trans, bint *m, bint *n, bint *nrhs, double *a, bint *lda, double *b, bint *ldb,
            double *work, bint *lwork, bint *info);

void dgelsy_(bint *m, bint *n, bint *nrhs, double *a, bint *lda, double *b, bint *ldb, bint *jpvt,
             double *rcond, bint *rank, double *work, bint *lwork, bint *info);

void dgeqp3_(bint *m, bint *n, double *a, bint *lda, bint *jpvt, double *tau, double *work, bint *lwork, bint *info);

void dorgqr_(bint *m, bint *n, bint *k, double *a, bint *lda, const double *tau, double *work, bint *lwork, bint *info);

void dpotrf_(char *uplo, bint *n, double *a, bint *lda, bint *info);

void dsyevr_(char *jobz, char *range, char *uplo, bint *n, double *a, bint *lda, double *vl, double *vu,
             bint *il, bint *iu, double *abstol, bint *m, double *w, double *z, bint *ldz, bint *isuppz,
             double *work, bint *lwork, bint *iwork, bint *liwork, bint *info);

void dsytrf_(char *uplo, bint *n, double *a, bint *lda, bint *ipiv, double *work, bint *lwork, bint *info);

void dsytri_(char *uplo, bint *n, double *a, bint *lda, const bint *ipiv, double *work, bint *info);

void dsycon_(char *uplo, bint *n, const double *a, bint *lda, const bint *ipiv, double *anorm, double *rcond,
             double *work, bint *iwork, bint *info);

double dlansy_(char *norm, char *uplo, bint *n, const double *a, bint *lda, double *work);


#ifdef __cplusplus
}
#endif // __cplusplus


double C_dnrm2(bint n, const double *x, bint incx);

double C_ddot(bint n, const double *x, bint incx, const double *y, bint incy);

int C_dscal(bint n, double alpha, double *x, bint incx);

int C_daxpy(bint n, double alpha, const double *x, bint incx, double *y, bint incy);

int C_dgemv(char trans, bint m, bint n, double alpha, const double *a, bint lda,
            const double *x, bint incx, double beta, double *y, bint incy);

int C_dgemm(char transa, char transb, bint m, bint n, bint k, double alpha, const double *a, bint lda,
            const double *b, bint ldb, double beta, double *c, bint ldc);

int C_dsyrk(char uplo, char trans, bint n, bint k, double alpha, const double *a, bint lda,
            double beta, double *c, bint ldc);

int C_dtrsv(char uplo, char trans, char diag, bint n, const double *a, bint lda, double *x, bint incx);

int C_dtrsm(char side, char uplo, char transa, char diag, bint m, bint n, double alpha,
            const double *a, bint lda, double *b, bint ldb);

int C_dgels(char trans, bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb);

int C_dgelsy(bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb, bint *jpvt, double rcond, bint *rank);

int C_dgeqp3(bint m, bint n, double *a, bint lda, bint *jpvt, double *tau);

int C_dorgqr(bint m, bint n, bint k, double *a, bint lda, const double *tau);

int C_dpotrf(char uplo, bint n, double *a, bint lda);

int C_dsyevr(char jobz, char range, char uplo, bint n, double *a, bint lda, double vl, double vu, bint il, bint iu,
             double abstol, bint *m, double *w, double *z, bint ldz, bint *isuppz);

int C_dsytrf(char uplo, bint n, double *a, bint lda, bint *ipiv);

int C_dsytri(char uplo, bint n, double *a, bint lda, const bint *ipiv);

int C_dsycon(char uplo, bint n, const double *a, bint lda, const bint *ipiv, double anorm, double *rcond);

double C_dlansy(char norm, char uplo, bint n, const double *a, bint lda);


// Symmetric Matrix Inverse
//
//   Input
//     n   the order of the symmetric matrix A
//     a   the n*n symmetric matrix A
//
int M_dsyinv(int n, double *a);


#endif // LAPACK_H
