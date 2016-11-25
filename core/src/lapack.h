#ifndef LAPACK_H
#define LAPACK_H

#ifdef BLAS_INT64
    #include <cstddef>
    using bint = std::ptrdiff_t;
#else
    using bint = int;
#endif

double la_dnrm2(bint n, const double *x, bint incx);

double la_ddot(bint n, const double *x, bint incx, const double *y, bint incy);

int la_dscal(bint n, double alpha, double *x, bint incx);

int la_daxpy(bint n, double alpha, const double *x, bint incx, double *y, bint incy);

int la_dgemv(char trans, bint m, bint n, double alpha, const double *a, bint lda, const double *x, bint incx,
             double beta, double *y, bint incy);

int la_dgemm(char transa, char transb, bint m, bint n, bint k, double alpha, const double *a, bint lda,
             const double *b, bint ldb, double beta, double *c, bint ldc);

int la_dsyrk(char uplo, char trans, bint n, bint k, double alpha, const double *a, bint lda,
             double beta, double *c, bint ldc);

int la_dgels(char trans, bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb);

int la_dgelsy(bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb,
              bint *jpvt, double rcond, bint *rank);

int la_dsyevr(char jobz, char range, char uplo, bint n, double *a, bint lda, double vl, double vu, bint il, bint iu,
              double abstol, bint *m, double *w, double *z, bint ldz, bint *isuppz);

int la_dtrtrs(char uplo, char trans, char diag, bint n, bint nrhs, const double *a, bint lda, double *b, bint ldb);

int la_dpotrf(char uplo, bint n, double *a, bint lda);

int la_dgeqrf(bint m, bint n, double *a, bint lda, double *tau);

int la_dgeqp3(bint m, bint n, double *a, bint lda, bint *jpvt, double *tau);

int la_dormqr(char side, char trans, bint m, bint n, bint k, const double *a, bint lda,
              const double *tau, double *c, bint ldc);

int la_dorgqr(bint m, bint n, bint k, double *a, bint lda, const double *tau);

int la_dsytrf(char uplo, bint n, double *a, bint lda, bint *ipiv);

int la_dsytri(char uplo, bint n, double *a, bint lda, const bint *ipiv);

#endif // LAPACK_H
