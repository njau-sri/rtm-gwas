#ifndef MATRIX_H
#define MATRIX_H

#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif // EIGEN_NO_DEBUG
#include <eigen3/Eigen/Core>

#include "types.h"

using mat = Eigen::MatrixXd;
using imat = Eigen::MatrixXi;


// Eigen                        // MATLAB
//
// mat a;
// mat b(m, n);
// mat x(m, 1);
//
// a.size();                    // length(a)
// a.rows();                    // size(a, 1)
// a.cols();                    // size(a, 2)
// a(i)                         // a(i+1)
// a(i,j)                       // a(i+1, j+1)
// a.resize(m,n)
// a << 1,2,3,4;
// a << b,b;
//
// mat::Identity(m, n)          // eye(m, n)
// mat::Zero(m, n)              // zeros(m, n)
// mat::Ones(m, n)              // ones(m, n)
// mat::Random(m, n)            // rand(m,n)
// mat::Constant(m, n, val)
//
// a.setIdentity(m, n)          // a = eye(m, n)
// a.setZero(m, n)              // a = zeros(m, n)
// a.setOnes(m, n)              // a = ones(m, n)
// a.setRandom(m, n)            // a = rand(m, n)
// a.setConstant(m, n, val)
//
// x.head(n)                    // x(1:n)
// x.tail(n)                    // x(end-n+1:end)
// x.segment(i,n)               // x(i+1:i+n)
//
// a.row(i)                     // a(i+1, :)
// a.col(j)                     // a(:, j+1)
// a.block(i, j, m, n)          // a(i+1:i+m, j+1:j+n)
//
// a.leftCols(n)                // a(:, 1:n)
// a.middleCols(j, n)           // a(:, j+1:j+n)
// a.rightCols(n)               // a(:, end-n+1:end)
//
// a.topRows(m)                 // a(1:m, :)
// a.middleRows(i, m)           // a(i+1:i+m, :)
// a.bottomRows(m)              // a(end-m+1:end, :)
//
// a.topLeftCorner(m, n)        // a(1:m, 1:n)
// a.topRightCorner(m, n)       // a(1:m, end-n+1:end)
// a.bottomLeftCorner(m, n)     // a(end-m+1:end, 1:n)
// a.bottomRightCorner(m, n)    // a(end-m+1:end, end-n+1:end)
//
// a.row(i) = b.col(j)          // a(i+1, :) = b(:, j+1)
// a.col(j1).swap(b.col(j2))    // a(:, [j1 j2]) = b(:, [j2 j1])
//
// a.transpose()                // a'
// a.transposeInPlace()         // a = a'
// a.diagonal()                 // diag(a)
// x.asDiagonal()               // diag(x)
// a.replicate(m, n)            // repmat(a, m, n)
//
// a.cwiseProduct(b)            // a .* b
// a.array() * b.array()        // a .* b
//
// a.sum()                      // sum(a(:))
// a.colwise().sum()            // sum(a)
//
// a.cast<double>()
//
// mat::Map(ptr, m, n)
//
// horitzontal concatenation
// c.resize(a.rows(), a.cols() + b.cols())  c = [a b]
// c << a, b                                c = [a b]
//
// vertical concatenation
// c.resize(a.rows() + b.rows(), a.cols())  c = [a; b]
// c << a, b                                c = [a; b]


template<typename Derived>
isize_t nrow(const Eigen::EigenBase<Derived> &a)
{
    return (isize_t) a.rows();
}

template<typename Derived>
isize_t ncol(const Eigen::EigenBase<Derived> &a)
{
    return (isize_t) a.cols();
}

template<typename Derived>
isize_t length(const Eigen::EigenBase<Derived> &a)
{
    return (isize_t) a.size();
}

// design matrix (grpidx, 0,1,2,3,...)
// 1: full rank, 0/1/-1, drop last level
// 2: full rank, 0/1, drop last level
// 3: overdetermined, 0/1
void design1(isize_t n, const isize_t *g, mat &x);
void design2(isize_t n, const isize_t *g, mat &x);
void design3(isize_t n, const isize_t *g, mat &x);

// vector operation
double nrm2(const Eigen::Ref<const mat> &x);
double dot(const Eigen::Ref<const mat> &x, const Eigen::Ref<const mat> &y);

int gemv(bool trans, double alpha, const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &x, double beta, mat &y);
int gemm(bool transa, bool transb, double alpha, const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b, double beta, mat &c);
int syrk(bool trans, double alpha, const Eigen::Ref<const mat> &a, double beta, mat &c);

// matrix-matrix product
//   mult, A*B
//   crossprod, A'*B, A'*A
//   tcrossprod, A*B', A*A'
mat mult(const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b);
mat crossprod(const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b);
mat tcrossprod(const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b);
mat crossprod(const Eigen::Ref<const mat> &a);
mat tcrossprod(const Eigen::Ref<const mat> &a);

// matrix norm
//   norm1, one norm, maximum absolute column sum
//   normi, infinity norm, maximum absolute row sum
//   normf, frobenius norm, square root of sum of squares
//   norm2, spectral norm, largest singular value
double norm1(const Eigen::Ref<const mat> &a);
double normi(const Eigen::Ref<const mat> &a);
double normf(const Eigen::Ref<const mat> &a);
double norm2(const Eigen::Ref<const mat> &a);

// condition number
double cond(const Eigen::Ref<const mat> &a);

// matrix inverse
//   inv_sym, symmetric matrix inverse
//   pinv, pseudoinverse
mat inv_sym(const Eigen::Ref<const mat> &a);
mat pinv(double rcond, const Eigen::Ref<const mat> &a);

// Eigen, A*V = V*D
int syev(const Eigen::Ref<const mat> &a, mat &eval);
int syev(const Eigen::Ref<const mat> &a, mat &eval, mat &evec);
int syevr(const Eigen::Ref<const mat> &a, mat &eval);
int syevr(const Eigen::Ref<const mat> &a, mat &eval, mat &evec);

// QR, A = Q*R
int qr(const Eigen::Ref<const mat> &a, mat &q, mat &r);
int qr_eco(const Eigen::Ref<const mat> &a, mat &q, mat &r);

// Pivoted QR, A[:,P] = Q*R
int qr(const Eigen::Ref<const mat> &a, mat &q, mat &r, imat &p);
int qr_eco(const Eigen::Ref<const mat> &a, mat &q, mat &r, imat &p);

// SVD, A = U*S*V'
int gesvd(const Eigen::Ref<const mat> &a, mat &s);
int gesvd(const Eigen::Ref<const mat> &a, mat &u, mat &s, mat &vt);
int gesvd_eco(const Eigen::Ref<const mat> &a, mat &u, mat &s, mat &vt);
int gesdd(const Eigen::Ref<const mat> &a, mat &s);
int gesdd(const Eigen::Ref<const mat> &a, mat &u, mat &s, mat &vt);
int gesdd_eco(const Eigen::Ref<const mat> &a, mat &u, mat &s, mat &vt);

// least norm solution, A * X = B
int gels(const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b, mat &x);
int gelsy(double rcond, const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b, mat &x);
int gelsd(double rcond, const Eigen::Ref<const mat> &a, const Eigen::Ref<const mat> &b, mat &x);

#endif // MATRIX_H
