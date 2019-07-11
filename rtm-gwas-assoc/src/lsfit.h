#ifndef LSFIT_H
#define LSFIT_H

#include <vector>


// Least Squares Fitting
//
// y = Xb + e
//
//   Input
//     y   the n*1 vector y
//     x   the n*q matrix X in column-major layout
//   Output
//     b     the q*1 vector b
//     dfe   the residual degrees of freedom
//     sse   the residual sum of squares
//
int lsfit(const std::vector<double> &y, const std::vector<double> &x, std::vector<double> &b, double &dfe, double &sse);


// Least Squares Fitting with Equality Constraints
//
// y = Xb + e  where  Zb = 0
//
//   Input
//     x   the n*q matrix X in column-major layout
//     y   the n*1 vector y
//     z   the p*q matrix Z in row-major layout
//   Output
//     b     q*1 vector b
//     dfe   residual degrees of freedom
//     sse   residual sum of squares
//
int lsfitc(const std::vector<double> &y, const std::vector<double> &x, std::vector<double> z,
           std::vector<double> &b, double &dfe, double &sse);


// Generalized Least Squares Fitting (GLS)
//
// b = argmin (y-Xb)' V^-1 (y-Xb)
//
// b = (X' * V^-1 * X)^-1 * X' * V^-1 * y

int glsfit(const std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &v,
           std::vector<double> &b);


// F-statistic for GLS
//
//   H0: Mb = 0
//
//        (Mb)' [M (X' V^-1 X)^-1 M']^-1 (Mb)
//   F = -------------------------------------
//                         p
//
//   b = (X' * V^-1 * X)^-1 * X' * V^-1 * y
//
//   Input
//     p   the number of columns of matrix X1
//     y   the n*1 vector y
//     x   the n*q matrix X = [X0 X1]
//     v   the n*n symmetric matrix V
//   Output
//     f   the F-statistic with df1 = p, df2 = n - q; f < 0 on error
//
double glsfstat(std::size_t p, const std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &v);


#endif // LSFIT_H
