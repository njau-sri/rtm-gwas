#ifndef STATUTIL_H
#define STATUTIL_H

#include <cmath>
#include "types.h"

// 1 - fcdf
double fpval(double x, double df1, double df2);

// floating-point relative accuracy
double eps(double x);

// a^2
template<typename T>
T square(T a)
{
    return a * a;
}

// a^3
template<typename T>
T cube(T a)
{
    return a * a * a;
}

// sum
template<typename T>
T sum(isize_t n, const T *x)
{
    T s = (T) 0;
    for (isize_t i = 0; i < n; ++i)
        s += x[i];
    return s;
}

// mean
template<typename T>
double mean(isize_t n, const T *x)
{
    return sum(n, x) / (double) n;
}

// sum of squares uncorrected
template<typename T>
T uss(isize_t n, const T *x)
{
    T s = (T) 0;
    for (isize_t i = 0; i < n; ++i)
        s += square(x[i]);
    return s;
}

// sum of squares corrected for mean
//    $$ \sum (x_i - \bar x)^2 $$
template<typename T>
double css(isize_t n, const T *x)
{
    double s = 0.0;
    double m = mean(n, x);
    for (isize_t i = 0; i < n; ++i)
        s += square(x[i] - m);
    return s;
}

// variance
template<typename T>
double var(isize_t n, const T *x)
{
    return css(n, x) / (n - 1);
}

// standard deviation
template<typename T>
double sd(isize_t n, const T *x)
{
    return std::sqrt(var(n, x));
}

// p-th percentile (0 <= p <= 100)
//    x is assumed in ascending order
//    proc univariate pctldef=5
template<typename T>
T pctl(int p, isize_t n, const T *x)
{
    if (p <= 0)
        return x[0];

    if (p >= 100)
        return x[n-1];

    double t = (double) p / 100;
    double np = n * t;
    double i = 0;
    double g = std::modf(np, &i);
    isize_t j = (isize_t) i;
    return g == 0.0 ? (x[j-1] + x[j]) / 2 : x[j];
}

#endif // STATUTIL_H
