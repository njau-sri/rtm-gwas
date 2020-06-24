#include "statutil.h"

#include <algorithm>

#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
#include <boost/math/distributions/fisher_f.hpp>

#include "lapack.h"
#include "vectorutil.h"

double fpval(double x, double df1, double df2)
{
    boost::math::fisher_f f(df1, df2);
    return boost::math::cdf(boost::math::complement(f, x));
}

double eps(double x)
{
    if (!std::isfinite(x))
        return std::numeric_limits<double>::quiet_NaN();

    if (std::fabs(x) < std::numeric_limits<double>::min())
        return std::numeric_limits<double>::denorm_min();

    int exponent;
    std::frexp(x, &exponent);

    int digits = std::numeric_limits<double>::digits;

    return std::ldexp(1, exponent - digits);
}
