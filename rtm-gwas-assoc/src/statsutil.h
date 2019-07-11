#ifndef STATSUTIL_H
#define STATSUTIL_H


#include <cmath>
#include <vector>
#include <algorithm>


// F Cumulative Distribution Function Upper Tail
double fpval(double x, double df1, double df2);


// Dummy Variables for Grouping Variable g (integers: 0,1,...)
//
// 0/-1/1 coding, full rank, drop last level
int idummy1(const std::vector<int> &g, std::vector< std::vector<double> > &x);

// 0/1 coding, full rank, drop last level
int idummy2(const std::vector<int> &g, std::vector< std::vector<double> > &x);

// 0/1 coding, overdetermined
int idummy3(const std::vector<int> &g, std::vector< std::vector<double> > &x);


// Calculate t-th Percentile (0 <= t <= 100)
//   x must be sorted in ascending order
//   SAS PROC UNIVARIATE PCTLDEF=5
//   j + g = np = n * t / 100
//   g = 0: (x[j] + x[j+1]) / 2
//   g > 0: x[j+1]
double percentile(int t, const std::vector<double> &x);


// generate factor variable
template<typename T>
std::vector<int> factor(const std::vector<T> &v)
{
    auto u = v;
    std::sort(u.begin(), u.end());
    u.erase(std::unique(u.begin(), u.end()), u.end());

    std::vector<int> gi;

    gi.reserve(v.size());
    for (auto &e : v) {
        auto itr = std::find(u.begin(), u.end(), e);
        gi.push_back(itr - u.begin());
    }

    return gi;
}


// generate factor variable with level names
template<typename T>
void factor(const std::vector<T> &v, std::vector<T> &gn, std::vector<int> &gi)
{
    auto u = v;
    std::sort(u.begin(), u.end());
    u.erase(std::unique(u.begin(), u.end()), u.end());

    gi.clear();
    gi.reserve(v.size());

    for (auto &e : v) {
        auto itr = std::find(u.begin(), u.end(), e);
        gi.push_back(itr - u.begin());
    }

    gn.swap(u);
}


// Sum
template<typename T>
double calc_sum(const std::vector<T> &x)
{
    double s = 0;
    for (auto e : x)
        s += e;
    return s;
}


// Mean
template<typename T>
double calc_mean(const std::vector<T> &x)
{
    return calc_sum(x) / x.size();
}


// Corrected Sum of Squares
template<typename T>
double calc_css(const std::vector<T> &x)
{
    auto mx = calc_mean(x);
    double css = 0;
    for (auto e : x)
        css += (e-mx)*(e-mx);
    return css;
}


// Variance
template<typename T>
double calc_var(const std::vector<T> &x)
{
    double n = x.size();
    return calc_css(x) / (n-1);
}


// Standard Deviation
template<typename T>
double calc_std(const std::vector<T> &x)
{
    return std::sqrt( calc_var(x) );
}


#endif // STATSUTIL_H
