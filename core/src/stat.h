#ifndef STAT_H
#define STAT_H

#include <cmath>
#include "main.h"
#include "util.h"

double fcdf(double x, double df1, double df2, bool lower_tail);

double zeroin(double a, double b, double fa, double fb, double (*func)(double x, void *info),
              void *info, double &tol, int &maxit);

void design(int method, const vector<int> &g, vector< vector<double> > &x);

double fdr_threshold(double a, const vector<double> &v);

template<typename T>
vector<int> factor(const vector<T> &v)
{
    vector<int> g;
    auto u = unique(v);
    g.reserve(v.size());
    for (auto &e : v)
        g.push_back(index(u,e));
    return g;
}

template<typename T>
void factor(const vector<T> &v, vector<T> &gn, vector<int> &g)
{
    auto u = unique(v);
    g.clear();
    g.reserve(v.size());
    for (auto &e : v)
        g.push_back(index(u,e));
    gn.swap(u);
}

template<typename T>
double sum(const vector<T> &v)
{
    double s = 0;
    for (double e : v)
        s += e;
    return s;
}

template<typename T>
double mean(const vector<T> &v)
{
    return sum(v) / v.size();
}

template<typename T>
double sumsq(const vector<T> &v)
{
    double ss = 0;
    for (double e : v)
        ss += e*e;
    return ss;
}

template<typename T>
double sumsqc(const vector<T> &v)
{
    double m = mean(v);
    double ss = 0;
    for (double e : v)
        ss += (e-m)*(e-m);
    return ss;
}

template<typename T>
double variance(const vector<T> &v)
{
    return sumsqc(v) / (v.size() - 1);
}

template<typename T>
double stddev(const vector<T> &v)
{
    return std::sqrt(variance(v));
}

#endif // STAT_H
