#ifndef VECTORUTIL_H
#define VECTORUTIL_H

#include <set>
#include <numeric>
#include <utility>
#include <algorithm>

#include "types.h"

template<typename T>
bool contains(const std::vector<T> &v, const T &a)
{
    return std::find(v.begin(), v.end(), a) != v.end();
}

template<typename T>
isize_t index(const std::vector<T> &v, const T &a)
{
    auto itr = std::find(v.begin(), v.end(), a);
    return itr != v.end() ? itr - v.begin() : -1;
}

template<typename T>
isize_t index_min(const std::vector<T> &v)
{
    auto itr = std::min_element(v.begin(), v.end());
    return itr != v.end() ? itr - v.begin() : -1;
}

template<typename T>
isize_t index_max(const std::vector<T> &v)
{
    auto itr = std::max_element(v.begin(), v.end());
    return itr != v.end() ? itr - v.begin() : -1;
}

template<typename T>
std::pair<isize_t, isize_t> index_minmax(const std::vector<T> &v)
{
    auto p = std::minmax_element(v.begin(), v.end());
    isize_t i = p.first != v.end() ? p.first - v.begin() : -1;
    isize_t j = p.second != v.end() ? p.second - v.begin() : -1;
    return { i, j };
}

template<typename T>
isize_t count(const std::vector<T> &v, const T &a)
{
    return std::count(v.begin(), v.end(), a);
}

template<typename T>
std::vector<isize_t> order_asc(const std::vector<T> &v)
{
    std::vector<isize_t> z(v.size());

    std::iota(z.begin(), z.end(), (isize_t) 0);

    std::sort(z.begin(), z.end(), [&v](isize_t i, isize_t j) { return v[i] < v[j]; });

    return z;
}

template<typename T>
std::vector<isize_t> order_desc(const std::vector<T> &v)
{
    std::vector<isize_t> z(v.size());

    std::iota(z.begin(), z.end(), (isize_t) 0);

    std::sort(z.begin(), z.end(), [&v](isize_t i, isize_t j) { return v[i] > v[j]; });

    return z;
}

template<typename T>
std::vector<T> unique(std::vector<T> v)
{
    std::sort(v.begin(), v.end());

    v.erase(std::unique(v.begin(), v.end()), v.end());

    return v;
}

template<typename T>
std::vector<T> stable_unique(std::vector<T> v)
{
    std::set<T> z;

    v.erase(std::remove_if(v.begin(), v.end(),
                           [&z](const T &a) { return !z.insert(a).second; }),
            v.end());

    return v;
}

template<typename T>
bool has_duplicate(const std::vector<T> &v)
{
    std::set<T> z;
    for (auto &e : v) {
        if (!z.insert(e).second)
            return true;
    }
    return false;
}

template<typename T>
std::vector<T> intersect(std::vector<T> a, std::vector<T> b)
{
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());

    std::vector<T> c;

    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(c));

    return c;
}

template<typename T1, typename T2>
std::vector<T1> subset(const std::vector<T1> &v, const std::vector<T2> &idx)
{
    std::vector<T1> z;

    z.reserve(idx.size());

    for (auto i : idx)
        z.push_back(v[i]);

    return z;
}

template<typename T>
std::vector<T> subset(const std::vector<T> &v, const std::vector<char> &mask)
{
    std::vector<T> z;

    z.reserve(mask.size() - count(mask, (char) 0));

    isize_t n = length(v);
    for (isize_t i = 0; i < n; ++i)
        if (mask[i])
            z.push_back(v[i]);

    return z;
}

// create index vector from grouping variable
template<typename T>
std::vector<isize_t> grpidx(const std::vector<T> &v)
{
    std::vector<isize_t> g;
    g.reserve(v.size());

    auto u = unique(v);

    for (auto &e : v)
        g.push_back(index(u, e));

    return g;
}

// create index vector from grouping variable with group names
template<typename T>
void grpidx(const std::vector<T> &v, std::vector<isize_t> &g, std::vector<T> &gn)
{
    g.clear();
    g.reserve(v.size());

    auto u = unique(v);

    for (auto &e : v)
        g.push_back(index(u, e));

    gn.swap(u);
}

// T f(const T &a)
template<typename T, typename F>
std::vector<T> apply(const std::vector<T> &v, F f)
{
    std::vector<T> w(v.size());
    std::transform(v.begin(), v.end(), w.begin(), f);
    return w;
}

template<typename T>
std::vector<T> rep(const std::vector<T> &v, isize_t n)
{
    std::vector<T> w;
    w.reserve(length(v)*n);
    for (isize_t i = 0; i < n; ++i)
        w.insert(w.end(), v.begin(), v.end());
    return w;
}

std::vector<isize_t> seq(isize_t a, isize_t b);

std::vector<isize_t> seq(isize_t a, isize_t b, isize_t by);

// design matrix (grpidx, 0,1,2,3,...)
// 1: full rank, 0/1/-1, drop last level
// 2: full rank, 0/1, drop last level
// 3: overdetermined, 0/1
void design1(const std::vector<isize_t> &g, std::vector< std::vector<double> > &x);
void design2(const std::vector<isize_t> &g, std::vector< std::vector<double> > &x);
void design3(const std::vector<isize_t> &g, std::vector< std::vector<double> > &x);

#endif // VECTORUTIL_H
