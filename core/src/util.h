#ifndef UTIL_H
#define UTIL_H

#include <set>
#include <numeric>
#include <algorithm>
#include "main.h"

template<typename T>
T number(const string &s, bool *ok = nullptr);

template<typename T1, typename T2>
bool contain(const vector<T1> &v, const T2 &a)
{
    return std::find(v.begin(), v.end(), a) != v.end();
}

template<typename T1, typename T2>
size_t index(const vector<T1> &v, const T2 &a)
{
    return std::distance(v.begin(), std::find(v.begin(), v.end(), a));
}

template<typename T>
vector<size_t> order(const vector<T> &v)
{
    vector<size_t> z(v.size());
    std::iota(z.begin(), z.end(), size_t(0));
    std::sort(z.begin(), z.end(), [&v](size_t i, size_t j) { return v[i] < v[j]; });
    return z;
}

template<typename T>
vector<T> unique(vector<T> v)
{
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

template<typename T>
vector<T> stable_unique(vector<T> v)
{
    std::set<T> seen;

    auto last = v.begin();
    for (auto itr = v.begin(); itr != v.end(); ++itr) {
        if (seen.insert(*itr).second) {
            if (last != itr)
                *last = *itr;
            ++last;
        }
    }

    v.erase(last, v.end());

    return v;
}

template<typename T1, typename T2>
vector<T1> subset(const vector<T1> &v, const vector<T2> &idx)
{
    vector<T1> z;
    z.reserve(idx.size());
    for (auto i : idx)
        z.push_back(v[i]);
    return z;
}

template<typename T>
vector<T> subset(const vector<T> &v, const vector<bool> &mask)
{
    vector<T> z;
    z.reserve( std::count(mask.begin(), mask.end(), true) );
    auto n = v.size();
    for (size_t i = 0; i < n; ++i) {
        if ( mask[i] )
            z.push_back(v[i]);
    }
    return z;
}

template<typename T>
vector<T> intersect(vector<T> a, vector<T> b)
{
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    vector<T> c;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(c));
    return c;
}

#endif // UTIL_H
