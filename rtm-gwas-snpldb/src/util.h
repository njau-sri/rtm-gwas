#ifndef UTIL_H
#define UTIL_H


#include <set>
#include <string>
#include <vector>
#include <numeric>
#include <iterator>
#include <algorithm>


bool starts_with(const std::string &s1, const std::string &s2);

bool ends_with(const std::string &s1, const std::string &s2);

std::vector<std::string> split(const std::string &str, const std::string &sep);

std::string join(const std::vector<std::string> &vs, const std::string &sep);


template<typename T1, typename T2>
std::size_t index(const std::vector<T1> &v, const T2 &a)
{
    std::size_t i = 0, n = v.size();

    for (i = 0; i < n; ++i)
        if (v[i] == a)
            break;

    return i;
}

template<typename T1, typename T2>
std::size_t count(const std::vector<T1> &v, const T2 &a)
{
    size_t c = 0;

    auto n = v.size();
    for (size_t i = 0; i < n; ++i) {
        if (v[i] == a)
            ++c;
    }

    return c;
}

template<typename T1, typename T2>
std::size_t count_if(const std::vector<T1> &v, const T2 &pred)
{
    size_t c = 0;

    auto n = v.size();
    for (size_t i = 0; i < n; ++i) {
        if ( pred(v[i]) )
            ++c;
    }

    return c;
}

template<typename T>
std::vector<std::size_t> order(const std::vector<T> &v, bool decreasing = false)
{
    std::vector<std::size_t> z(v.size());

    std::iota(z.begin(), z.end(), std::size_t(0));

    if ( decreasing )
        std::sort(z.begin(), z.end(), [&v](std::size_t i, std::size_t j) { return v[i] > v[j]; });
    else
        std::sort(z.begin(), z.end(), [&v](std::size_t i, std::size_t j) { return v[i] < v[j]; });

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

    v.erase(std::remove_if(v.begin(), v.end(), [&z] (const T &a) { return ! z.insert(a).second; }), v.end());

    return v;
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
std::vector<T1> subset(const std::vector<T1> &vec, const std::vector<T2> &idx)
{
    std::vector<T1> out;

    out.reserve(idx.size());

    for (auto i : idx)
        out.push_back(vec[i]);

    return out;
}

template<typename T>
std::vector<T> subset(const std::vector<T> &vec, const std::vector<char> &mask)
{
    std::vector<T> out;

    out.reserve( mask.size() - count(mask, 0) );

    auto n = vec.size();
    for (std::size_t i = 0; i < n; ++i)
        if ( mask[i] )
            out.push_back(vec[i]);

    return out;
}


#endif // UTIL_H
