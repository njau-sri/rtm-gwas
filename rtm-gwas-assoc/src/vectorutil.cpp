#include "vectorutil.h"

std::vector<isize_t> seq(isize_t a, isize_t b)
{
    std::vector<isize_t> v;

    if (a < b) {
        isize_t n = b - a + 1;
        v.reserve(n);
        for (isize_t i = a; i <= b; ++i)
            v.push_back(i);
    }
    else {
        isize_t n = a - b + 1;
        v.reserve(n);
        for (isize_t i = a; i >= b; --i)
            v.push_back(i);
    }

    return v;
}

std::vector<isize_t> seq(isize_t a, isize_t b, isize_t by)
{
    std::vector<isize_t> v;

    if (a < b) {
        isize_t n = (b - a + by) / by;
        v.reserve(n);
        for (isize_t i = a; i <= b; i += by)
            v.push_back(i);
    }
    else {
        isize_t n = (a - b + by) / by;
        v.reserve(n);
        for (isize_t i = a; i >= b; i -= by)
            v.push_back(i);
    }

    return v;
}

void design1(const std::vector<isize_t> &g, std::vector< std::vector<double> > &x)
{
    auto p = std::minmax_element(g.begin(), g.end());

    if (*p.first != 0 || *p.second <= 0) {
        x.clear();
        return;
    }

    isize_t n = length(g);
    isize_t m = *p.second;

    x.assign(m, std::vector<double>(n, 0.0));

    for (isize_t i = 0; i < n; ++i) {
        isize_t k = g[i];
        if (k != m)
            x[k][i] = 1.0;
        else
            for (isize_t j = 0; j < m; ++j)
                x[j][i] = -1.0;
    }
}

void design2(const std::vector<isize_t> &g, std::vector< std::vector<double> > &x)
{
    auto p = std::minmax_element(g.begin(), g.end());

    if (*p.first != 0 || *p.second <= 0) {
        x.clear();
        return;
    }

    isize_t n = length(g);
    isize_t m = *p.second;

    x.assign(m, std::vector<double>(n, 0.0));

    for (isize_t i = 0; i < n; ++i) {
        isize_t k = g[i];
        if (k != m)
            x[k][i] = 1.0;
    }
}

void design3(const std::vector<isize_t> &g, std::vector< std::vector<double> > &x)
{
    auto p = std::minmax_element(g.begin(), g.end());

    if (*p.first != 0 || *p.second <= 0) {
        x.clear();
        return;
    }

    isize_t n = length(g);
    isize_t m = *p.second;

    x.assign(m+1, std::vector<double>(n, 0.0));

    for (isize_t i = 0; i < n; ++i)
        x[g[i]][i] = 1.0;
}
