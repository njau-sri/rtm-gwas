#ifndef LSTSQR_H
#define LSTSQR_H

#include "main.h"

class LstSqr
{
public:
    struct Stats
    {
        vector<double> b;
        double dfe = 0.0;
        double sse = 0.0;
    };

    const Stats& stats() const
    {
        return m_st;
    }

    // X*b = y
    int solve(size_t p, const vector<double> &y, const vector<double> &x);

    // X*b = y,  subject to C*b = 0,  (ct = C')
    int solve(size_t p, size_t q, const vector<double> &y, const vector<double> &x, vector<double> ct);

private:
    Stats m_st;
    vector<double> m_y;
    vector<double> m_x;
};

#endif // LSTSQR_H
