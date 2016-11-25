#include "stat.h"

extern "C"
{

// R src/library/stats/src/zeroin.c, http://www.R-project.org/
double R_zeroin2(double ax, double bx, double fa, double fb,
                 double (*f)(double x, void *info), void *info,
                 double *Tol, int *Maxit);

// R standalone Rmath library, http://www.R-project.org/
double pf(double x, double df1, double df2, int lower_tail, int log_p);

} // extern "C"

double fcdf(double x, double df1, double df2, bool lower_tail)
{
    return pf(x, df1, df2, lower_tail ? 1 : 0, 0);
}

double zeroin(double a, double b, double fa, double fb, double (*func)(double x, void *info),
              void *info, double &tol, int &maxit)
{
    return R_zeroin2(a, b, fa, fb, func, info, &tol, &maxit);
}

// method = 1: 0/-1/1 coding, full rank, drop first level
// method = 2: 0/1 coding, full rank, drop first level
// method = 3: 0/1 coding, overdetermined

void design(int method, const vector<int> &g, vector< vector<double> > &x)
{
    int n = g.size();
    int m = g.empty() ? 0 : *std::max_element(g.begin(),g.end());

    if (method < 1 || method > 3)
        method = 3;

    if ( method == 3 && ! g.empty() )
        m += 1;

    x.assign(m, vector<double>(n, 0.0));

    if (method == 3) {
        for (int i = 0; i < n; ++i)
            x[g[i]][i] = 1.0;
    }
    else {
        for (int i = 0; i < n; ++i) {
            int k = g[i];
            if (k != 0) {
                x[k-1][i] = 1.0;
            }
            else if (method != 2) {
                for (int j = 0; j < m; ++j)
                    x[j][i] = -1.0;
            }
        }
    }
}

// Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate - a practical and powerful approach
//   to multiple testing. J R Stat Soc B 57, 289-300 (1995).
double fdr_threshold(double a, const vector<double> &v)
{
    auto n = v.size();
    for (size_t i = n; i > 0; --i) {
        double p = a*i/n;
        if (v[i-1] <= p)
            return p;
    }
    return a/n;
}
