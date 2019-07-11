#include <iostream>
#include "stepreg.h"
#include "lsfit.h"
#include "statsutil.h"


using std::ptrdiff_t;
using std::size_t;


namespace {

void forward_impl(int mtc, double sle, double maxrsq, const std::vector<double> &y, const std::vector<double> &x0,
                  const std::vector<const double*> &xs, const std::vector<std::size_t> &cs,
                  std::vector<double> &ps, std::vector<std::size_t> &in)
{
    auto n = y.size();
    auto m = xs.size();

    double dfe0 = 0, sse0 = 0;
    std::vector<double> b0;
    lsfit(y, x0, b0, dfe0, sse0);

    double sst = calc_css(y);
    double rsq = 1 - sse0 / sst;
    std::cerr << "INFO: forward selection step 0: " << rsq << "\n";

    auto x = x0;
    auto nq0 = x.size();
    std::vector<double> b1;
    std::vector<char> ignore(m,0);

    for (size_t step = 0; step < m; ++step) {
        size_t idx = 0;
        double pval = 1;
        double dfe = 0, sse = 0;

        for (size_t j = 0; j < m; ++j) {
            if ( ignore[j] )
                continue;

            x.resize(nq0);
            auto nq1 = n * cs[j];
            x.insert(x.end(), xs[j], xs[j] + nq1);

            double dfe1 = 0, sse1 = 0;
            lsfit(y, x, b1, dfe1, sse1);

            if (dfe1 > 0.0 && sse1 > 0.0 && dfe0 > dfe1 && sse0 > sse1) {
                auto f = ((sse0 - sse1) / (dfe0 - dfe1)) / (sse1 / dfe1);
                auto p = fpval(f, dfe0 - dfe1, dfe1);
                if (p < pval) {
                    idx = j;
                    pval = p;
                    dfe = dfe1;
                    sse = sse1;
                }
                ps[j] = p;
            }
        }

        double alpha = sle;
        if (mtc > 0) {
            if (mtc == 1)
                alpha /= m;
            else if (mtc == 2)
                alpha *= (step + 1) / m;
            else if (mtc == 3)
                alpha /= (m - step);
            else
                alpha /= m;
        }

        if (pval > alpha)
            break;

        rsq = 1 - sse / sst;
        if (rsq > maxrsq)
            break;

        in.push_back(idx);
        ignore[idx] = 1;

        x.resize(nq0);
        auto nq1 = n * cs[idx];
        x.insert(x.end(), xs[idx], xs[idx] + nq1);
        nq0 = x.size();
        dfe0 = dfe;
        sse0 = sse;

        std::cerr << "INFO: forward selection step " << step+1 << ": " << rsq << " " << pval << "\n";
    }
}

void forward_omp_impl(int mtc, double sle, double maxrsq, const std::vector<double> &y, const std::vector<double> &x0,
                      const std::vector<const double*> &xs, const std::vector<std::size_t> &cs,
                      std::vector<double> &ps, std::vector<std::size_t> &in)
{
    auto n = y.size();
    auto m = xs.size();

    double dfe0 = 0, sse0 = 0;
    std::vector<double> b0;
    lsfit(y, x0, b0, dfe0, sse0);

    double sst = calc_css(y);
    double rsq = 1 - sse0 / sst;
    std::cerr << "INFO: forward selection step 0: " << rsq << "\n";

    auto x = x0;
    std::vector<char> ignore(m,0);
    std::vector<double> dfe(m);
    std::vector<double> sse(m);
    std::vector<double> pval(m);

    for (size_t step = 0; step < m; ++step) {
        pval.assign(m, 1);

        // in earlier OpenMP specifications (<3.0), unsigned integer is not allowed in loop construct
        auto m2 = static_cast<ptrdiff_t>(m);

        #pragma omp parallel for
        for (ptrdiff_t j2 = 0; j2 < m2; ++j2) {
            auto j = static_cast<size_t>(j2);
            if ( ignore[j] )
                continue;

            auto xx = x;
            auto nq = n * cs[j];
            xx.insert(xx.end(), xs[j], xs[j] + nq);

            std::vector<double> b1;
            double dfe1 = 0, sse1 = 0;
            lsfit(y, xx, b1, dfe1, sse1);

            if (dfe1 > 0.0 && sse1 > 0.0 && dfe0 > dfe1 && sse0 > sse1) {
                auto f = ((sse0 - sse1) / (dfe0 - dfe1)) / (sse1 / dfe1);
                ps[j] = pval[j] = fpval(f, dfe0 - dfe1, dfe1);
            }

            dfe[j] = dfe1;
            sse[j] = sse1;
        }

        double alpha = sle;
        if (mtc > 0) {
            if (mtc == 1)
                alpha /= m;
            else if (mtc == 2)
                alpha *= (step + 1) / m;
            else if (mtc == 3)
                alpha /= (m - step);
            else
                alpha /= m;
        }

        auto itr = std::min_element(pval.begin(), pval.end());
        if (*itr > alpha)
            break;

        auto j = static_cast<size_t>( itr - pval.begin() );

        rsq = 1 - sse[j] / sst;
        if (rsq > maxrsq)
            break;

        in.push_back(j);
        ignore[j] = 1;

        auto nq = n * cs[j];
        x.insert(x.end(), xs[j], xs[j] + nq);
        dfe0 = dfe[j];
        sse0 = sse[j];

        std::cerr << "INFO: forward selection step " << step+1 << ": " << rsq << " " << *itr << "\n";
    }
}


void backward_impl(double sls, const std::vector<double> &y, const std::vector<double> &x0,
                   const std::vector<const double*> &xs, const std::vector<std::size_t> &cs,
                   std::vector<double> &ps, std::vector<std::size_t> &in)
{
    auto n = y.size();

    auto x = x0;
    auto nq0 = x.size();

    int step = 0;
    std::vector<double> b0, b1;

    while ( ! in.empty() ) {
        x.resize(nq0);
        for (size_t i : in) {
            auto nq1 = n * cs[i];
            x.insert(x.end(), xs[i], xs[i] + nq1);
        }

        double dfe1 = 0, sse1 = 0;
        lsfit(y, x, b1, dfe1, sse1);

        if (dfe1 == 0.0 || sse1 == 0.0) {
            std::cerr << "ERROR: invalid initial model for backward elimination\n";
            return;
        }

        std::vector<double> pval;

        for (size_t i : in) {
            x.resize(nq0);
            for (size_t j : in) {
                if (j != i) {
                    auto nq1 = n * cs[j];
                    x.insert(x.end(), xs[j], xs[j] + nq1);
                }
            }

            double dfe0 = 0, sse0 = 0;
            lsfit(y, x, b0, dfe0, sse0);

            if (dfe1 > 0.0 && sse1 > 0.0 && dfe0 > dfe1 && sse0 > sse1) {
                auto f = ((sse0 - sse1) / (dfe0 - dfe1)) / (sse1 / dfe1);
                auto p = fpval(f, dfe0 - dfe1, dfe1);
                pval.push_back(p);
            }
        }

        auto itr = std::max_element(pval.begin(), pval.end());
        if (*itr <= sls)
            break;

        auto j = itr - pval.begin();
        ps[ in[j] ] = *itr;

        in.erase(in.begin() + j);
        std::cerr << "INFO: backward elimination step " << ++step << ": " << *itr << "\n";
    }
}

} // namespace


void StepReg::forward(const std::vector<double> &x0, const std::vector<double> &x1, const std::vector<std::size_t> &c1,
                      const std::vector<double> &y, std::vector<double> &ps, std::vector<std::size_t> &in)
{
    auto n = y.size();
    auto m = c1.size();

    in.clear();
    ps.assign(m, 1);

    std::vector<const double*> xs(m);
    xs[0] = x1.data();
    for (size_t i = 1; i < m; ++i) {
        auto nq = n * c1[i-1];
        xs[i] = xs[i-1] + nq;
    }

    forward_impl(mtc, sle, rsq, y, x0, xs, c1, ps, in);
}


void StepReg::forward_omp(const std::vector<double> &x0, const std::vector<double> &x1, const std::vector<std::size_t> &c1,
                          const std::vector<double> &y, std::vector<double> &ps, std::vector<std::size_t> &in)
{
    auto n = y.size();
    auto m = c1.size();

    in.clear();
    ps.assign(m, 1);

    std::vector<const double*> xs(m);
    xs[0] = x1.data();
    for (size_t i = 1; i < m; ++i) {
        auto nq = n * c1[i-1];
        xs[i] = xs[i-1] + nq;
    }

    forward_omp_impl(mtc, sle, rsq, y, x0, xs, c1, ps, in);
}


void StepReg::backward(const std::vector<double> &x0, const std::vector<double> &x1, const std::vector<std::size_t> &c1,
                       const std::vector<double> &y, std::vector<double> &ps, std::vector<std::size_t> &in)
{
    auto n = y.size();
    auto m = c1.size();

    if (ps.size() != m)
        ps.assign(m, 1);

    if ( in.empty() )
        return;

    std::vector<const double*> xs(m);
    xs[0] = x1.data();
    for (size_t i = 1; i < m; ++i) {
        auto nq = n * c1[i-1];
        xs[i] = xs[i-1] + nq;
    }

    backward_impl(sls, y, x0, xs, c1, ps, in);
}
