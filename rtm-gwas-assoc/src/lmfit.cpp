#include "lmfit.h"

#include <memory>
#include <algorithm>

#include "print.h"
#include "lapack.h"
#include "matrix.h"
#include "statutil.h"

int lmfit(isize_t n, isize_t p, const double *y, const double *x, LmFit *out)
{
    auto mx = mat::Map(x, n, p);
    auto my = mat::Map(y, n, 1);

    mat b;
    int rank = gelsy(out->par.rcond, mx, my, b);

    if (rank < 0) {
        out->b.assign(p, 0.0);
        out->dfe = out->sse = 0.0;
        eprint("ERROR: gelsy failed in lmfit (%s:%d): %d\n", __FILE__, __LINE__, rank);
        return 1;
    }

    out->b.assign(b.data(), b.data() + length(b));

    if (rank < n) {
        out->dfe = (double) n - rank;

        b = my;
        gemv(false, -1.0, mx, mat::Map(out->b.data(), p, 1), 1.0, b);

        double r = nrm2(b);
        out->sse = r * r;
    }
    else
        out->dfe = out->sse = 0.0;

    return 0;
}

// Ake Bjorck (2015) Numerical Methods in Matrix Computations, Springer Cham, DOI: 10.1007/978-3-319-05089-8

int lmfit_eq(isize_t n, isize_t p, isize_t q, const double *y, const double *x, const double *z, LmFitEq *out)
{
    auto mx = mat::Map(x, n, p);
    auto mz = mat::Map(z, q, p);

    mat qz, rz;
    imat pz;
    qr(mz.transpose(), qz, rz, pz);

    double tol = out->par.tol;
    isize_t k = (rz.diagonal().cwiseAbs().array() > tol).count();

    LmFit lm;
    lm.par.rcond = out->par.rcond;

    if (k < p) {
        mat xz = mult(mx, qz.rightCols(p-k));

        int info = lmfit(n, p-k, y, xz.data(), &lm);

        if (info != 0) {
            out->b.assign(p, 0.0);
            out->dfe = out->sse = 0.0;
            eprint("ERROR: lmfit failed in lmfit_eq (%s:%d): %d\n", __FILE__, __LINE__, info);
            return 1;
        }

        auto bz = mat::Map(lm.b.data(), length(lm.b), 1);

        mat b;
        gemv(false, 1.0, qz.rightCols(p-k), bz, 0.0, b);

        out->b.assign(b.data(), b.data() + b.size());
    }
    else {
        eprint("WARNING: null space of constraints is empty in lmfit_eq\n");

        int info = lmfit(n, p, y, x, &lm);

        if (info != 0) {
            out->b.assign(p, 0.0);
            out->dfe = out->sse = 0.0;
            eprint("ERROR: lmfit failed in lmfit_eq (%s:%d): %d\n", __FILE__, __LINE__, info);
            return 2;
        }

        out->b = lm.b;
    }

    out->dfe = lm.dfe;
    out->sse = lm.sse;

    return 0;
}

namespace {

    void select(isize_t n, isize_t q, isize_t p, const double *y, const double *w, const double *x, const isize_t *g,
                StepwiseFit *out)
    {
        const auto& par = out->par;

        bool known_error = par.dfe > 0.0 && par.mse > 0.0;

        if (out->ps.empty())
            out->ps.assign(p, 1.0);

        std::vector<char> ignore(p, 0);

        std::vector<const double*> xptr(p);
        xptr[0] = x;
        for (isize_t i = 1; i < p; ++i)
            xptr[i] = xptr[i-1] + n * g[i-1];

        static const isize_t reserve = 1000;

        isize_t p0 = q;
        mat X(n, p0 + reserve);
        X.leftCols(q) = mat::Map(w, n, q);

        LmFit lm;
        lm.par.rcond = par.rcond;

        int info = lmfit(n, p0, y, X.data(), &lm);

        if (info != 0) {
            eprint("ERROR: lmfit failed in select (%s:%d): %d\n", __FILE__, __LINE__, info);
            return;
        }

        double dfe0 = lm.dfe;
        double sse0 = lm.sse;
        double sst = css(n, y);
        eprint("INFO: R2 = %g\n", 1.0 - sse0 / sst);

        if (known_error && dfe0 <= par.dfe)
            return;

        for (isize_t i = 0; i < p; ++i) {
            isize_t next = 0;
            double pval_next = 1.0, dfe_next = 0.0, sse_next = 0.0;
            for (isize_t j = 0; j < p; ++j) {
                if (ignore[j])
                    continue;

                isize_t p1 = p0 + g[j];
                if (p1 > ncol(X)) {
                    isize_t grow = p1 - ncol(X);
                    if (grow < reserve)
                        grow = reserve;
                    X.conservativeResize(Eigen::NoChange, ncol(X) + grow);
                }

                X.middleCols(p0, g[j]) = mat::Map(xptr[j], n, g[j]);

                info = lmfit(n, p1, y, X.data(), &lm);

                if (info != 0) {
                    eprint("ERROR: lmfit failed in select (%s:%d): %d\n", __FILE__, __LINE__, info);
                    return;
                }

                double dfe1 = lm.dfe;
                double sse1 = lm.sse;

                if (known_error && dfe1 <= par.dfe)
                    continue;

                if (dfe1 > 0.0 && sse1 > 0.0 && dfe0 > dfe1 && sse0 > sse1) {
                    double f = (sse0 - sse1) / (dfe0 - dfe1);
                    double df = dfe1;
                    if (known_error) {
                        f /= par.mse;
                        df = par.dfe;
                    }
                    else
                        f /= (sse1 / dfe1);
                    double pval = fpval(f, dfe0 - dfe1, df);
                    if (pval < pval_next) {
                        next = j;
                        pval_next = pval;
                        dfe_next = dfe1;
                        sse_next = sse1;
                    }
                    out->ps[j] = pval;
                }
            }

            double sle = par.sle;
            if (par.mtc == 1)
                sle /= p;
            else if (par.mtc == 2)
                sle = sle * (i + 1) / p;
            else if (par.mtc == 3)
                sle /= (p - i);

            if (pval_next > sle)
                break;

            double rsq = 1 - sse_next / sst;
            if (rsq > par.rsq)
                break;

            out->in.push_back(next);
            ignore[next] = 1;

            X.middleCols(p0, g[next]) = mat::Map(xptr[next], n, g[next]);
            p0 += g[next];
            dfe0 = dfe_next;
            sse0 = sse_next;

            eprint("INFO: selection %td: %td, R2 = %g, P = %g\n", i+1, next+1, rsq, pval_next);
        }
    }

    void select_omp(isize_t n, isize_t q, isize_t p, const double *y, const double *w,
                    const double *x, const isize_t *g, StepwiseFit *out)
    {
        const auto& par = out->par;

        bool known_error = par.dfe > 0.0 && par.mse > 0.0;

        if (out->ps.empty())
            out->ps.assign(p, 1.0);

        std::vector<char> ignore(p, 0);

        std::vector<const double*> xptr(p);
        xptr[0] = x;
        for (isize_t i = 1; i < p; ++i)
            xptr[i] = xptr[i-1] + n * g[i-1];

        static const isize_t reserve = 1000;

        isize_t p0 = q;
        mat x0(n, p0 + reserve);
        x0.leftCols(q) = mat::Map(w, n, q);

        LmFit lm0;
        lm0.par.rcond = par.rcond;

        int info = lmfit(n, p0, y, x0.data(), &lm0);

        if (info != 0) {
            eprint("ERROR: lmfit failed in select_omp (%s:%d): %d\n", __FILE__, __LINE__, info);
            return;
        }

        double dfe0 = lm0.dfe;
        double sse0 = lm0.sse;
        double sst = css(n, y);
        eprint("INFO: R2 = %g\n", 1.0 - sse0 / sst);

        if (known_error && dfe0 <= par.dfe)
            return;

        std::vector<double> dfe(p);
        std::vector<double> sse(p);
        std::vector<double> pval(p);

        for (isize_t i = 0; i < p; ++i) {
            dfe.assign(p, 0.0);
            sse.assign(p, 0.0);
            pval.assign(p, 1.0);

            #pragma omp parallel for
            for (isize_t j = 0; j < p; ++j) {
                if (ignore[j])
                    continue;

                isize_t p1 = p0 + g[j];
                mat x1(n, p1);
                x1.leftCols(p0) = x0;
                x1.rightCols(g[j]) = mat::Map(xptr[j], n, g[j]);

                LmFit lm1;
                lm1.par.rcond = par.rcond;

                info = lmfit(n, p1, y, x1.data(), &lm1);

                if (info != 0) {
                    eprint("ERROR: lmfit failed in select_omp (%s:%d): %d\n", __FILE__, __LINE__, info);
                    continue;
                }

                double dfe1 = lm1.dfe;
                double sse1 = lm1.sse;

                if (known_error && dfe1 <= par.dfe)
                    continue;

                if (dfe1 > 0.0 && sse1 > 0.0 && dfe0 > dfe1 && sse0 > sse1) {
                    double f = (sse0 - sse1) / (dfe0 - dfe1);
                    double df = dfe1;
                    if (known_error) {
                        f /= par.mse;
                        df = par.dfe;
                    }
                    else
                        f /= (sse1 / dfe1);
                    out->ps[j] = pval[j] = fpval(f, dfe0 - dfe1, df);
                    dfe[j] = dfe1;
                    sse[j] = sse1;
                }
            }

            double sle = par.sle;
            if (par.mtc == 1)
                sle /= p;
            else if (par.mtc == 2)
                sle = sle * (i + 1) / p;
            else if (par.mtc == 3)
                sle /= (p - i);

            auto itr = std::min_element(pval.begin(), pval.end());
            if (*itr > sle)
                break;

            isize_t next = itr - pval.begin();

            double rsq = 1 - sse[next] / sst;
            if (rsq > par.rsq)
                break;

            out->in.push_back(next);
            ignore[next] = 1;

            isize_t p1 = p0 + g[next];
            if (p1 > ncol(x0)) {
                isize_t grow = p1 - ncol(x0);
                if (grow < reserve)
                    grow = reserve;
                x0.conservativeResize(Eigen::NoChange, ncol(x0) + grow);
            }

            x0.middleCols(p0, g[next]) = mat::Map(xptr[next], n, g[next]);
            p0 += g[next];
            dfe0 = dfe[next];
            sse0 = sse[next];

            eprint("INFO: selection %td: %td, R2 = %g, P = %g\n", i+1, next+1, rsq, pval[next]);
        }
    }

    void eliminate(isize_t n, isize_t q, isize_t p, const double *y, const double *w,
                   const double *x, const isize_t *g, StepwiseFit *out)
    {
        const auto& par = out->par;

        bool known_error = par.dfe > 0.0 && par.mse > 0.0;

        if (out->ps.empty())
            out->ps.assign(p, 1.0);

        std::vector<const double*> xptr(p);
        xptr[0] = x;
        for (isize_t i = 1; i < p; ++i)
            xptr[i] = xptr[i-1] + n * g[i-1];

        isize_t p1 = q;
        for (auto i : out->in)
            p1 += g[i];

        if (p1 > n) {
            eprint("ERROR: P (%td) > N (%td) is not allowed in eliminate\n", p1, n);
            return;
        }

        mat X(n, p1);
        X.leftCols(q) = mat::Map(w, n, q);

        LmFit lm;
        lm.par.rcond = par.rcond;

        int step = 0;

        while (!out->in.empty()) {
            // full model
            p1 = q;
            for (auto i : out->in) {
                X.middleCols(p1, g[i]) = mat::Map(xptr[i], n, g[i]);
                p1 += g[i];
            }

            int info = lmfit(n, p1, y, X.data(), &lm);

            if (info != 0) {
                eprint("ERROR: lmfit failed in eliminate (%s:%d): %d\n", __FILE__, __LINE__, info);
                return;
            }

            double dfe1 = lm.dfe;
            double sse1 = lm.sse;

            // reduced model
            std::vector<double> ps;

            for (auto i : out->in) {
                isize_t p0 = q;
                for (auto j : out->in) {
                    if (j != i) {
                        X.middleCols(p0, g[j]) = mat::Map(xptr[j], n, g[j]);
                        p0 += g[j];
                    }
                }

                info = lmfit(n, p0, y, X.data(), &lm);

                if (info != 0) {
                    eprint("ERROR: lmfit failed in eliminate (%s:%d): %d\n", __FILE__, __LINE__, info);
                    return;
                }

                double dfe0 = lm.dfe;
                double sse0 = lm.sse;

                double pval = 1.0;
                if (dfe1 > 0.0 && sse1 > 0.0 && dfe0 > dfe1 && sse0 > sse1) {
                    double f = (sse0 - sse1) / (dfe0 - dfe1);
                    double df = dfe1;
                    if (known_error) {
                        f /= par.mse;
                        df = par.dfe;
                    }
                    else
                        f /= (sse1 / dfe1);
                    pval = fpval(f, dfe0 - dfe1, df);
                }
                ps.push_back(pval);
                out->ps[i] = pval;
            }

            auto itr = std::max_element(ps.begin(), ps.end());
            if (*itr <= out->par.sls)
                break;

            auto i = itr - ps.begin();
            auto ii = out->in[i];
            out->in.erase(out->in.begin() + i);

            eprint("INFO: elimination %d: %td, P = %g\n", ++step, ii+1, *itr);
        }
    }

} // namespace

void forwardfit(isize_t n, isize_t q, isize_t p, const double *y, const double *w, const double *x, const isize_t *g,
                StepwiseFit *out)
{
    out->in.clear();
    out->ps.clear();

    if (out->par.thread > 0)
        select_omp(n, q, p, y, w, x, g, out);
    else
        select(n, q, p, y, w, x, g, out);
}

void backwardfit(isize_t n, isize_t q, isize_t p, const double *y, const double *w, const double *x, const isize_t *g,
                 StepwiseFit *out)
{
    out->in.resize(p);
    for (isize_t i = 0; i < p; ++i)
        out->in[i] = i;

    out->ps.clear();

    eliminate(n, q, p, y, w, x, g, out);
}

void stepwisefit(isize_t n, isize_t q, isize_t p, const double *y, const double *w, const double *x, const isize_t *g,
                 StepwiseFit *out)
{
    out->in.clear();
    out->ps.clear();

    if (out->par.thread > 0)
        select_omp(n, q, p, y, w, x, g, out);
    else
        select(n, q, p, y, w, x, g, out);

    eliminate(n, q, p, y, w, x, g, out);
}
