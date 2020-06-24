#ifndef LMFIT_H
#define LMFIT_H

#include "types.h"

// linear model fitting
//
// y = X*b + e
//
// y   n*1 vector
// x   n*p matrix

struct LmFit
{
    struct Parameter
    {
        double rcond = 1e-8;
    };
    Parameter par;

    std::vector<double> b;
    double dfe = 0.0;
    double sse = 0.0;
};

int lmfit(isize_t n, isize_t p, const double *y, const double *x, LmFit *out);

// linear model fitting with equality constraints
//
// y = X*b + e  and  Z*b = 0
//
// y   n*1 vector
// x   n*p matrix
// z   q*p matrix

struct LmFitEq
{
    struct Parameter
    {
        double tol = 1e-8;
        double rcond = 1e-8;
    };
    Parameter par;

    std::vector<double> b;
    double dfe = 0.0;
    double sse = 0.0;
};

int lmfit_eq(isize_t n, isize_t p, isize_t q, const double *y, const double *x, const double *z, LmFitEq *out);

// stepwise regression fitting
//
// y = W*a + X*b +e
//
// y   n*1 vector
// w   n*q covariate matrix, including intercept
// x   n*sum(g) predictor matrix
// g   p*1 vector, number of columns of each predictor

struct StepwiseFit
{
    struct Parameter
    {
        double sle = 0.05;  // the significance level for entry into the model
        double sls = 0.05;  // the significance level for staying in the model
        double rsq = 0.95;  // maximum model R^2
        double dfe = 0.0;   // error df
        double mse = 0.0;   // error mean squares
        double rcond = 1e-8;
        int mtc = 0;        // multiple testing correction, 0=NONE, 1=BON, 2=FDR, 3=HOLM
        int thread = 0;
    };
    Parameter par;

    std::vector<double> ps;     // m*1 p-value vector
    std::vector<isize_t> in;    // index of predictors in model
};

void forwardfit (isize_t n, isize_t q, isize_t p, const double *y, const double *w, const double *x, const isize_t *g,
                 StepwiseFit *out);

void backwardfit(isize_t n, isize_t q, isize_t p, const double *y, const double *w, const double *x, const isize_t *g,
                 StepwiseFit *out);

void stepwisefit(isize_t n, isize_t q, isize_t p, const double *y, const double *w, const double *x, const isize_t *g,
                 StepwiseFit *out);

#endif // LMFIT_H
