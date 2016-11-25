#ifndef EMMA_H
#define EMMA_H

#include "main.h"

class EMMA
{
public:
    struct Params
    {
        double tol = 1e-8;
        double llim = -10.0;
        double ulim = 10.0;
        int ngrids = 100;
        int maxit = 100;
    };

    struct VarComp
    {
        double REML = 0.0;
        double delta = 0.0;
        double vg = 0.0;
        double ve = 0.0;
    };

    const VarComp& varcomp() const { return m_vc; }

    const Params& params() const { return m_par; }

    void set_params(const Params &opt) { m_par = opt; }

    void solve(const vector<double> &y, const vector< vector<double> > &X, const vector< vector<double> > &K);

    void ftest(size_t p, const vector< vector<double> > &X, double &fval, double &pval);

    void ftest(size_t p, const vector<double> &y, const vector< vector<double> > &X,
               const vector< vector<double> > &K, double &fval, double &pval);

private:
    void compute();

    int  cholfact();

private:
    Params m_par;
    VarComp m_vc;

    vector<double> m_y;
    vector<double> m_X;
    vector<double> m_K;

    vector<double> m_L;
    vector<double> m_Ly;
    vector<double> m_LX;

    size_t m_n = 0;
    size_t m_q = 0;
};

#endif // EMMA_H
