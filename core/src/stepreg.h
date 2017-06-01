#ifndef STEPREG_H
#define STEPREG_H

#include "main.h"

class StepReg
{
public:
    struct Params
    {
        double maxrsq = 0.99;
        double sle = 0.05;
        double sls = 0.05;
        int maxstep = 0;
        int multtest = 0;  // 1:BON, 2:FDR
    };

    void set_params(const Params &par) { m_par = par; }

    void init_model(const vector<size_t> &init) { m_model = init; }

    const vector<size_t>& model() const { return m_model; }

    const vector<double>& ps() const { return m_ps; }

    void add_covariate(const vector<double> &v)
    {
        m_cofs.insert(m_cofs.end(), v.begin(), v.end());
    }

    void add_effect(const vector<double> &v)
    {
        m_effs.insert(m_effs.end(), v.begin(), v.end());
        m_cols.push_back(1);
    }

    void add_effect(const vector< vector<double> > &x)
    {
        for (auto &v : x)
            m_effs.insert(m_effs.end(), v.begin(), v.end());
        m_cols.push_back(x.size());
    }

    void forward(const vector<double> &y);

    void backward(const vector<double> &y);

    void stepwise(const vector<double> &y);

private:
    Params m_par;
    vector<double> m_ps;
    vector<double> m_cofs;
    vector<double> m_effs;
    vector<size_t> m_cols;
    vector<size_t> m_model;
};

#endif // STEPREG_H
