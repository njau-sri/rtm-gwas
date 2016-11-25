#ifndef GLM_H
#define GLM_H

#include <memory>
#include "main.h"

class GLM
{
public:
    struct AnovaTable
    {
        vector<string> names;
        vector<double> df;
        vector<double> ss;
        vector<double> ms;
        vector<double> f;
        vector<double> p;
        vector<double> error;
        vector<double> total;
    };

    struct Coefficient
    {
        vector<string> params;
        vector<double> coeffs;
    };

    // Regression effect, X
    void add_reg(const string &name, const vector<double> &x);

    // Main effect, A
    void add_main(const string &name, const vector<string> &a);

    // Crossed effect, A*B
    void add_crossed(const string &name, const vector<string> &a, const vector<string> &b);

    // Nested effect, A(B)
    void add_nested(const string &name, const vector<string> &a, const vector<string> &b);

    // Continuous-by-Class effect, X*A
    void add_continuous_by_class(const string &name, const vector<double> &x, const vector<string> &a);

    // Continuous-Nesting-Class effect, X(A)
    void add_continuous_nesting_class(const string &name, vector<double> &x, const vector<string> &a);

    AnovaTable solve1(const vector<double> &y) const;

    AnovaTable solve3(const vector<double> &y) const;

    Coefficient coeff(const vector<double> &y) const;

    void reset();

private:
    struct ModelTerm
    {
        string name;
        vector<double> coeffs;
        vector<string> coeffnames;
        vector< vector<double> > data;
        vector< vector<double> > constr;
    };

    vector< std::shared_ptr<ModelTerm> > m_terms;
};

#endif // GLM_H
