#ifndef ANOVA_H
#define ANOVA_H

#include <memory>
#include <string>
#include <vector>


class ANOVA
{
public:

    // ANOVA Table
    //
    //   src    source names
    //   df     degrees of freedom
    //   ss     sum of squares
    //   ms     mean squares
    //   f      F value
    //   p      Pr > F
    //   error  the error info, containing dfe, sse, mse
    //   total  the total info, containing dft, sst
    //
    struct Table
    {
        std::vector<std::string> src;
        std::vector<double> df;
        std::vector<double> ss;
        std::vector<double> ms;
        std::vector<double> f;
        std::vector<double> p;
        std::vector<double> error;
        std::vector<double> total;
        std::string to_string() const;
    };

    // Parameter Estimates
    //
    //   par  parameter names
    //   est  parameter estimates
    //
    struct Solution
    {
        std::vector<std::string> par;
        std::vector<double> est;
        std::string to_string() const;
    };


public:

    // Add regression effect, X
    void add_reg(const std::string &name, const std::vector<double> &x);

    // Add main effect, A
    void add_main(const std::string &name, const std::vector<std::string> &a);

    // Add crossed effect, A*B
    void add_crossed(const std::string &name, const std::vector<std::string> &a, const std::vector<std::string> &b);

    // Add nested effect, A(B)
    void add_nested(const std::string &name, const std::vector<std::string> &a, const std::vector<std::string> &b);

    // Solve with Type I sum of squares
    Table solve1(const std::vector<double> &y) const;

    // Solve with Type III sum of squares
    Table solve3(const std::vector<double> &y) const;

    // Parameter estimates
    Solution solution(const std::vector<double> &y) const;


private:

    // Model Term
    //
    //   name   term name
    //   par    parameter names
    //   dat    design matrix
    //   contr  contrasts
    //
    struct Term
    {
        std::string name;
        std::vector<std::string> par;
        std::vector< std::vector<double> > dat;
        std::vector< std::vector<double> > contr;
    };

private:

    std::vector< std::shared_ptr<Term> > tms_;
};


#endif // ANOVA_H
