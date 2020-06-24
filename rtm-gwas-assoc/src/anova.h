#ifndef ANOVA_H
#define ANOVA_H

#include <memory>
#include <string>
#include <vector>

class ANOVA
{
public:

    struct Table
    {
        std::vector<std::string> src;   // model term names
        std::vector<double> df;         // degrees of freedom
        std::vector<double> ss;         // sum of squares
        std::vector<double> ms;         // mean squares
        std::vector<double> f;          // F value
        std::vector<double> p;          // Pr > F
        std::vector<double> error;      // dfe, sse, mse
        std::vector<double> total;      // dfe, sst
        std::string to_string() const;
    };

    struct Solution
    {
        std::vector<std::string> par;   // names
        std::vector<double> est;        // estimates
        std::string to_string() const;
    };

public:
    // std::unique_ptr with incomplete type
    // https://en.cppreference.com/w/cpp/language/pimpl
    // https://stackoverflow.com/questions/9954518/stdunique-ptr-with-an-incomplete-type-wont-compile
    ANOVA();
    ~ANOVA();

    ANOVA(ANOVA&&) = delete;
    ANOVA(const ANOVA&) = delete;
    ANOVA& operator=(ANOVA&&) = delete;
    ANOVA& operator=(const ANOVA&) = delete;

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

    Solution solution_wtsum(const std::vector<double> &y) const;

private:
    struct Term;
    std::vector< std::unique_ptr<Term> > tms_;
};

#endif // ANOVA_H
