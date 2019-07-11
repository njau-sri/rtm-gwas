#ifndef STEPREG_H
#define STEPREG_H


#include <vector>


// Stepwise Regression
//
//   Option
//     sle   the significance level for entry into the model
//     sls   the significance level for staying in the model
//     rsq   maximum model R square
//     mtc   multiple testing correction, 0=NONE;1=BON;2=FDR;3=HOLM
//
//   Input
//     x0    n*p covariates matrix, including intercept
//     x1    n*[sum(c1)] design matrix x1, m independent variables with varying degrees of freedom
//     c1    m*1 vector, containing the number of columns for design matrix of each independent variable
//     y     n*1 dependent variable y
//
//   Input (Backward)
//     ps   m*1 vector p-values
//     in   index of independent variables for initial model
//
//   Output
//     ps    m*1 p-value vector
//     in    index of selected independent variables
//
//   Output (Backward)
//     ps   m*1 vector p-values, updated
//     in   index of selected independent variables
//
struct StepReg
{
    double sle = 0.05;
    double sls = 0.05;
    double rsq = 0.95;
    int mtc = 0;


    void forward(const std::vector<double> &x0, const std::vector<double> &x1, const std::vector<std::size_t> &c1,
                 const std::vector<double> &y, std::vector<double> &ps, std::vector<std::size_t> &in);


    void forward_omp(const std::vector<double> &x0, const std::vector<double> &x1, const std::vector<std::size_t> &c1,
                     const std::vector<double> &y, std::vector<double> &ps, std::vector<std::size_t> &in);


    void backward(const std::vector<double> &x0, const std::vector<double> &x1, const std::vector<std::size_t> &c1,
                  const std::vector<double> &y, std::vector<double> &ps, std::vector<std::size_t> &in);
};


#endif // STEPREG_H
