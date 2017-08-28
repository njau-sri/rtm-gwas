#include <limits>
#include <memory>
#include <utility>
#include <algorithm>
#include "assocrtm.h"
#include "stepreg.h"
#include "stat.h"

void assoc_RTM(int multtest, int maxstep, double alpha, double maxrsq,
               bool haploid, const vector< vector<allele_t> > &geno, const vector<size_t> &obs,
               const vector<double> &pheno,
               const vector< vector<double> > &addcov,
               const vector< vector<double> > &intcov,
               vector<size_t> &result, vector<double> &ps)
{
    auto stepreg = std::make_shared<StepReg>();

    StepReg::Params par;
    par.multtest = multtest;
    par.maxstep = maxstep;
    par.maxrsq = maxrsq;
    par.sle = alpha;
    par.sls = alpha;

    stepreg->set_params(par);

    for (auto &e : addcov)
        stepreg->add_covariate(e);

    auto m = geno.size();
    vector<allele_t> ga;
    vector<std::pair<allele_t,allele_t> > gab;
    vector<size_t> idx;

    for (size_t j = 0; j < m; ++j) {
        vector< vector<double> > xg;

        // complete genotype
        if ( haploid ) {
            ga.clear();
            for (auto i : obs)
                ga.push_back(geno[j][i]);
            design(1, factor(ga), xg);
        }
        else {
            gab.clear();
            for (auto i : obs) {
                auto a = geno[j][2*i], b = geno[j][2*i+1];
                if (a > b)
                    std::swap(a, b);
                gab.emplace_back(a, b);
            }
            design(1, factor(gab), xg);
        }

        if ( xg.empty() )
            continue;

        if ( ! intcov.empty() ) {
            auto ng = xg.size();
            for (size_t k = 0; k < ng; ++k) {
                for (auto v : intcov) {
                    std::transform(v.begin(), v.end(), xg[k].begin(), v.begin(), std::multiplies<double>());
                    xg.push_back(v);
                }
            }
        }

        stepreg->add_effect(xg);

        idx.push_back(j);
    }

    stepreg->stepwise(pheno);

    result = stepreg->model();
    subset(idx,result).swap(result);

    ps.assign(m, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < idx.size(); ++i)
        ps[idx[i]] = stepreg->ps()[i];
}
