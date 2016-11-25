#include <memory>
#include <iostream>
#include "assoclmm.h"
#include "emma.h"
#include "util.h"
#include "stat.h"

namespace
{

void design(const vector<allele_t> &a, const vector<allele_t> &b, vector< vector<double> > &x)
{
    vector<int> group;

    if ( b.empty() || a == b ) {
        factor(a).swap(group);
    }
    else {
        vector< std::pair<allele_t,allele_t> > ab;
        auto n = a.size();
        for (size_t i = 0; i < n; ++i) {
            if (a[i] < b[i])
                ab.emplace_back(a[i],b[i]);
            else
                ab.emplace_back(b[i],a[i]);
        }
        factor(ab).swap(group);
    }

    ::design(2, group, x);
}

} // namespace

// Kang, H.M. et al. Variance component model to account for sample structure in genome-wide association studies.
//   Nat Genet 42, 348-54 (2010). doi: 10.1038/ng.548

void assoc_LMM(bool haploid, const vector< vector<allele_t> > &geno, const vector<size_t> &obs,
               const vector<double> &pheno,
               const vector< vector<double> > &addcov,
               const vector< vector<double> > &intcov,
               const vector< vector<double> > &kin,
               vector<double> &result)
{
    auto nobs = pheno.size();
    auto nmar = geno.size();

    vector< vector<double> > X0;
    X0.push_back(vector<double>(nobs,1.0));
    for (auto &e : addcov)
        X0.push_back(e);

    vector< vector<double> > K;
    for (auto &e : subset(kin,obs))
        K.push_back(subset(e,obs));

    auto emma = std::make_shared<EMMA>();
    emma->solve(pheno, X0, K);

    auto vc = emma->varcomp();
    std::cerr << "INFO: REML=" << vc.REML << ", delta=" << vc.delta << ", Vg=" << vc.vg << ", Ve=" << vc.ve << "\n";

    result.assign(2*nmar, 0.0);

    for (size_t j = 0; j < nmar; ++j) {
        result[nmar+j] = 1.0;

        vector<allele_t> ga, gb;
        vector<size_t> vobs;

        if ( haploid ) {
            for (size_t i = 0; i < nobs; ++i) {
                auto a = geno[j][obs[i]];
                if (a != 0) {
                    vobs.push_back(i);
                    ga.push_back(a);
                }
            }
        }
        else {
            for (size_t i = 0; i < nobs; ++i) {
                auto k = obs[i];
                auto a = geno[j][2*k];
                auto b = geno[j][2*k+1];
                if (a != 0 && b != 0) {
                    vobs.push_back(i);
                    ga.push_back(a);
                    gb.push_back(b);
                }
            }
        }

        if ( vobs.empty() )
            continue;

        vector< vector<double> > xg;
        design(ga, gb, xg);
        if (xg.empty())
            continue;

        auto nvobs = vobs.size();

        vector< vector<double> > X;
        X.push_back(vector<double>(nvobs,1.0));

        if (nvobs == nobs)
            std::copy(addcov.begin(), addcov.end(), std::back_inserter(X));
        else
            for (auto &e : addcov)
                X.push_back(subset(e,vobs));

        std::copy(xg.begin(), xg.end(), std::back_inserter(X));

        for (auto &eg : xg) {
            for (auto ei : intcov) {
                if (nvobs != nobs)
                    subset(ei,vobs).swap(ei);
                std::transform(eg.begin(), eg.end(), ei.begin(), ei.begin(), std::multiplies<double>());
                X.push_back(ei);
            }
        }

        auto p = xg.size() + xg.size() * intcov.size();
        double fval, pval;

        if (nvobs != nobs) {
            vector< vector<double> > Kc;
            for (auto &e : subset(K,vobs))
                Kc.push_back(subset(e,vobs));
            emma->ftest(p, subset(pheno,vobs), X, Kc, fval, pval);
        }
        else
            emma->ftest(p, X, fval, pval);

        result[j] = p;
        result[nmar+j] = pval;
    }
}
