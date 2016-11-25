#include <memory>
#include <utility>
#include "assoclm.h"
#include "lstsqr.h"
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

void assoc_LM(bool haploid, const vector< vector<allele_t> > &geno,  const vector<size_t> &obs,
              const vector<double> &pheno,
              const vector< vector<double> > &addcov,
              const vector< vector<double> > &intcov,
              vector<double> &result)
{
    auto nobs = pheno.size();
    auto nmar = geno.size();

    if ( intcov.empty() ) {
        result.assign(3*nmar, 0.0);
        for (size_t j = 0; j < nmar; ++j)
            result[nmar+j] = 1.0;
    }
    else {
        result.assign(6*nmar, 0.0);
        for (size_t j = 0; j < nmar; ++j)
            result[nmar+j] = result[4*nmar+j] = 1.0;
    }

    auto lstsqr = std::make_shared<LstSqr>();

    vector<size_t> vobs;
    vector<allele_t> ga, gb;
    vector< vector<double> > xg;
    vector<double> y, x;

    auto sst = sumsqc(pheno);

    size_t q0 = 1 + addcov.size();
    vector<double> x0(nobs,1.0);
    for (auto &e : addcov)
        x0.insert(x0.end(), e.begin(), e.end());

    lstsqr->solve(q0, pheno, x0);

    auto dfe0 = lstsqr->stats().dfe;
    auto sse0 = lstsqr->stats().sse;

    for (size_t j = 0; j < nmar; ++j) {
        vobs.clear();
        ga.clear();
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
            gb.clear();
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

        design(ga, gb, xg);
        if ( xg.empty() )
            continue;

        auto nvobs = vobs.size();

        auto jsst = sst;
        auto jdfe0 = dfe0;
        auto jsse0 = sse0;

        if (nvobs != nobs) {
            y = subset(pheno, vobs);
            jsst = sumsqc(y);
            x.assign(nvobs,1.0);
            for (auto &e : addcov) {
                auto v = subset(e, vobs);
                x.insert(x.end(), v.begin(), v.end());
            }
            lstsqr->solve(q0, y, x);
            jdfe0 = lstsqr->stats().dfe;
            jsse0 = lstsqr->stats().sse;
        }
        else {
            y = pheno;
            x = x0;
        }

        auto q1 = q0 + xg.size();
        for (auto &e : xg)
            x.insert(x.end(), e.begin(), e.end());

        lstsqr->solve(q1, y, x);
        auto jdfe1 = lstsqr->stats().dfe;
        auto jsse1 = lstsqr->stats().sse;

        auto dfx = jdfe0 - jdfe1;
        auto ssx = jsse0 - jsse1;
        auto dfe = jdfe1;
        auto sse = jsse1;

        if ( ! intcov.empty() ) {
            auto q2 = q1 + xg.size()*intcov.size();
            for (auto &e1 : xg) {
                for (auto e2 : intcov) {
                    if (nvobs != nobs) subset(e2, vobs).swap(e2);
                    std::transform(e1.begin(), e1.end(), e2.begin(), e2.begin(), std::multiplies<double>());
                    x.insert(x.end(), e2.begin(), e2.end());
                }
            }

            lstsqr->solve(q2, y, x);
            dfe = lstsqr->stats().dfe;
            sse = lstsqr->stats().sse;

            auto dfxi = jdfe1 - dfe;
            auto ssxi = jsse1 - sse;

            auto f = (ssxi/dfxi) / (sse/dfe);
            result[3*nmar+j] = dfxi;
            result[4*nmar+j] = fcdf(f, dfxi, dfe, false);
            result[5*nmar+j] = ssxi / jsst;
        }

        if (dfx > 0.0 && ssx > 0.0) {
            auto f = (ssx/dfx) / (sse/dfe);
            result[j] = dfx;
            result[nmar+j] = fcdf(f, dfx, dfe, false);
            result[2*nmar+j] = ssx / jsst;
        }
    }
}
