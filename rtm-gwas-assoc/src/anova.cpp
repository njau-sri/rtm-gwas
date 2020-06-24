#include "anova.h"

#include "print.h"
#include "statutil.h"
#include "vectorutil.h"
#include "matrix.h"
#include "lmfit.h"

struct ANOVA::Term
{
    std::string name;                // term name
    std::vector<std::string> param;  // parameter names
    mat dat;                         // design matrix
    mat constr;                      // sum to zero constraints
    mat constr_wt;                   // weighted sum to zero constraints
};

ANOVA::ANOVA() = default;
ANOVA::~ANOVA() = default;

void ANOVA::add_reg(const std::string &name, const std::vector<double> &x)
{
    tms_.emplace_back(new Term);
    auto &ptr = tms_.back();

    ptr->name = name;
    ptr->param.push_back(name);
    ptr->dat = mat::Map(x.data(), length(x), 1);
}

void ANOVA::add_main(const std::string &name, const std::vector<std::string> &a)
{
    tms_.emplace_back(new Term);
    auto &ptr = tms_.back();

    std::vector<isize_t> gi;
    std::vector<std::string> gn;
    grpidx(a, gi, gn);

    design3(length(gi), gi.data(), ptr->dat);

    ptr->name = name;
    for (auto &e : gn)
        ptr->param.push_back(name + " " + e);

    ptr->constr.setOnes(1, length(gn));
    ptr->constr_wt = ptr->dat.colwise().mean();
}

void ANOVA::add_crossed(const std::string &name, const std::vector<std::string> &a, const std::vector<std::string> &b)
{
    tms_.emplace_back(new Term);
    auto &ptr = tms_.back();

    std::vector<isize_t> gia;
    std::vector<std::string> gna;
    grpidx(a, gia, gna);

    std::vector<isize_t> gib;
    std::vector<std::string> gnb;
    grpidx(b, gib, gnb);

    mat xa, xb;
    design3(length(gia), gia.data(), xa);
    design3(length(gib), gib.data(), xb);

    ptr->name = name;

    isize_t na = length(gna);
    isize_t nb = length(gnb);

    ptr->dat.resize(length(a), na*nb);
    ptr->constr.setZero(na+nb, na*nb);
    ptr->constr_wt.setZero(na+nb, na*nb);

    // A * B : A1B1 A1B2 A1B3 A2B1 A2B2 A2B3
    for (isize_t i = 0, k = 0; i < na; ++i) {
        for (isize_t j = 0; j < nb; ++j) {
            ptr->param.push_back(name + " " + gna[i] + "*" + gnb[j]);
            ptr->dat.col(k) = xa.col(i).cwiseProduct(xb.col(j));
            ++k;
        }
    }

    isize_t k = 0;
    mat z = ptr->dat.colwise().mean();

    // constraints for A
    for (isize_t j = 0; j < nb; ++j) {
        for (isize_t i = 0; i < na; ++i)
            ptr->constr(k, i*nb+j) = 1.0;
        for (isize_t i = 0; i < na; ++i)
            ptr->constr_wt(k, i*nb+j) = z(i*nb+j);
        ++k;
    }

    // constraints for B
    for (isize_t i = 0; i < na; ++i) {
        ptr->constr.row(k).segment(i*nb, nb).setOnes();
        ptr->constr_wt.row(k).segment(i*nb, nb) = z;
        ++k;
    }
}

void ANOVA::add_nested(const std::string &name, const std::vector<std::string> &a, const std::vector<std::string> &b)
{
    tms_.emplace_back(new Term);
    auto &ptr = tms_.back();

    std::vector<isize_t> gia;
    std::vector<std::string> gna;
    grpidx(a, gia, gna);

    std::vector<isize_t> gib;
    std::vector<std::string> gnb;
    grpidx(b, gib, gnb);

    mat xa, xb;
    design3(length(gia), gia.data(), xa);
    design3(length(gib), gib.data(), xb);

    ptr->name = name;

    isize_t na = length(gna);
    isize_t nb = length(gnb);

    ptr->dat.resize(length(a), na*nb);
    ptr->constr.setZero(nb, na*nb);
    ptr->constr_wt.setZero(nb, na*nb);

    // A (B) : A1B1 A2B1 A3B1 A1B2 A2B2 A3B2
    for (isize_t j = 0, k = 0; j < nb; ++j) {
        for (isize_t i = 0; i < na; ++i) {
            ptr->param.push_back(name + " " + gna[i] + "(" + gnb[j] + ")");
            ptr->dat.col(k) = xa.col(i).cwiseProduct(xb.col(j));
            ++k;
        }
    }

    mat z = ptr->dat.colwise().mean();

    // constraints for A
    for (isize_t j = 0; j < nb; ++j) {
        ptr->constr.row(j).segment(j*na, na).setOnes();
        ptr->constr_wt.row(j).segment(j*na, na) = z;
    }

    std::vector<isize_t> idx;
    isize_t n = ncol(ptr->dat);
    for (isize_t j = 0; j < n; ++j)
        if ((ptr->dat.col(j).array() != 0.0).any())
            idx.push_back(j);

    isize_t n1 = length(idx);
    if (n1 < n) {
        for (isize_t j = 0; j < n1; ++j) {
            isize_t jj = idx[j];
            if (jj != j) {
                ptr->dat.col(j) = ptr->dat.col(jj);
                ptr->constr.col(j) = ptr->constr.col(jj);
                ptr->constr_wt.col(j) = ptr->constr_wt.col(jj);
            }
        }
        ptr->dat.conservativeResize(Eigen::NoChange, n1);
        ptr->constr.conservativeResize(Eigen::NoChange, n1);
        ptr->constr_wt.conservativeResize(Eigen::NoChange, n1);
    }

    subset(ptr->param, idx).swap(ptr->param);
}

ANOVA::Table ANOVA::solve1(const std::vector<double> &y) const
{
    Table tbl;

    isize_t n = length(y);
    
    double dfe = n - 1.0;
    double sse = css(n, y.data());

    tbl.total.push_back(dfe);
    tbl.total.push_back(sse);

    isize_t p = 1;
    for (auto &ptr : tms_)
        p += ncol(ptr->dat);

    mat x = mat::Zero(n, p);
    x.col(0).setOnes();

    p = 1;
    for (auto &ptr : tms_) {
        x.middleCols(p, ncol(ptr->dat)) = ptr->dat;
        p += ncol(ptr->dat);

        LmFit lm;
        lmfit(n, p, y.data(), x.data(), &lm);

        tbl.src.push_back(ptr->name);
        tbl.df.push_back(dfe - lm.dfe);
        tbl.ss.push_back(sse - lm.sse);

        dfe = lm.dfe;
        sse = lm.sse;
    }

    tbl.error.push_back(dfe);
    tbl.error.push_back(sse);
    tbl.error.push_back(sse / dfe);

    isize_t nterms = length(tms_);
    tbl.ms.resize(nterms);
    tbl.f.resize(nterms);
    tbl.p.assign(nterms, std::numeric_limits<double>::quiet_NaN());
    for (isize_t j = 0; j < nterms; ++j) {
        tbl.ms[j] = tbl.ss[j] / tbl.df[j];
        tbl.f[j] = tbl.ms[j] / tbl.error[2];
        if (tbl.df[j] > 0.0 && tbl.ss[j] > 0.0)
            tbl.p[j] = fpval(tbl.f[j], tbl.df[j], tbl.error[0]);
    }

    return tbl;
}

ANOVA::Table ANOVA::solve3(const std::vector<double> &y) const
{
    Table tbl;

    isize_t n = length(y);

    double dft = n - 1.0;
    double sst = css(n, y.data());

    tbl.total.push_back(dft);
    tbl.total.push_back(sst);

    isize_t p = 1;
    isize_t q = 0;
    for (auto &ptr : tms_) {
        p += ncol(ptr->dat);
        q += nrow(ptr->constr);
    }

    mat x(n, p);
    x.col(0).setOnes();
    p = 1;
    for (auto &ptr : tms_) {
        x.middleCols(p, ncol(ptr->dat)) = ptr->dat;
        p += ncol(ptr->dat);
    }

    mat z = mat::Zero(q, p);
    isize_t r = 0, c = 1;
    for (auto &ptr : tms_) {
        isize_t nr = nrow(ptr->constr);
        isize_t nc = ncol(ptr->constr);
        if (nr > 0) {
            z.block(r, c, nr, nc) = ptr->constr;
            r += nr;
            c += nc;
        }
        else
            c += ncol(ptr->dat);
    }

    LmFitEq lm;
    lmfit_eq(n, p, q, y.data(), x.data(), z.data(), &lm);

    tbl.error.push_back(lm.dfe);
    tbl.error.push_back(lm.sse);
    tbl.error.push_back(lm.sse / lm.dfe);

    isize_t nterms = length(tms_);

    if (nterms == 1) {
        tbl.src.push_back(tms_[0]->name);
        tbl.df.push_back(dft - lm.dfe);
        tbl.ss.push_back(sst - lm.sse);
    }
    else {
        for (isize_t i = 0; i < nterms; ++i) {
            isize_t p0 = p - ncol(tms_[i]->dat);
            mat x0(n, p0);

            isize_t nleft = 1;
            isize_t nright = 0;
            for (isize_t j = 0; j < nterms; ++j) {
                if (j < i)
                    nleft += ncol(tms_[j]->dat);
                else if (j > i)
                    nright += ncol(tms_[j]->dat);
            }

            x0.leftCols(nleft) = x.leftCols(nleft);
            if (nright > 0)
                x0.rightCols(nright) = x.rightCols(nright);

            isize_t q0 = q - nrow(tms_[i]->constr);
            mat z0 = mat::Zero(q0, p0);

            r = 0, c = 1;
            for (isize_t j = 0; j < nterms; ++j) {
                if (j == i)
                    continue;
                isize_t nr = nrow(tms_[j]->constr);
                isize_t nc = ncol(tms_[j]->constr);
                if (nr > 0) {
                    z0.block(r, c, nr, nc) = tms_[j]->constr;
                    r += nr;
                    c += nc;
                }
                else
                    c += ncol(tms_[j]->dat);
            }

            LmFitEq lm0;
            lmfit_eq(n, p0, q0, y.data(), x0.data(), z0.data(), &lm0);

            tbl.src.push_back(tms_[i]->name);
            tbl.df.push_back(lm0.dfe - lm.dfe);
            tbl.ss.push_back(lm0.sse - lm.sse);
        }
    }

    tbl.ms.resize(nterms);
    tbl.f.resize(nterms);
    tbl.p.assign(nterms, std::numeric_limits<double>::quiet_NaN());
    for (isize_t j = 0; j < nterms; ++j) {
        tbl.ms[j] = tbl.ss[j] / tbl.df[j];
        tbl.f[j] = tbl.ms[j] / tbl.error[2];
        if (tbl.df[j] > 0.0 && tbl.ss[j] > 0.0)
            tbl.p[j] = fpval(tbl.f[j], tbl.df[j], tbl.error[0]);
    }

    return tbl;
}

ANOVA::Solution ANOVA::solution(const std::vector<double> &y) const
{
    isize_t n = length(y);

    isize_t p = 1;
    isize_t q = 0;
    for (auto &ptr : tms_) {
        p += ncol(ptr->dat);
        q += nrow(ptr->constr);
    }

    mat x(n, p);
    std::vector<std::string> param;

    x.col(0).setOnes();
    param.push_back("Constant");

    p = 1;
    for (auto &ptr : tms_) {
        param.insert(param.end(), ptr->param.begin(), ptr->param.end());
        x.middleCols(p, ncol(ptr->dat)) = ptr->dat;
        p += ncol(ptr->dat);
    }

    mat z = mat::Zero(q, p);
    isize_t r = 0, c = 1;
    for (auto &ptr : tms_) {
        isize_t nr = nrow(ptr->constr);
        isize_t nc = ncol(ptr->constr);
        if (nr > 0) {
            z.block(r, c, nr, nc) = ptr->constr;
            r += nr;
            c += nc;
        }
        else
            c += ncol(ptr->dat);
    }

    LmFitEq lm;
    lmfit_eq(n, p, q, y.data(), x.data(), z.data(), &lm);

    Solution sol;
    sol.par.swap(param);
    sol.est = lm.b;

    return sol;
}

ANOVA::Solution ANOVA::solution_wtsum(const std::vector<double> &y) const
{
    isize_t n = length(y);

    isize_t p = 1;
    isize_t q = 0;
    for (auto &ptr : tms_) {
        p += ncol(ptr->dat);
        q += nrow(ptr->constr_wt);
    }

    mat x(n, p);
    std::vector<std::string> param;

    x.col(0).setOnes();
    param.push_back("Constant");

    p = 1;
    for (auto &ptr : tms_) {
        param.insert(param.end(), ptr->param.begin(), ptr->param.end());
        x.middleCols(p, ncol(ptr->dat)) = ptr->dat;
        p += ncol(ptr->dat);
    }

    mat z = mat::Zero(q, p);
    isize_t r = 0, c = 1;
    for (auto &ptr : tms_) {
        isize_t nr = nrow(ptr->constr_wt);
        isize_t nc = ncol(ptr->constr_wt);
        if (nr > 0) {
            z.block(r, c, nr, nc) = ptr->constr_wt;
            r += nr;
            c += nc;
        }
        else
            c += ncol(ptr->dat);
    }

    LmFitEq lm;
    lmfit_eq(n, p, q, y.data(), x.data(), z.data(), &lm);

    Solution sol;
    sol.par.swap(param);
    sol.est = lm.b;

    return sol;
}

std::string ANOVA::Table::to_string() const
{
    std::string s;

    s += "Source\tDF\tSS\tMS\tF\tp\n";

    isize_t n = length(src);
    for (isize_t i = 0; i < n; ++i)
        s += sprint("%s\t%g\t%g\t%g\t%g\t%g\n", src[i], df[i], ss[i], ms[i], f[i], p[i]);

    s += "Error";
    for (auto &e : error)
        s += sprint("\t%g", e);
    s += '\n';

    s += "Total";
    for (auto &e : total)
        s += sprint("\t%g", e);
    s += '\n';

    return s;
}

std::string ANOVA::Solution::to_string() const
{
    std::string s;

    s += "Parameter\tEstimate\n";

    isize_t n = length(par);
    for (isize_t i = 0; i < n; ++i)
        s += sprint("%s\t%g\n", par[i], est[i]);

    return s;
}
