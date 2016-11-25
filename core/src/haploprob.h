#ifndef HAPLOPROB_H
#define HAPLOPROB_H

#include <cmath>

// Genotype / Haplotype
//      BB      Bb      bb
// -----------------------------
// AA   AABB    AABb    AAbb
//      AB/AB   AB/Ab   Ab/Ab
// -----------------------------
// Aa   AaBB    AaBb    Aabb
//      AB/aB   AB/ab   Ab/ab
//              Ab/aB
// -----------------------------
// aa   aaBB    aaBb    aabb
//      aB/aB   aB/ab   ab/ab
// -----------------------------

struct HaploProb
{
    double p_AB = 0.0;
    double p_Ab = 0.0;
    double p_aB = 0.0;
    double p_ab = 0.0;

    double tol = 1e-7;

    int n_AABB = 0;
    int n_AABb = 0;
    int n_AAbb = 0;
    int n_AaBB = 0;
    int n_AaBb = 0;
    int n_Aabb = 0;
    int n_aaBB = 0;
    int n_aaBb = 0;
    int n_aabb = 0;

    int maxit = 1000;

    void solve()
    {
        double n = n_AABB + n_AABb + n_AAbb +
                   n_AaBB + n_AaBb + n_Aabb +
                   n_aaBB + n_aaBb + n_aabb;
        n *= 2;

        const auto c_AB = (n_AABB*2 + n_AABb + n_AaBB) / n;
        const auto c_Ab = (n_AAbb*2 + n_AABb + n_Aabb) / n;
        const auto c_aB = (n_aaBB*2 + n_AaBB + n_aaBb) / n;
        const auto c_ab = (n_aabb*2 + n_aaBb + n_Aabb) / n;

        auto p_AB_ab = n_AaBb/n / 2;
        auto p_Ab_aB = p_AB_ab;

        double p_AB_ab_new = 0.0;

        for (int i = 0; i < maxit; ++i) {
            p_AB = c_AB + p_AB_ab;
            p_Ab = c_Ab + p_Ab_aB;
            p_aB = c_aB + p_Ab_aB;
            p_ab = c_ab + p_AB_ab;

            p_AB_ab_new = n_AaBb * p_AB*p_ab / (p_AB*p_ab + p_Ab*p_aB) / n;

            if (std::fabs(p_AB_ab_new - p_AB_ab) < tol)
                break;

            p_AB_ab = p_AB_ab_new;
            p_Ab_aB = n_AaBb/n - p_AB_ab;
        }
    }
};

#endif // HAPLOPROB_H
