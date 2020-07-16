#include "popgen.h"

#include <cmath>
#include <algorithm>

#include "print.h"

// Genotype
//         BB     Bb     bb
//   AA   AABB   AABb   AAbb      g[0][0]  g[0][1]  g[0][2]
//   Aa   AaBB   AaBb   Aabb      g[1][0]  g[1][1]  g[1][2]
//   aa   aaBB   aaBb   aabb      g[2][0]  g[2][1]  g[2][2]

//          BB       Bb       bb
//   AA   AB//AB   AB//Ab   Ab//Ab
//   Aa   AB//aB            Ab//ab
//   aa   aB//aB   aB//ab   ab//ab

// AaBb = AB//ab + Ab//aB

// David Siegmund & Benjamin Yakir. 2007. The Statistics of Gene Mapping. Springer. doi: 10.1007/978-0-387-49686-3
//   Pages 307-322: Inferring Haplotypes from Genotypes and Testing Association

// x = p1*p4 / (p1*p4 + p2*p3);

// AB    = 2*AABB + AABb + AaBB
// Ab    = 2*AAbb + AABb + Aabb
// aB    = 2*aaBB + AaBB + aaBb
// ab    = 2*aabb + aaBb + Aabb
// AaBb  = AaBb

double est_het(int maxit, double tol, int AB, int Ab, int aB, int ab, int AaBb)
{
    double x = 0.5;
    double conv = 0.0;

    for (int i = 0; i < maxit; ++i) {
        double y1 = AaBb * x;
        double y2 = AaBb - y1;

        double n1 = AB + y1;
        double n2 = Ab + y2;
        double n3 = aB + y2;
        double n4 = ab + y1;

        double x_new = n1 * n4 / (n1 * n4 + n2 * n3);

        conv = std::fabs(x - x_new);
        if (conv < tol)
            return x_new;

        x = x_new;
    }

    eprint("WARNING: maximum number of iterations exceeded in est_het: "
           "maxit = %d, tol = %g, conv = %g\n", maxit, tol, conv);

    return x;
}

void est_hap_freq(int maxit, double tol, int g[3][3], double p[4])
{
    int AB = 2 * g[0][0] + g[0][1] + g[1][0];
    int Ab = 2 * g[0][2] + g[0][1] + g[1][2];
    int aB = 2 * g[2][0] + g[1][0] + g[2][1];
    int ab = 2 * g[2][2] + g[2][1] + g[1][2];

    int AaBb = g[1][1];

    double x = AaBb > 0 ? est_het(maxit, tol, AB, Ab, aB, ab, AaBb) : 0.0;

    double y1 = AaBb * x;
    double y2 = AaBb - y1;

    int n = AB + Ab + aB + ab + 2 * AaBb;

    p[0] = (AB + y1) / n;
    p[1] = (Ab + y2) / n;
    p[2] = (aB + y2) / n;
    p[3] = (ab + y1) / n;
}

void calc_dprime_rsq(double pa, double pb, double pab, double &dprime, double &rsq)
{
    double qa = 1.0 - pa;
    double qb = 1.0 - pb;

    double d = pab - pa * pb;

    // Lewontin (1964). Genetics 49:49-67
    double dmax = d > 0 ? std::min(pa * qb, qa * pb) : std::min(pa * pb, qa * qb);
    dprime = d / dmax;

    // Hill & Robertson (1968). Theor Appl Genet 38:226-231
    rsq = d * d / (pa * qa * pb * qb);
}

// 王建康, 李慧慧, 张鲁燕. 2014. 基因定位与育种设计. 北京: 科学出版社
//   第 38-86 页: 连锁分析和遗传图谱构建
//
// F2: P1(AABB) * P2(aabb) -> F1(AaBb) -> selfing for one generation
//
// Gamete        AB (1-r)/2   Ab r/2     aB r/2     ab (1-r)/2
// AB (1-r)/2    (1-r)^2/4    r(1-r)/4   r(1-r)/4   (1-r)^2/4
// Ab  r/2        r(1-r)/4     r^2/4      r^2/4      r(1-r)/4
// aB  r/2        r(1-r)/4     r^2/4      r^2/4      r(1-r)/4
// ab (1-r)/2    (1-r)^2/4    r(1-r)/4   r(1-r)/4   (1-r)^2/4
//
//
// Genotype      BB              Bb                  bb
// ----------------------------------------------------------------
// AA            AB//AB          AB//Ab,Ab//AB       Ab/Ab
//               (1-r)^2/4       r(1-r)/2            r^2/4
// ----------------------------------------------------------------
// Aa            AB//aB,aB//AB   AB//ab,ab//AB       Ab//ab,ab//Ab
//                               Ab//aB,aB//Ab
//               r(1-r)/2        (1-r)^2/2 + r^2/2   r(1-r)/2
// ----------------------------------------------------------------
// aa            aB//aB          aB//ab,ab//aB       ab//ab
//               r^2/4           r(1-r)/2            (1-r)^2/4
// ----------------------------------------------------------------
//
//         BB     Bb     bb
//   AA   AABB   AABb   AAbb      g[0][0]  g[0][1]  g[0][2]
//   Aa   AaBB   AaBb   Aabb      g[1][0]  g[1][1]  g[1][2]
//   aa   aaBB   aaBb   aabb      g[2][0]  g[2][1]  g[2][2]
//
// one recombinant gamete (AABb AaBB Aabb aaBb)
// m1 = g[0][1] + g[1][0] + g[1][2] + g[2][1]
//
// two recombinant gametes (AAbb aaBB)
// m2 = g[0][2] + g[2][0]
//
// double heterozygote (AaBb)
// zero recombinant gamete: AB//ab
// two recombinant gametes: Ab//aB
// E-step: p = (r^2/2) / [(1-r)^2/2 + r^2/2]
//
// recombination fraction
// M-step: r = (m1 + 2*m2 + p*2*g[1][1]) / 2n

double est_rec_frac_f2(int maxit, double tol, int g[3][3])
{
    // number of recombinant gametes
    int m = g[0][1] + g[1][0] + g[1][2] + g[2][1] + 2 * g[0][2] + 2 * g[2][0];

    // total gametes
    int n = g[0][0] + g[0][1] + g[0][2] + g[1][0] + g[1][1] + g[1][2] + g[2][0] + g[2][1] + g[2][2];
    n *= 2;

    // double heterozygote
    int h = 2 * g[1][1];

    double r = 0.25;
    double conv = 0.0;

    for (int i = 0; i < maxit; ++i) {
        double p = r * r / (1 - 2 * r + 2 * r * r);
        double r_new = (m + p * h) / n;

        conv = std::fabs(r - r_new);
        if (conv < tol)
            return r_new;

        r = r_new;
    }

    eprint("WARNING: maximum number of iterations exceeded in est_rec_frac_f2: "
           "maxit = %d, tol = %g, conv = %g\n", maxit, tol, conv);

    return r;
}

// BC: P1(AABB) * P2(aabb) -> F1(AaBb) -> * P1(AABB)
//
// Gamete   AB        Ab    aB    ab
// AB       (1-r)/2   r/2   r/2   (1-r)/2
//
// Genotype   BB        Bb
// AA         (1-r)/2   r/2
// Aa         r/2       (1-r)/2
//
//         BB      Bb
//   AA   AB//AB   AB//Ab      g[0][0]  g[0][1]
//   Aa   AB//aB   AB//ab      g[1][0]  g[1][1]

double est_rec_frac_bc(int g[2][2])
{
    int m = g[0][1] + g[1][0];
    int n = g[0][0] + g[0][1] + g[1][0] + g[1][1];
    return (double) m / n;
}

// RIL: P1(AABB) * P2(aabb) -> F1(AaBb) -> repeated selfing for at least seven generations
//
//         BB      bb
//   AA   AB//AB   Ab//Ab      g[0][0]  g[0][1]
//   aa   aB//aB   ab//ab      g[1][0]  g[1][1]

double est_rec_frac_riself(int g[2][2])
{
    double R = est_rec_frac_bc(g);
    return R / 2 / (1 - R);
}
