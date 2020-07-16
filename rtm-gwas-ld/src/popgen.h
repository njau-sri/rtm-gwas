#ifndef POPGEN_H
#define POPGEN_H

// proportion of haplotype AB//ab in double-heterozygote
double est_het(int maxit, double tol, int AB, int Ab, int aB, int ab, int AaBb);

// haplotype frequency, p = { p_AB, p_Ab, p_aB, p_ab }
void est_hap_freq(int maxit, double tol, int g[3][3], double p[4]);

// linkage disequilibrium
void calc_dprime_rsq(double pa, double pb, double pab, double &dprime, double &rsq);

// recombination fraction in F2 population
double est_rec_frac_f2(int maxit, double tol, int g[3][3]);

// recombination fraction in BC population
double est_rec_frac_bc(int g[2][2]);

// recombination fraction in RIL (selfing) population
double est_rec_frac_riself(int g[2][2]);

#endif // POPGEN_H
