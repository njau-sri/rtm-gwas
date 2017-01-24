#include <map>
#include <memory>
#include <fstream>
#include <iostream>
#include "appsum.h"
#include "cmdline.h"
#include "vcfio.h"
#include "plinkio.h"
#include "hapmapio.h"
#include "rtmio.h"

namespace {

std::map<allele_t,int> count_allele(const vector<allele_t> &v)
{
    std::map<allele_t,int> ac;

    for (auto e : v) {
        if (ac.find(e) == ac.end())
            ac[e] = 0;
        ++ac[e];
    }

    return ac;
}

} // namespace

int AppSum::run(int argc, char *argv[])
{
    auto cmd = std::make_shared<CmdLine>("sum [options]");

    cmd->add("--vcf", "VCF file", "");
    cmd->add("--ped", "PLINK PED/MAP file prefix", "");
    cmd->add("--hmp", "HapMap file", "");
    cmd->add("--geno", "genotype file", "");
    cmd->add("--out", "output file", "appsum.out");
    cmd->add("--allele", "output allele statistics");
    cmd->add("--indiv", "output individual statistics");

    if (argc < 2) {
        cmd->help();
        return 1;
    }

    cmd->parse(argc, argv);

    m_par.vcf = cmd->get("--vcf");
    m_par.ped = cmd->get("--ped");
    m_par.hmp = cmd->get("--hmp");
    m_par.geno = cmd->get("--geno");
    m_par.out = cmd->get("--out");

    m_par.allele = cmd->has("--allele");
    m_par.indiv = cmd->has("--indiv");

    int info = perform();

    return info;
}

int AppSum::perform()
{
    load_genotype();

    if (m_par.allele)
        calc_save_stats_allele();
    else if (m_par.indiv)
        calc_save_stats_indiv();
    else
        calc_save_stats_locus();

    return 0;
}

void AppSum::load_genotype()
{
    if ( m_par.vcf.empty() && m_par.ped.empty() && m_par.hmp.empty() && m_par.geno.empty() )
        return;

    std::cerr << "INFO: reading genotype file...\n";

    int info = 0;

    if ( ! m_par.vcf.empty() )
        info = read_vcf(m_par.vcf, m_gt);
    else if ( ! m_par.ped.empty() )
        info = read_plink(m_par.ped, m_gt);
    else if ( ! m_par.hmp.empty() )
        info = read_hapmap(m_par.hmp, m_gt);
    else if ( ! m_par.geno.empty() )
        info = read_genotype(m_par.geno, m_gt);

    if (info != 0) {
        m_gt.loc.clear();
        m_gt.ind.clear();
        m_gt.dat.clear();
    }

    std::cerr << "INFO: " << m_gt.ind.size() << " individuals and " << m_gt.loc.size() << " loci were observed\n";
}

void AppSum::calc_save_stats_locus() const
{
    std::ofstream ofs(m_par.out + ".locus.txt");

    ofs << "Locus\tAlleleNumber\tMinAlleleCount\tMinAlleleFreq\tMaxAlleleCount\tMaxAlleleFreq\t"""
           "TotalAlleleCount\tTotalAlleleFreq\n";

    auto m = m_gt.loc.size();

    for (size_t j = 0; j < m; ++j) {
        auto n = m_gt.dat[j].size();
        auto ac = count_allele(m_gt.dat[j]);

        int ac_min = n, ac_max = 0, ac_tot = 0;
        for (auto e : ac) {
            if ( ! e.first )
                continue;
            if (e.second < ac_min)
                ac_min = e.second;
            if (e.second > ac_max)
                ac_max = e.second;
            ac_tot += e.second;
        }

        ofs << m_gt.loc[j] << "\t" << m_gt.allele[j].size() << "\t"
            << ac_min << "\t" << (double) ac_min / ac_tot << "\t"
            << ac_max << "\t" << (double) ac_max / ac_tot << "\t"
            << ac_tot << "\t" << (double) ac_tot / n << "\n";
    }
}

void AppSum::calc_save_stats_indiv() const
{
    std::ofstream ofs(m_par.out + ".indiv.txt");

    ofs << "Indiv\tTotalAlleleCount\tTotalAlleleFreq\n";

    auto m = m_gt.loc.size(), n = m_gt.ind.size();
    bool haploid = m_gt.ploidy == 1;
    auto mt = haploid ? m : m*2;

    for (size_t i = 0; i < n; ++i) {
        size_t mv = 0;

        if ( haploid ) {
            for (size_t j = 0; j < m; ++j) {
                if ( m_gt.dat[j][i] )
                    ++mv;
            }
        }
        else {
            for (size_t j = 0; j < m; ++j) {
                if ( m_gt.dat[j][i*2] )
                    ++mv;
                if ( m_gt.dat[j][i*2+1] )
                    ++mv;
            }
        }

        ofs << m_gt.ind[i] << "\t" << mv << "\t" << (double) mv / mt << "\n";
    }
}

void AppSum::calc_save_stats_allele() const
{
    std::ofstream ofs(m_par.out + ".allele.txt");

    ofs << "Locus\tAllelel\tCount\tFreq\n";

    auto m = m_gt.loc.size();

    for (size_t j = 0; j < m; ++j) {
        auto ac = count_allele(m_gt.dat[j]);

        size_t tc = 0;
        for (auto e : ac) {
            if ( e.first )
                tc += e.second;
        }

        for (auto e : ac) {
            if ( ! e.first )
                continue;
            ofs << m_gt.loc[j] << "\t" << m_gt.allele[j][e.first-1] << "\t"
                << e.second << "\t" << (double) e.second / tc << "\n";
        }
    }
}
