#include <memory>
#include <fstream>
#include <iostream>
#include <iterator>
#include "appdata.h"
#include "cmdline.h"
#include "util.h"
#include "number.h"
#include "vcf.h"
#include "plink.h"
#include "hapmap.h"
#include "data.h"

int AppData::run(int argc, char *argv[])
{
    auto cmd = std::make_shared<CmdLine>("data [options]");

    cmd->add("--vcf", "VCF file", "");
    cmd->add("--ped", "PLINK PED/MAP file prefix", "");
    cmd->add("--hmp", "HapMap file", "");
    cmd->add("--geno", "genotype file", "");
    cmd->add("--out", "output file", "appdata.out");
    cmd->add("--format", "output file format", "vcf");
    cmd->add("--loc", "include locus list file", "");
    cmd->add("--ind", "include individual list file", "");
    cmd->add("--maf", "minimum minor allele frequency", "0");
    cmd->add("--het", "maximum heterozygosity", "1");
    cmd->add("--mind", "maximum per-individual missing", "1");
    cmd->add("--mloc", "maximum per-locus missing", "1");
    cmd->add("--min-alleles", "minimum number of alleles", "0");
    cmd->add("--max-alleles", "maximum number of alleles", "255");

    cmd->add("--diploid", "force output diploid genotype");

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
    m_par.format = cmd->get("--format");
    m_par.loc = cmd->get("--loc");
    m_par.ind = cmd->get("--ind");

    m_par.maf = number<double>(cmd->get("--maf"));
    m_par.het = number<double>(cmd->get("--het"));
    m_par.mloc = number<double>(cmd->get("--mloc"));
    m_par.mind = number<double>(cmd->get("--mind"));

    m_par.min_alleles = number<int>(cmd->get("--min-alleles"));
    m_par.max_alleles = number<int>(cmd->get("--max-alleles"));

    m_par.diploid = cmd->has("--diploid");

    int info = perform();

    return info;
}

int AppData::perform()
{
    load_genotype();

    apply_filter();

    save_genotype();

    return 0;
}

void AppData::load_genotype()
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

void AppData::apply_filter()
{
    m_keep_ind.assign(m_gt.ind.size(), true);

    m_keep_loc.assign(m_gt.loc.size(), true);

    filter_ind();

    filter_mis_ind();

    filter_loc();

    filter_allele_number();

    filter_mis_loc();

    filter_maf();

    filter_het();

    if (std::count(m_keep_ind.begin(), m_keep_ind.end(), false) != 0) {
        if ( m_gt.dat.empty() || m_gt.dat[0].size() == m_gt.ind.size() ) {
            for (auto &e : m_gt.dat)
                subset(e, m_keep_ind).swap(e);
        }
        else {
            auto n = m_keep_ind.size();
            vector<bool> keep;
            for (size_t i = 0; i < n; ++i) {
                keep.push_back(m_keep_ind[i]);
                keep.push_back(m_keep_ind[i]);
            }
            for (auto &e : m_gt.dat)
                subset(e, keep).swap(e);
        }
        subset(m_gt.ind, m_keep_ind).swap(m_gt.ind);
    }

    if (std::count(m_keep_loc.begin(), m_keep_loc.end(), false) != 0) {
        subset(m_gt.loc, m_keep_loc).swap(m_gt.loc);
        subset(m_gt.chr, m_keep_loc).swap(m_gt.chr);
        subset(m_gt.pos, m_keep_loc).swap(m_gt.pos);
        subset(m_gt.dist, m_keep_loc).swap(m_gt.dist);
        subset(m_gt.dat, m_keep_loc).swap(m_gt.dat);
        subset(m_gt.allele, m_keep_loc).swap(m_gt.allele);
    }
}

void AppData::save_genotype() const
{
    std::cerr << "INFO: writing genotype to file...\n";

    int info = 1;

    if (m_par.format == "vcf") {
        info = write_vcf(m_gt, m_par.out + ".vcf", m_par.diploid);
    }
    else if (m_par.format == "ped") {
        if ( is_compatible_plink(m_gt) )
            info = write_plink(m_gt, m_par.out);
        else
            std::cerr << "ERROR: can't convert to PLINK PED/MAP format\n";
    }
    else if (m_par.format == "hmp") {
        if ( is_compatible_hapmap(m_gt) )
            info = write_hapmap(m_gt, m_par.out + ".hmp.txt");
        else
            std::cerr << "ERROR: can't convert to HapMap format\n";
    }
    else if (m_par.format == "geno") {
        info = write_genotype(m_gt, m_par.out + ".geno.txt");
    }
    else {
        if ( ! m_par.format.empty() )
            std::cerr << "WARN: invalid output format: " << m_par.format << "\n";
        info = write_vcf(m_gt, m_par.out, m_par.diploid);
    }

    if (info == 0)
        std::cerr << "INFO: genotype data were successfully saved in " << m_par.format << " format\n";
}

void AppData::filter_ind()
{
    if ( m_par.ind.empty() )
        return;

    std::ifstream ifs(m_par.ind);

    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file: " << m_par.ind << "\n";
        return;
    }

    vector<string> vs;
    std::copy(std::istream_iterator<string>(ifs), std::istream_iterator<string>(), std::back_inserter(vs));

    std::sort(vs.begin(), vs.end());

    if ( vs.empty() )
        return;

    auto n = m_gt.ind.size();
    for (size_t i = 0; i < n; ++i) {
        if ( m_keep_ind[i] && ! std::binary_search(vs.begin(), vs.end(), m_gt.ind[i]) )
            m_keep_ind[i] = false;
    }

    auto nv = std::count(m_keep_ind.begin(), m_keep_ind.end(), true);
    std::cerr << "INFO: after filtering individual list, there are " << nv << " individuals\n";
}

void AppData::filter_mis_ind()
{
    if (m_par.mind >= 1.0)
        return;

    auto m = m_gt.loc.size(), n = m_gt.ind.size();
    bool haploid = m_gt.ploidy == 1;

    size_t mt = std::count(m_keep_loc.begin(), m_keep_loc.end(), true);
    if ( ! haploid )
        mt *= 2;

    for (size_t i = 0; i < n; ++i) {
        if ( ! m_keep_ind[i] )
            continue;

        size_t mv = 0;
        if ( haploid ) {
            for (size_t j = 0; j < m; ++j) {
                if ( m_keep_loc[j] && m_gt.dat[j][i] )
                    ++mv;
            }
        }
        else {
            for (size_t j = 0; j < m; ++j) {
                if ( ! m_keep_loc[j] )
                    continue;
                if ( m_gt.dat[j][i*2] )
                    ++mv;
                if ( m_gt.dat[j][i*2+1] )
                    ++mv;
            }
        }

        if ((double) mv / mt > m_par.mind)
            m_keep_ind[i] = false;
    }

    auto nv = std::count(m_keep_ind.begin(), m_keep_ind.end(), true);
    std::cerr << "INFO: after filtering MIND<=" << m_par.mind << ", there are " << nv << " individuals\n";
}

void AppData::filter_loc()
{
    if ( m_par.loc.empty() )
        return;

    std::ifstream ifs(m_par.loc);

    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file: " << m_par.loc << "\n";
        return;
    }

    vector<string> vs;
    std::copy(std::istream_iterator<string>(ifs), std::istream_iterator<string>(), std::back_inserter(vs));

    std::sort(vs.begin(), vs.end());

    if ( vs.empty() )
        return;

    auto m = m_gt.loc.size();
    for (size_t j = 0; j < m; ++j) {
        if ( m_keep_loc[j] && ! std::binary_search(vs.begin(), vs.end(), m_gt.loc[j]) )
            m_keep_loc[j] = false;
    }

    auto mv = std::count(m_keep_loc.begin(), m_keep_loc.end(), true);
    std::cerr << "INFO: after filtering locus list, there are " << mv << " loci\n";
}

void AppData::filter_allele_number()
{
    if (m_par.min_alleles <= 0 && m_par.max_alleles >= 255)
        return;

    auto m = m_gt.loc.size();
    for (size_t j = 0; j < m; ++j) {
        if ( ! m_keep_loc[j] )
            continue;
        int k = m_gt.allele[j].size();
        if (k < m_par.min_alleles || k > m_par.max_alleles)
            m_keep_loc[j] = false;
    }

    auto mv = std::count(m_keep_loc.begin(), m_keep_loc.end(), true);
    std::cerr << "INFO: after filtering allele number, there are " << mv << " loci\n";
}

void AppData::filter_mis_loc()
{
    if (m_par.mloc >= 1.0)
        return;

    auto m = m_gt.loc.size(), n = m_gt.ind.size();
    bool haploid = m_gt.ploidy == 1;

    for (size_t j = 0; j < m; ++j) {
        if ( ! m_keep_loc[j] )
            continue;
        size_t ns = 0, nt = 0;
        if ( haploid ) {
            for (size_t i = 0; i < n; ++i) {
                if ( ! m_keep_ind[i] )
                    continue;
                if (m_gt.dat[j][i] == 0)
                    ++ns;
                ++nt;
            }
        }
        else {
            for (size_t i = 0; i < n; ++i) {
                if ( ! m_keep_ind[i] )
                    continue;
                if (m_gt.dat[j][i*2] == 0 || m_gt.dat[j][i*2+1] == 0)
                    ++ns;
                ++nt;
            }
        }

        if (ns == nt || (double) ns / nt > m_par.mloc)
            m_keep_loc[j] = false;
    }

    auto mv = std::count(m_keep_loc.begin(), m_keep_loc.end(), true);
    std::cerr << "INFO: after filtering MLOC<=" << m_par.mloc << ", there are " << mv << " loci\n";
}

void AppData::filter_maf()
{
    if (m_par.maf <= 0.0)
        return;

    auto m = m_gt.loc.size(), n = m_gt.ind.size();
    bool haploid = m_gt.ploidy == 1;

    for (size_t j = 0; j < m; ++j) {
        if ( ! m_keep_loc[j] )
            continue;

        std::map<allele_t,int> ac;

        if ( haploid ) {
            for (size_t i = 0; i < n; ++i) {
                if ( ! m_keep_ind[i] )
                    continue;
                auto a = m_gt.dat[j][i];
                if (ac.find(a) == ac.end())
                    ac[a] = 0;
                ++ac[a];
            }
        }
        else {
            for (size_t i = 0; i < n; ++i) {
                if ( ! m_keep_ind[i] )
                    continue;
                auto a = m_gt.dat[j][i*2];
                auto b = m_gt.dat[j][i*2+1];
                if (ac.find(a) == ac.end())
                    ac[a] = 0;
                if (ac.find(b) == ac.end())
                    ac[b] = 0;
                ++ac[a];
                ++ac[b];
            }
        }

        int ac_min = n*2, ac_tot = 0;
        for (auto e : ac) {
            if ( ! e.first )
                continue;
            if (e.second < ac_min)
                ac_min = e.second;
            ac_tot += e.second;
        }

        if (ac_tot == 0 || ac_min == ac_tot || (double) ac_min / ac_tot < m_par.maf)
            m_keep_loc[j] = false;
    }

    auto mv = std::count(m_keep_loc.begin(), m_keep_loc.end(), true);
    std::cerr << "INFO: after filtering MAF>=" << m_par.maf << ", there are " << mv << " loci\n";
}

void AppData::filter_het()
{
    if (m_par.het >= 1.0)
        return;

    if (m_gt.dat.empty() || m_gt.dat[0].size() == m_gt.ind.size())
        return;

    auto m = m_gt.loc.size(), n = m_gt.ind.size();
    for (size_t j = 0; j < m; ++j) {
        if ( ! m_keep_loc[j] )
            continue;
        size_t ns = 0, nt = 0;
        for (size_t i = 0; i < n; ++i) {
            if ( ! m_keep_ind[i] )
                continue;
            auto a = m_gt.dat[j][i*2];
            auto b = m_gt.dat[j][i*2+1];
            if (a != 0 && b != 0) {
                ++nt;
                if (a != b)
                    ++ns;
            }
        }

        if (nt != 0 && (double) ns / nt > m_par.het)
            m_keep_loc[j] = false;
    }

    auto mv = std::count(m_keep_loc.begin(), m_keep_loc.end(), true);
    std::cerr << "INFO: after filtering HET<=" << m_par.het << ", there are " << mv << " loci\n";
}
