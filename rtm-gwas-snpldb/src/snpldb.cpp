#include <numeric>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "cmdline.h"
#include "vcf.h"
#include "split.h"
#include "util.h"
#include "gabriel.h"
#include "haplotype.h"


#ifndef RTM_GWAS_VERSION
#define RTM_GWAS_VERSION  "2019.4.dev"
#endif


using std::size_t;


namespace {


struct Parameter
{
    std::string vcf;
    std::string block;
    std::string gene;
    std::string out;
    std::string fam;
    double maf = 0.01;
    double inform = 0.95;
    int maxlen = 100000;
    int llim = 70;
    int ulim = 98;
    int recomb = 90;
    bool openmp = false;
} par ;


int read_block(const std::string &filename, std::vector<std::string> &chrom,
               std::vector<int> &start, std::vector<int> &stop)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    size_t ln = 0;
    for (std::string line; std::getline(ifs, line); ) {
        ++ln;

        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if ( vs.empty() )
            continue;

        if (vs.size() < 3) {
            std::cerr << "ERROR: expected at least 3 columns at line " << ln << "\n";
            return 1;
        }

        auto pos1 = std::stoi(vs[1]);
        auto pos2 = std::stoi(vs[2]);

        if (pos2 <= pos1) {
            std::cerr << "ERROR: invalid block position at line " << ln << ": " << vs[1] << " " << vs[2] << "\n";
            return 1;
        }

        chrom.push_back(vs[0]);
        start.push_back(pos1);
        stop.push_back(pos2);
    }

    auto n = chrom.size();
    bool require_sort = false;

    std::vector<size_t> ord;
    ord.reserve(n);

    for (auto &e: stable_unique(chrom)) {
        std::vector<size_t> idx;
        for (size_t i = 0; i < n; ++i) {
            if (chrom[i] == e)
                idx.push_back(i);
        }
        auto pos = subset(start, idx);
        if ( ! std::is_sorted(pos.begin(), pos.end()) ) {
            require_sort = true;
            subset(idx,order(pos)).swap(idx);
        }

        size_t j = 0;
        int prev_stop = -1;
        for (auto i : idx) {
            if (start[i] <= prev_stop) {
                std::cerr << "ERROR: overlapping block is not allowed:\n"
                          << "   " << chrom[j] << " " << start[j] << " " << stop[j] << "\n"
                          << "   " << chrom[i] << " " << start[i] << " " << stop[i] << "\n";
                return 1;
            }
            j = i;
            prev_stop = stop[i];
        }

        ord.insert(ord.end(), idx.begin(), idx.end());
    }

    if ( require_sort ) {
        subset(chrom, ord).swap(chrom);
        subset(start, ord).swap(start);
        subset(stop, ord).swap(stop);
    }

    return 0;
}

int read_gene(const std::string &filename, std::vector<std::string> &gene, std::vector<std::string> &chrom,
              std::vector<int> &start, std::vector<int> &stop)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    size_t ln = 0;
    for (std::string line; std::getline(ifs, line); ) {
        ++ln;

        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if ( vs.empty() )
            continue;

        if (vs.size() < 4) {
            std::cerr << "ERROR: expected at least 4 columns at line: " << ln << "\n";
            return 1;
        }

        auto pos1 = std::stoi(vs[2]);
        auto pos2 = std::stoi(vs[3]);

        if (pos2 <= pos1) {
            std::cerr << "ERROR: invalid gene position at line " << ln << ": " << vs[2] << " " << vs[3] << "\n";
            return 1;
        }

        gene.push_back(vs[0]);
        chrom.push_back(vs[1]);
        start.push_back(pos1);
        stop.push_back(pos2);
    }

    auto n = chrom.size();
    bool require_sort = false;

    std::vector<size_t> ord;
    ord.reserve(n);

    for (auto &e: stable_unique(chrom)) {
        std::vector<size_t> idx;
        for (size_t i = 0; i < n; ++i) {
            if (chrom[i] == e)
                idx.push_back(i);
        }
        auto pos = subset(start, idx);
        if ( ! std::is_sorted(pos.begin(), pos.end()) ) {
            require_sort = true;
            subset(idx,order(pos)).swap(idx);
        }

        size_t j = 0;
        int prev_stop = -1;
        for (auto i : idx) {
            if (start[i] <= prev_stop) {
                std::cerr << "ERROR: overlapping gene is not allowed:\n"
                          << "   " << gene[j] << " " << chrom[j] << " " << start[j] << " " << stop[j] << "\n"
                          << "   " << gene[i] << " " << chrom[i] << " " << start[i] << " " << stop[i] << "\n";
                return 1;
            }
            j = i;
            prev_stop = stop[i];
        }

        ord.insert(ord.end(), idx.begin(), idx.end());
    }

    if ( require_sort ) {
        subset(gene, ord).swap(gene);
        subset(chrom, ord).swap(chrom);
        subset(start, ord).swap(start);
        subset(stop, ord).swap(stop);
    }

    return 0;
}

void calc_intergenic_region(const std::vector<std::string> &gene_chrom,
                            const std::vector<int> &gene_start,
                            const std::vector<int> &gene_stop,
                            std::vector<std::string> &igr_chrom,
                            std::vector<int> &igr_start,
                            std::vector<int> &igr_stop)
{
    auto n = gene_chrom.size();
    for (auto &e : stable_unique(gene_chrom)) {
        int start = 1, stop = 1;
        for (size_t i = 0; i < n; ++i) {
            if (gene_chrom[i] != e)
                continue;
            if (gene_start[i] > start) {
                stop = gene_start[i] - 1;
                igr_chrom.push_back(e);
                igr_start.push_back(start);
                igr_stop.push_back(stop);
            }
            start = gene_stop[i] + 1;
        }
        stop = std::numeric_limits<int>::max();
        igr_chrom.push_back(e);
        igr_start.push_back(start);
        igr_stop.push_back(stop);
    }
}

std::vector< std::vector<size_t> > index_snp(const Genotype &gt, const std::vector<std::string> &chr)
{
    auto k = chr.size();

    std::map<std::string,size_t> chridx;
    for (size_t i = 0; i < k; ++i)
        chridx[chr[i]] = i;

    auto m = gt.loc.size();

    std::vector< std::vector<size_t> > sidx(k);
    for (size_t i = 0; i < m; ++i) {
        auto j = chridx[gt.chr[i]];
        sidx[j].push_back(i);
    }

    for (auto &v : sidx) {
        auto pos = subset(gt.pos, v);
        if ( ! std::is_sorted(pos.begin(), pos.end()) ) {
            auto ord = order(pos);
            subset(v,ord).swap(v);
        }
    }

    return sidx;
}

std::vector< std::vector<size_t> > index_family(const std::vector<size_t> &fam)
{
    std::vector< std::vector<size_t> > idx;

    size_t start = 0;
    for (auto n : fam) {
        std::vector<size_t> v(n);
        std::iota(v.begin(), v.end(), start);
        idx.push_back(v);
        start += n;
    }

    return idx;
}

void recode_block_allele(std::vector< std::vector<std::string> > &allele)
{
    for (auto &v : allele) {
        bool require = false;
        for (auto &e : v) {
            if (e.size() > 1) {
                require = true;
                break;
            }
        }
        if ( require ) {
            auto n = v.size();
            for (size_t i = 0; i < n; ++i)
                v[i] = std::to_string(i);
        }
    }
}

void write_block_allele(const Genotype &gt, const std::vector< std::vector<std::string> > &allele)
{
    std::ofstream ofs(par.out + ".allele");
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << par.out << ".allele\n";
        return;
    }

    auto m = gt.loc.size();

    ofs << "Locus\tCode\tAllele\n";

    for (size_t j = 0; j < m; ++j) {
        auto n = gt.allele[j].size();
        for (size_t k = 0; k < n; ++k)
            ofs << gt.loc[j] << "\t" << allele[j][k] << "\t" << gt.allele[j][k] << "\n";
    }
}

// AA,Aa,aa,NN -> 0,1,2,3
void recode_012(const Genotype &gt, const std::vector<size_t> &sidx,
                std::vector<int> &pos, std::vector< std::vector<char> > &dat)
{
    auto n = gt.ind.size();
    auto m = sidx.size();

    pos.resize(m);

    dat.resize(m);
    for (auto &v : dat)
        v.resize(n);

    size_t k = 0;

    if (gt.ploidy == 2) {
        for (auto j : sidx) {
            for (size_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][2*i];
                auto b = gt.dat[j][2*i+1];
                char c = 3;  // missing genotype
                if (a == 1) {
                    if (b == 1)
                        c = 0;
                    else if (b == 2)
                        c = 1;
                }
                else if (a == 2) {
                    if (b == 1)
                        c = 1;
                    else if (b == 2)
                        c = 2;
                }
                dat[k][i] = c;
            }
            pos[k++] = gt.pos[j];
        }
    }
    else {
        for (auto j : sidx) {
            for (size_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i];
                char c = 3;  // missing genotype
                if (a == 1)
                    c = 0;
                else if (a == 2)
                    c = 2;
                dat[k][i] = c;
            }
            pos[k++] = gt.pos[j];
        }
    }
}

int define_block(const Genotype &gt, std::vector<std::string> &chrom, std::vector<int> &start, std::vector<int> &stop)
{
    BlockGabriel gab;
    gab.maxlen = par.maxlen;
    gab.llim = par.llim;
    gab.ulim = par.ulim;
    gab.recomb = par.recomb;
    gab.frac = par.inform;

    auto chrid = stable_unique(gt.chr);
    auto sidx = index_snp(gt, chrid);

    auto nchr = chrid.size();

    for (size_t i = 0; i < nchr; ++i) {
        std::cerr << "INFO: finding block on " << chrid[i] << "\n";

        std::vector<int> pos;
        std::vector< std::vector<char> > dat;
        recode_012(gt, sidx[i], pos, dat);

        std::vector< std::pair<int,int> > bpos;

        int ret = par.openmp ? find_block_omp(gab, pos, dat, bpos) :
                               find_block(gab, pos, dat, bpos);
        if (ret != 0)
            return 1;

        chrom.insert(chrom.end(), bpos.size(), chrid[i]);
        for (auto &e : bpos) {
            start.push_back(e.first);
            stop.push_back(e.second);
        }
    }

    return 0;
}

int define_gblock(const Genotype &gt, std::vector<std::string> &name, std::vector<std::string> &chrom,
                  std::vector<int> &start, std::vector<int> &stop)
{
    auto ng = name.size();
    auto m = gt.loc.size();

    std::vector<char> ignore(m, 0);
    for (size_t i = 0; i < m; ++i) {
        auto chr = gt.chr[i];
        auto pos = gt.pos[i];
        for (size_t j = 0; j < ng; ++j) {
            if (chrom[j] == chr && start[j] <= pos && stop[j] >= pos) {
                ignore[i] = 1;
                break;
            }
        }
    }

    std::vector<std::string> igr_chrom;
    std::vector<int> igr_start, igr_stop;
    calc_intergenic_region(chrom, start, stop, igr_chrom, igr_start, igr_stop);

    BlockGabriel gab;
    gab.maxlen = par.maxlen;
    gab.llim = par.llim;
    gab.ulim = par.ulim;
    gab.recomb = par.recomb;
    gab.frac = par.inform;

    auto chrid = stable_unique(gt.chr);
    auto sidx = index_snp(gt, chrid);

    for (size_t i = 0; i < chrid.size(); ++i) {
        std::cerr << "INFO: finding block on " << chrid[i] << "\n";

        std::string prefix = "LDB_" + chrid[i];

        std::vector<int> igr_start_chr, igr_stop_chr;
        for (size_t k = 0; k < igr_chrom.size(); ++k) {
            if (igr_chrom[k] == chrid[i]) {
                igr_start_chr.push_back(igr_start[k]);
                igr_stop_chr.push_back(igr_stop[k]);
            }
        }

        std::vector<int> pos;
        std::vector< std::vector<char> > dat;

        for (size_t k = 0; k < igr_start_chr.size(); ++k) {
            std::vector<size_t> idx;
            for (auto j : sidx[i]) {
                if (ignore[j])
                    continue;
                if (gt.pos[j] >= igr_start_chr[k] && gt.pos[j] <= igr_stop_chr[k]) {
                    idx.push_back(j);
                    ignore[j] = 1;
                }
            }
            if (idx.size() > 1) {
                recode_012(gt, idx, pos, dat);

                std::vector< std::pair<int,int> > bpos;
                int ret = par.openmp ? find_block_omp(gab, pos, dat, bpos) :
                                       find_block(gab, pos, dat, bpos);
                if (ret != 0)
                    return 1;

                chrom.insert(chrom.end(), bpos.size(), chrid[i]);
                for (auto &e : bpos) {
                    start.push_back(e.first);
                    stop.push_back(e.second);
                    name.push_back(prefix + "_" + std::to_string(e.first) + "_" + std::to_string(e.second));
                }
            }
        }
    }

    return 0;
}

int separate_genotype(size_t n, Genotype &gt, Genotype &pgt)
{
    bool diploid = gt.ploidy == 2;

    for (size_t i = 0; i < n; ++i)
        pgt.ind.push_back( gt.ind[i] );

    for (size_t i = 0; i < n; ++i)
        gt.ind.erase( gt.ind.begin() );

    for (auto &v : gt.dat) {
        std::vector<allele_t> w(n);

        if ( diploid ) {
            // NOTE: presuming that parents are homozygous
            for (size_t i = 0; i < n; ++i) {
                if (v[i*2] != v[i*2+1]) {
                    std::cerr << "ERROR: parent genotype is not homozygous\n";
                    return 1;
                }
                w[i] = v[i*2];
            }
        }
        else {
            for (size_t i = 0; i < n; ++i)
                w[i] = v[i];
        }

        pgt.dat.push_back(w);

        for (size_t i = 0; i < n; ++i)
            v.erase( v.begin() );
    }

    return 0;
}

// RIL layout: p1 p2 ind1 ind2 ...
// NAM layout: cp p1 p2 p3 f1-ind1 f1-ind2 ... f2-ind1 f2-ind2 ... f3-ind1 f3-ind2 ...
int rtm_gwas_snpldb_fam()
{
    size_t tfam = 0;
    std::vector<size_t> fam;
    for (auto &e : split(par.fam," \t,")) {
        auto a = std::stoi(e);
        if (a < 1) {
            std::cerr << "ERROR: invalid family size: " << e << "\n";
            return 1;
        }
        fam.push_back( static_cast<size_t>(a) );
        tfam += fam.back();
    }

    if ( fam.empty() ) {
        std::cerr << "ERROR: invalid arguments for --fam: " << par.fam << "\n";
        return 1;
    }

    std::cerr << "INFO: " << fam.size() << " " << (fam.size() <= 1 ? "family was" : "families were")
              << " observed in the population\n";

    Genotype gt;

    std::cerr << "INFO: reading genotype file...\n";
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    if (gt.ploidy != 1 && gt.ploidy != 2) {
        std::cerr << "ERROR: unsupported genotype ploidy: " << gt.ploidy << "\n";
        return 1;
    }

    size_t np = fam.size() + 1;

    if (gt.ind.size() != tfam + np) {
        std::cerr << "ERROR: inconsistent number of individuals: " << gt.ind.size() << " " << tfam + np << "\n";
        return 1;
    }

    Genotype pgt;
    if (separate_genotype(np, gt, pgt) != 0)
        return 1;

    // Define block
    int ret = 0;
    std::vector<std::string> blk_name, blk_chr;
    std::vector<int> blk_start, blk_stop;

    if ( ! par.gene.empty() ) {
        std::cerr << "INFO: reading gene list file...\n";
        ret = read_gene(par.gene, blk_name, blk_chr, blk_start, blk_stop);
        if (ret != 0)
            return 1;
        std::cerr << "INFO: " << blk_name.size() << " genes\n";

        if ( blk_name.empty() ) {
            std::cerr << "ERROR: no valid gene could be found\n";
            return 1;
        }

        ret = define_gblock(gt, blk_name, blk_chr, blk_start, blk_stop);
        if (ret != 0)
            return 1;
    }
    else if ( ! par.block.empty() ) {
        std::cerr << "INFO: reading predefined block file...\n";
        ret = read_block(par.block, blk_chr, blk_start, blk_stop);
        if (ret != 0)
            return 1;
        std::cerr << "INFO: " << blk_start.size() << " blocks\n";
    }
    else {
        ret = define_block(gt, blk_chr, blk_start, blk_stop);
        if (ret != 0)
            return 1;
    }

    if ( blk_start.empty() ) {
        std::cerr << "ERROR: no blocks are found\n";
        return 1;
    }

    auto m = gt.loc.size();
    auto nb = blk_start.size();

    Genotype ibgt;
    Genotype pbgt;
    std::vector<char> inblock(m,0);

    std::vector<int> blk_length(nb,0);
    std::vector<size_t> blk_size(nb,0);
    std::vector<size_t> blk_rec(nb,0);

    auto chrid = stable_unique(gt.chr);
    auto sidx = index_snp(gt, chrid);
    auto fidx = index_family(fam);

    auto nchr = chrid.size();

    for (size_t i = 0; i < nchr; ++i) {
        Genotype iht;
        Genotype pht;

        for (size_t k = 0; k < nb; ++k) {
            if (blk_chr[k] != chrid[i])
                continue;
            std::vector<size_t> jidx;
            for (auto j : sidx[i]) {
                if (gt.pos[j] < blk_start[k] || gt.pos[j] > blk_stop[k])
                    continue;
                inblock[j] = true;
                jidx.push_back(j);
            }
            blk_length[k] = blk_stop[k] - blk_start[k] + 1;
            blk_size[k] = jidx.size();
            if ( jidx.empty() ) {
                std::cerr << "WARNING: no SNPs were found in block: " << chrid[i] << " "
                          << blk_start[k] << " " << blk_stop[k] << "\n";
                continue;
            }

            blk_rec[k] = infer_haplotype_ril(gt, pgt, jidx, fidx, iht, pht);

            iht.chr.push_back(chrid[i]);

            if ( blk_name.empty() ) {  // block-info-based locus name
                std::string loc = "LDB_";
                loc += blk_chr[k];
                loc += "_";
                loc += std::to_string(blk_start[k]);
                loc += "_";
                loc += std::to_string(blk_stop[k]);
                iht.loc.push_back(loc);
            }
            else
                iht.loc.push_back(blk_name[k]);
        }

        for (auto j : sidx[i]) {
            if ( inblock[j] ) {
                auto k = index(iht.pos, gt.pos[j]);
                if (k != iht.pos.size()) {
                    ibgt.loc.push_back(iht.loc[k]);
                    ibgt.chr.push_back(iht.chr[k]);
                    ibgt.pos.push_back(iht.pos[k]);
                    ibgt.dat.push_back(iht.dat[k]);
                    ibgt.allele.push_back(iht.allele[k]);
                    pbgt.dat.push_back(pht.dat[k]);
                }
            }
            else {
                // unified locus name
                std::string loc = "LDB_" + gt.chr[j] + "_" + std::to_string(gt.pos[j]);
                ibgt.loc.push_back(loc);
                ibgt.chr.push_back(gt.chr[j]);
                ibgt.pos.push_back(gt.pos[j]);
                ibgt.dat.push_back(gt.dat[j]);
                ibgt.allele.push_back(gt.allele[j]);
                pbgt.dat.push_back(pgt.dat[j]);
            }
        }
    }

    std::ofstream ofs(par.out + ".block");

    if ( ! ofs )
        std::cerr << "ERROR: can't open file for writing: " << par.out << ".block\n";
    else {
        if ( ! blk_name.empty() )
            ofs << "Block\t";
        ofs << "Chromosome\tStart\tStop\tLength\tSNPs\tRecombination\n";
        for (size_t i = 0; i < nb; ++i) {
            if ( ! blk_name.empty() )
                ofs << blk_name[i] << "\t";
            ofs << blk_chr[i] << "\t" << blk_start[i] << "\t" << blk_stop[i] << "\t"
                << blk_length[i] << "\t" << blk_size[i] << "\t" << blk_rec[i] << "\n";
        }
    }

    ibgt.ind = gt.ind;
    ibgt.ploidy = gt.ploidy;

    auto allele = ibgt.allele;
    recode_block_allele(allele);
    write_block_allele(ibgt, allele);

    ibgt.allele.swap(allele);

    if (write_vcf(ibgt, par.out + ".vcf") != 0)
        return 3;

    ibgt.ploidy = 1;
    ibgt.ind = pgt.ind;
    ibgt.dat = pbgt.dat;

    if (write_vcf(ibgt, par.out + ".parent.vcf") != 0)
        return 4;

    return 0;
}


} // namespace


int rtm_gwas_snpldb(int argc, char *argv[])
{
    std::cerr << "RTM-GWAS " RTM_GWAS_VERSION " SNPLDB (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd;

    cmd.add("--vcf", "VCF file", "");
    cmd.add("--block", "predefined block file", "");
    cmd.add("--gene", "gene coordinate file", "");
    cmd.add("--out", "output file prefix", "snpldb.out");
    cmd.add("--maf", "minimum minor haplotype frequency", "0.01");
    cmd.add("--maxlen", "maximum length of blocks", "100000");

    cmd.add("--llim", "lower limit CI for strong LD", "70");
    cmd.add("--ulim", "upper limit CI for string LD", "98");
    cmd.add("--recomb", "upper limit CI for strong recombination", "90");
    cmd.add("--inform", "minimum fraction of informative strong LD", "0.95");

    cmd.add("--fam", "sample sizes in RIL/NAM, n1,n2,n3,...", "");

    cmd.add("--openmp", "enable OpenMP multithreading");

    cmd.parse(argc, argv);

    if (argc < 2) {
        cmd.show();
        return 1;
    }

    par.vcf = cmd.get("--vcf");
    par.block = cmd.get("--block");
    par.gene = cmd.get("--gene");
    par.out = cmd.get("--out");

    par.maf = std::stod(cmd.get("--maf"));
    par.maxlen = std::stoi(cmd.get("--maxlen"));

    par.llim = std::stoi(cmd.get("--llim"));
    par.ulim = std::stoi(cmd.get("--ulim"));
    par.recomb = std::stoi(cmd.get("--recomb"));
    par.inform = std::stod(cmd.get("--inform"));

    par.fam = cmd.get("--fam");

    par.openmp = cmd.has("--openmp");

    if ( ! par.fam.empty() )
        return rtm_gwas_snpldb_fam();

    // Load genotype

    Genotype gt;
    std::cerr << "INFO: reading genotype file...\n";
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    if (gt.ploidy != 1 && gt.ploidy != 2) {
        std::cerr << "ERROR: unsupported genotype ploidy: " << gt.ploidy << "\n";
        return 1;
    }

    // Define block

    int ret = 0;
    std::vector<std::string> blk_name, blk_chr;
    std::vector<int> blk_start, blk_stop;

    if ( ! par.gene.empty() ) {
        std::cerr << "INFO: reading gene list file...\n";
        ret = read_gene(par.gene, blk_name, blk_chr, blk_start, blk_stop);
        if (ret != 0)
            return 1;
        std::cerr << "INFO: " << blk_name.size() << " genes\n";

        if ( blk_name.empty() ) {
            std::cerr << "ERROR: no valid gene could be found\n";
            return 1;
        }

        ret = define_gblock(gt, blk_name, blk_chr, blk_start, blk_stop);
        if (ret != 0)
            return 1;
    }
    else if ( ! par.block.empty() ) {
        std::cerr << "INFO: reading predefined block file...\n";
        ret = read_block(par.block, blk_chr, blk_start, blk_stop);
        if (ret != 0)
            return 1;
        std::cerr << "INFO: " << blk_start.size() << " blocks\n";
    }
    else {
        ret = define_block(gt, blk_chr, blk_start, blk_stop);
        if (ret != 0)
            return 1;
    }

    if ( blk_start.empty() ) {
        std::cerr << "WARNING: no blocks could be found, no results generated\n";
        return 0;
    }

    // Group SNPs within block

    auto m = gt.loc.size();
    auto nb = blk_start.size();

    Genotype bgt;
    std::vector<char> inblock(m,0);

    std::vector<int> blk_length(nb,0);
    std::vector<size_t> blk_size(nb,0);
    std::vector<size_t> blk_rec(nb,0);

    auto chrid = stable_unique(gt.chr);
    auto sidx = index_snp(gt, chrid);

    auto nchr = chrid.size();

    for (size_t i = 0; i < nchr; ++i) {
        Genotype ht;

        for (size_t k = 0; k < nb; ++k) {
            if (blk_chr[k] != chrid[i])
                continue;

            std::vector<size_t> jidx;
            for (auto j : sidx[i]) {
                if (gt.pos[j] < blk_start[k] || gt.pos[j] > blk_stop[k])
                    continue;
                inblock[j] = true;
                jidx.push_back(j);
            }

            blk_length[k] = blk_stop[k] - blk_start[k] + 1;
            blk_size[k] = jidx.size();

            if ( jidx.empty() ) {
                std::cerr << "WARNING: no SNPs were found in block: " << chrid[i] << " "
                          << blk_start[k] << " " << blk_stop[k] << "\n";
                continue;
            }

            blk_rec[k] = infer_haplotype(par.maf, gt, jidx, ht);

            ht.chr.push_back(chrid[i]);

            if ( blk_name.empty() ) {  // block-info-based locus name
                std::string loc = "LDB_";
                loc += blk_chr[k];
                loc += "_";
                loc += std::to_string(blk_start[k]);
                loc += "_";
                loc += std::to_string(blk_stop[k]);
                ht.loc.push_back(loc);
            }
            else
                ht.loc.push_back(blk_name[k]);
        }

        for (auto j : sidx[i]) {
            if ( inblock[j] ) {
                auto k = index(ht.pos, gt.pos[j]);
                if (k != ht.pos.size()) {
                    bgt.loc.push_back(ht.loc[k]);
                    bgt.chr.push_back(ht.chr[k]);
                    bgt.pos.push_back(ht.pos[k]);
                    bgt.dat.push_back(ht.dat[k]);
                    bgt.allele.push_back(ht.allele[k]);
                }
            }
            else {
                // unified locus name
                std::string loc = "LDB_" + gt.chr[j] + "_" + std::to_string(gt.pos[j]);
                bgt.loc.push_back(loc);
                bgt.chr.push_back(gt.chr[j]);
                bgt.pos.push_back(gt.pos[j]);
                bgt.dat.push_back(gt.dat[j]);
                bgt.allele.push_back(gt.allele[j]);
            }
        }
    }

    std::ofstream ofs(par.out + ".block");

    if ( ! ofs )
        std::cerr << "ERROR: can't open file for writing: " << par.out << ".block\n";
    else {
        if ( ! blk_name.empty() )
            ofs << "Block\t";
        ofs << "Chromosome\tStart\tStop\tLength\tSNPs\tRecombination\n";
        for (size_t i = 0; i < nb; ++i) {
            if (blk_size[i] == 0)
                continue;
            if ( ! blk_name.empty() )
                ofs << blk_name[i] << "\t";
            ofs << blk_chr[i] << "\t" << blk_start[i] << "\t" << blk_stop[i] << "\t"
                << blk_length[i] << "\t" << blk_size[i] << "\t" << blk_rec[i] << "\n";
        }
    }

    bgt.ind = gt.ind;
    bgt.ploidy = gt.ploidy;

    auto allele = bgt.allele;
    recode_block_allele(allele);
    write_block_allele(bgt, allele);

    bgt.allele.swap(allele);

    if (write_vcf(bgt, par.out + ".vcf") != 0)
        return 3;

    return 0;
}
