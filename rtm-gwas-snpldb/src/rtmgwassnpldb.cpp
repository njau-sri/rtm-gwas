#include "rtmgwassnpldb.h"

#include "types.h"
#include "print.h"
#include "openmp.h"
#include "cmdline.h"
#include "vcf.h"
#include "cfile.h"
#include "stringutil.h"
#include "vectorutil.h"
#include "hapblock.h"
#include "haplotype.h"

#ifndef RTM_GWAS_VERSION
#define RTM_GWAS_VERSION "unknown"
#endif // RTM_GWAS_VERSION

namespace {

    struct Parameter
    {
        std::string vcf;
        std::string block;
        std::string gene;
        std::string out;
        std::string nam;
        double maf = 0.01;
        double identity = 0.0;
        double inform = 0.95;
        int maxlen = 100000;
        int llim = 70;
        int ulim = 98;
        int recomb = 90;
        int thread = 0;
        bool ril = false;
    } par;

    void show_help()
    {
        const char *help =
            "Usage: rtm-gwas-snpldb [options...]\n"
            "Options:\n"
            "  --vcf <>            SNP genotype data file in VCF format\n"
            "  --out <snpldb.out>  output file prefix\n"
            "  --maxlen <100000>   maximum length of haplotype block (bp)\n"
            "  --maf <0.01>        minimum minor haplotype frequency [0,1]\n"
            "  --identity <0>      minimum identity for haplotype replacement [0,1]\n"
            "  --gene <>           gene coordinate file\n"
            "  --block <>          reference haplotype block file\n"
            "  --llim <70>         lower limit CI for strong LD [0,100]\n"
            "  --ulim <98>         upper limit CI for string LD [0,100]\n"
            "  --recomb <90>       upper limit CI for strong recombination [0,100]\n"
            "  --inform <0.95>     minimum fraction of informative strong LD [0,1]\n"
            "  --nam <>            NAM population mode, specify the number of lines\n"
            "                      in each family, n1,n2,n3,...\n"
            "  --ril               RIL population mode, the first two individuals\n"
            "                      must be parental lines\n"
            "  --thread <0>        set the number of threads\n"
            "\n";
        eprint(help);
    }

    void parse_cmdline(int argc, char *argv[])
    {
        CmdLine cmd;

        cmd.add("--vcf", "");
        cmd.add("--out", "snpldb.out");
        cmd.add("--maxlen", "100000");
        cmd.add("--maf", "0.01");
        cmd.add("--identity", "0");
        cmd.add("--gene", "");
        cmd.add("--block", "");
        cmd.add("--llim", "70");
        cmd.add("--ulim", "98");
        cmd.add("--recomb", "90");
        cmd.add("--inform",  "0.95");
        cmd.add("--nam", "");
        cmd.add("--ril");
        cmd.add("--thread", "0");

        cmd.parse(argc, argv);

        par.vcf = cmd.get("--vcf");
        par.block = cmd.get("--block");
        par.gene = cmd.get("--gene");
        par.out = cmd.get("--out");

        par.maf = std::stod(cmd.get("--maf"));
        par.identity = std::stod(cmd.get("--identity"));
        par.maxlen = std::stoi(cmd.get("--maxlen"));

        par.llim = std::stoi(cmd.get("--llim"));
        par.ulim = std::stoi(cmd.get("--ulim"));
        par.recomb = std::stoi(cmd.get("--recomb"));
        par.inform = std::stod(cmd.get("--inform"));

        par.nam = cmd.get("--nam");
        par.ril = cmd.has("--ril");

        par.thread = std::stoi(cmd.get("--thread"));
    }

    int read_block(const std::string &filename, std::vector<std::string> &chrom,
                   std::vector<int> &start, std::vector<int> &stop)
    {
        CFileLineReader file(filename);

        if (!file) {
            eprint("ERROR: can't open file for reading: %s\n", filename);
            return 1;
        }

        isize_t ln = 0;
        for (std::string line; file.read(line); ) {
            ++ln;

            auto vs = split(line, " \t");
            if (vs.empty())
                continue;

            if (length(vs) < 3) {
                eprint("ERROR: expected at least 3 columns at line %td\n", ln);
                return 1;
            }

            int pos1 = std::stoi(vs[1]);
            int pos2 = std::stoi(vs[2]);

            if (pos2 <= pos1) {
                eprint("ERROR: invalid block position at line %td: %s %s\n", ln, vs[1], vs[2]);
                return 1;
            }

            chrom.push_back(vs[0]);
            start.push_back(pos1);
            stop.push_back(pos2);
        }

        isize_t n = length(chrom);
        bool require_sort = false;

        std::vector<isize_t> ord;
        ord.reserve(n);

        for (auto &e: stable_unique(chrom)) {
            std::vector<isize_t> idx;
            for (isize_t i = 0; i < n; ++i) {
                if (chrom[i] == e)
                    idx.push_back(i);
            }
            auto pos = subset(start, idx);
            if (!std::is_sorted(pos.begin(), pos.end())) {
                require_sort = true;
                subset(idx, order_asc(pos)).swap(idx);
            }

            isize_t j = 0;
            int prev_stop = -1;
            for (auto i : idx) {
                if (start[i] <= prev_stop) {
                    eprint("ERROR: overlapping block is not allowed:\n");
                    eprint("       %s %d %d\n", chrom[j], start[j], stop[j]);
                    eprint("       %s %d %d\n", chrom[i], start[i], stop[i]);
                    return 1;
                }
                j = i;
                prev_stop = stop[i];
            }

            ord.insert(ord.end(), idx.begin(), idx.end());
        }

        if (require_sort) {
            subset(chrom, ord).swap(chrom);
            subset(start, ord).swap(start);
            subset(stop, ord).swap(stop);
        }

        return 0;
    }

    int read_gene(const std::string &filename, std::vector<std::string> &gene, std::vector<std::string> &chrom,
                  std::vector<int> &start, std::vector<int> &stop)
    {
        CFileLineReader file(filename);

        if (!file) {
            eprint("ERROR: can't open file for reading: %s\n", filename);
            return 1;
        }

        isize_t ln = 0;
        for (std::string line; file.read(line); ) {
            ++ln;

            auto vs = split(line, " \t");
            if (vs.empty())
                continue;

            if (length(vs) < 4) {
                eprint("ERROR: expected at least 4 columns at line %td\n", ln);
                return 1;
            }

            int pos1 = std::stoi(vs[2]);
            int pos2 = std::stoi(vs[3]);

            if (pos2 <= pos1) {
                eprint("ERROR: invalid gene position at line %td: %s %s\n", ln, vs[2], vs[3]);
                return 1;
            }

            gene.push_back(vs[0]);
            chrom.push_back(vs[1]);
            start.push_back(pos1);
            stop.push_back(pos2);
        }

        isize_t n = length(chrom);
        bool require_sort = false;

        std::vector<isize_t> ord;
        ord.reserve(n);

        for (auto &e: stable_unique(chrom)) {
            std::vector<isize_t> idx;
            for (isize_t i = 0; i < n; ++i) {
                if (chrom[i] == e)
                    idx.push_back(i);
            }
            auto pos = subset(start, idx);
            if (!std::is_sorted(pos.begin(), pos.end())) {
                require_sort = true;
                subset(idx, order_asc(pos)).swap(idx);
            }

            isize_t j = 0;
            int prev_stop = -1;
            for (auto i : idx) {
                if (start[i] <= prev_stop) {
                    eprint("ERROR: overlapped gene is not allowed:\n");
                    eprint("       %s %s %d %d\n", gene[j], chrom[j], start[j], stop[j]);
                    eprint("       %s %s %d %d\n", gene[i], chrom[i], start[i], stop[i]);
                    return 1;
                }
                j = i;
                prev_stop = stop[i];
            }

            ord.insert(ord.end(), idx.begin(), idx.end());
        }

        if (require_sort) {
            subset(gene, ord).swap(gene);
            subset(chrom, ord).swap(chrom);
            subset(start, ord).swap(start);
            subset(stop, ord).swap(stop);
        }

        return 0;
    }

    std::string rename_block(const std::string &chr, int start, int stop)
    {
        std::string s = "BLK";
        s.append("_").append(chr);
        s.append("_").append(std::to_string(start));
        s.append("_").append(std::to_string(stop));
        return s;
    }

    std::vector< std::vector<isize_t> > index_loc(const Genotype &gt, const std::vector<std::string> &chr)
    {
        isize_t n = length(chr);

        std::map<std::string, isize_t> chridx;
        for (isize_t i = 0; i < n; ++i)
            chridx[chr[i]] = i;

        isize_t m = length(gt.loc);

        std::vector< std::vector<isize_t> > idx(n);
        for (isize_t i = 0; i < m; ++i) {
            isize_t j = chridx[gt.chr[i]];
            idx[j].push_back(i);
        }

        for (auto &v : idx) {
            auto pos = subset(gt.pos, v);
            if (!std::is_sorted(pos.begin(), pos.end())) {
                auto ord = order_asc(pos);
                subset(v, ord).swap(v);
            }
        }

        return idx;
    }

    // AA:0, Aa:1, aa:2, missing:-1
    void recode_snp_012(const Genotype &gt, const std::vector<isize_t> &sidx,
                        std::vector<int> &pos, std::vector< std::vector<int8_t> > &dat)
    {
        isize_t n = length(gt.ind);
        isize_t m = length(sidx);

        pos.resize(m);
        dat.assign(m, std::vector<int8_t>(n, -1));

        if (gt.ploidy == 2) {
            for (isize_t j = 0; j < m; ++j) {
                isize_t jj = sidx[j];
                pos[j] = gt.pos[jj];
                for (isize_t i = 0; i < n; ++i) {
                    auto a = gt.dat[jj][i*2];
                    auto b = gt.dat[jj][i*2+1];
                    if (a == 1) {
                        if (b == 1)
                            dat[j][i] = 0;
                        else if (b == 2)
                            dat[j][i] = 1;
                    }
                    else if (a == 2) {
                        if (b == 1)
                            dat[j][i] = 1;
                        else if (b == 2)
                            dat[j][i] = 2;
                    }
                }
            }
        }
        else {
            for (isize_t j = 0; j < m; ++j) {
                isize_t jj = sidx[j];
                pos[j] = gt.pos[jj];
                for (isize_t i = 0; i < n; ++i) {
                    auto a = gt.dat[jj][i];
                    if (a == 1)
                        dat[j][i] = 0;
                    else if (a == 2)
                        dat[j][i] = 2;
                }
            }
        }
    }

    int define_block(const Genotype &gt, std::vector<std::string> &chrom, std::vector<int> &start, std::vector<int> &stop)
    {
        HapBlockGabriel hb;
        hb.par.thread = par.thread;
        hb.par.maxlen = par.maxlen;
        hb.par.llim = par.llim;
        hb.par.ulim = par.ulim;
        hb.par.recomb = par.recomb;
        hb.par.frac = par.inform;

        auto chrid = stable_unique(gt.chr);
        auto sidx = index_loc(gt, chrid);

        isize_t nchr = length(chrid);

        for (isize_t i = 0; i < nchr; ++i) {
            if (length(sidx[i]) < 2)
                continue;

            eprint("INFO: finding block on chromosome %s\n", chrid[i]);

            std::vector<int> pos;
            std::vector< std::vector<int8_t> > dat;
            recode_snp_012(gt, sidx[i], pos, dat);

            if (find_hapblock_gabriel(pos, dat, &hb) != 0)
                return 1;

            for (auto &e : hb.block) {
                chrom.push_back(chrid[i]);
                start.push_back(e.first);
                stop.push_back(e.second);
            }
        }

        return 0;
    }

    int define_igr_block_chr(const std::vector<int> &pos, const std::vector< std::vector<int8_t> > &dat,
                             const std::vector<int> &gstart, const std::vector<int> &gstop,
                             std::vector< std::pair<int,int> > &block)
    {
        HapBlockGabriel hb;
        hb.par.thread = par.thread;
        hb.par.maxlen = par.maxlen;
        hb.par.llim = par.llim;
        hb.par.ulim = par.ulim;
        hb.par.recomb = par.recomb;
        hb.par.frac = par.inform;

        if (gstart.empty()) {
            if (find_hapblock_gabriel(pos, dat, &hb) != 0)
                return 1;
            block.insert(block.end(), hb.block.begin(), hb.block.end());
            return 0;
        }

        isize_t m = length(pos);
        isize_t g = length(gstart);

        std::vector<char> ignore(m, 0);
        for (isize_t j = 0; j < m; ++j) {
            for (isize_t i = 0; i < g; ++i) {
                if (pos[j] >= gstart[i] && pos[j] <= gstop[i]) {
                    ignore[j] = 1;
                    break;
                }
            }
        }

        int istart = 1, istop = 1;
        for (isize_t i = 0; i < g; ++i) {
            if (gstart[i] > istart) {
                istop = gstart[i] - 1;

                std::vector<isize_t> idx;
                for (isize_t j = 0; j < m; ++j) {
                    if (ignore[j])
                        continue;
                    if (pos[j] >= istart && pos[j] <= istop) {
                        idx.push_back(j);
                        ignore[j] = 1;
                    }
                }

                if (find_hapblock_gabriel(subset(pos,idx), subset(dat,idx), &hb))
                    return 1;

                block.insert(block.end(), hb.block.begin(), hb.block.end());
            }

            istart = gstop[i] + 1;
        }

        istop = std::numeric_limits<int>::max();

        std::vector<isize_t> idx;
        for (isize_t j = 0; j < m; ++j) {
            if (ignore[j])
                continue;
            if (pos[j] >= istart && pos[j] <= istop) {
                idx.push_back(j);
                ignore[j] = 1;
            }
        }

        if (find_hapblock_gabriel(subset(pos,idx), subset(dat,idx), &hb))
            return 1;

        block.insert(block.end(), hb.block.begin(), hb.block.end());

        return 0;
    }

    int define_igr_block(const Genotype &gt, std::vector<std::string> &name, std::vector<std::string> &chrom,
                         std::vector<int> &start, std::vector<int> &stop)
    {
        isize_t g = length(name);

        auto chrid = stable_unique(gt.chr);
        auto idx = index_loc(gt, chrid);

        isize_t nchr = length(chrid);

        for (isize_t i = 0; i < nchr; ++i) {
            if (length(idx[i]) < 2)
                continue;

            eprint("INFO: finding block on chromosome %s\n", chrid[i]);

            std::vector<int> gstart, gstop;
            for (isize_t j = 0; j < g; ++j) {
                if (chrom[j] == chrid[i]) {
                    gstart.push_back(start[j]);
                    gstop.push_back(stop[j]);
                }
            }

            std::vector<int> pos;
            std::vector< std::vector<int8_t> > dat;
            recode_snp_012(gt, idx[i], pos, dat);

            std::vector< std::pair<int, int> > block;
            int info = define_igr_block_chr(pos, dat, gstart, gstop, block);

            if (info != 0)
                return 1;

            for (auto &e : block) {
                name.push_back( rename_block(chrid[i], e.first, e.second) );
                chrom.push_back(chrid[i]);
                start.push_back(e.first);
                stop.push_back(e.second);
            }
        }

        return 0;
    }

    //void recode_haplotype(std::vector< std::vector<std::string> > &haps)
    //{
    //    for (auto &v : haps) {
    //        bool require = false;
    //        for (auto &e : v) {
    //            if (e.size() > 1) {
    //                require = true;
    //                break;
    //            }
    //        }
    //        if (require) {
    //            isize_t n = length(v);
    //            for (isize_t i = 0; i < n; ++i)
    //                v[i] = std::to_string(i);
    //        }
    //    }
    //}

    void recode_haplotype(std::vector< std::vector<std::string> > &haps)
    {
        for (auto &v : haps) {
            isize_t n = length(v);
            for (isize_t i = 0; i < n; ++i)
                v[i] = std::to_string(i);
        }
    }

    void write_block_allele(const Genotype &gt, const std::vector< std::vector<std::string> > &allele)
    {
        CFile file(par.out + ".allele", "w");

        if (!file) {
            eprint("ERROR: can't open file for writing: %s.allele\n", par.out);
            return;
        }

        isize_t m = length(gt.loc);

        fprint(file, "Locus\tCode\tAllele\n");

        for (isize_t j = 0; j < m; ++j) {
            isize_t n = length(gt.allele[j]);
            for (isize_t k = 0; k < n; ++k)
                fprint(file, "%s\t%s\t%s\n", gt.loc[j], allele[j][k], gt.allele[j][k]);
        }
    }

    int rtm_gwas_snpldb_ril(Genotype &gt)
    {
        eprint("ERROR: function not implemented\n");
        return 1;
    }

    int rtm_gwas_snpldb_nam(Genotype &gt)
    {
        eprint("ERROR: function not implemented\n");
        return 1;
    }

} // namespace

int rtm_gwas_snpldb(int argc, char *argv[])
{
    eprint("RTM-GWAS %s SNPLDB (Built on %s %s)\n", RTM_GWAS_VERSION, __DATE__, __TIME__);

    if (argc < 2) {
        show_help();
        return 1;
    }

    parse_cmdline(argc, argv);

    if (par.thread > 0)
        call_omp_set_num_threads(par.thread);

    // load genotype

    Genotype gt;
    eprint("INFO: reading genotype file...\n");
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    eprint("INFO: %td individuals, %td loci\n", length(gt.ind), length(gt.loc));

    if (gt.ploidy != 1 && gt.ploidy != 2) {
        eprint("ERROR: unsupported ploidy of genotype data: %d\n", gt.ploidy);
        return 1;
    }

    if (par.ril)
        return rtm_gwas_snpldb_ril(gt);

    if (!par.nam.empty())
        return rtm_gwas_snpldb_nam(gt);

    // define block

    std::vector<std::string> bname, bchr;
    std::vector<int> bstart, bstop;

    if (!par.gene.empty()) {
        eprint("INFO: reading gene list file...\n");

        int info = read_gene(par.gene, bname, bchr, bstart, bstop);

        if (info != 0)
            return 1;

        eprint("INFO: %td genes\n", length(bname));

        if (bname.empty()) {
            eprint("ERROR: no valid gene could be found\n");
            return 1;
        }

        info = define_igr_block(gt, bname, bchr, bstart, bstop);

        if (info != 0)
            return 1;
    }
    else if (!par.block.empty()) {
        eprint("INFO: reading predefined block file...\n");

        int info = read_block(par.block, bchr, bstart, bstop);

        if (info != 0)
            return 1;

        eprint("INFO: %td blocks\n", length(bstart));
    }
    else {
        int info = define_block(gt, bchr, bstart, bstop);

        if (info != 0)
            return 1;
    }

    if (bstart.empty()) {
        eprint("WARNING: no blocks could be found, no results generated\n");
        return 0;
    }

    isize_t nb = length(bchr);

    if (bname.empty()) {
        for (isize_t i = 0; i < nb; ++i)
            bname.push_back( rename_block(bchr[i], bstart[i], bstop[i]) );
    }

    // form haplotype

    Haplotype ht;
    ht.par.maf = par.maf;
    ht.par.identity = par.identity;

    auto chrid = stable_unique(gt.chr);
    auto sidx = index_loc(gt, chrid);

    isize_t m = length(gt.loc);
    isize_t nchr = length(chrid);

    std::vector<char> inblock(m, 0);
    std::vector<int> bsize(nb, 0);
    std::vector<int> bfilter(nb, 0);

    Genotype snpldb;

    CFile file(par.out + ".block", "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s.block\n", par.out);
        return 1;
    }

    fprint(file, "Block\tChromosome\tStart\tStop\tLength\tSNPs\tFilter\n");

    isize_t nempty = 0;
    isize_t nsingle = 0;

    for (isize_t i = 0; i < nchr; ++i) {
        Genotype bgt;
        for (isize_t k = 0; k < nb; ++k) {
            if (bchr[k] != chrid[i])
                continue;

            std::vector<isize_t> idx;
            for (auto j : sidx[i]) {
                if (gt.pos[j] >= bstart[k] && gt.pos[j] <= bstop[k]) {
                    inblock[j] = 1;
                    idx.push_back(j);
                }
            }

            if (idx.empty()) {
                ++nempty;
                continue;
            }

            if (length(idx) == 1) {
                inblock[idx[0]] = 0;
                gt.loc[idx[0]].append("_").append(bname[k]);
                ++nsingle;
                continue;
            }

            form_haplotype(gt, idx, &ht);

            int amax = * std::max_element(ht.dat.begin(), ht.dat.end());
            if (amax > std::numeric_limits<allele_t>::max()) {
                eprint("WARNING: *drop* block exceeding maximum number (%d) of alleles: %s %d %d %d\n",
                       std::numeric_limits<allele_t>::max(), chrid[i], bstart[k], bstop[k], amax);
                for (auto j : idx) inblock[j] = 0;
                // TODO: repartition block
                continue;
            }

            fprint(file, "%s\t%s\t%d\t%d\t%d\t%d\t%d\n", bname[k], ht.chr, bstart[k], bstop[k],
                   bstop[k] - bstart[k] + 1, ht.size, ht.filter);

            bgt.loc.push_back(bname[k]);
            bgt.chr.push_back(ht.chr);
            bgt.pos.push_back(bstart[k]);

            bgt.dat.emplace_back(ht.dat.size());
            std::transform(ht.dat.begin(), ht.dat.end(), bgt.dat.back().begin(),
                           [](int a) { return static_cast<allele_t>(a); });

            bgt.allele.push_back(ht.hap);
        }

        for (auto j : sidx[i]) {
            if (inblock[j]) {
                isize_t k = index(bgt.pos, gt.pos[j]);
                if (k != -1) {
                    snpldb.loc.push_back(bgt.loc[k]);
                    snpldb.chr.push_back(bgt.chr[k]);
                    snpldb.pos.push_back(bgt.pos[k]);
                    snpldb.dat.push_back(bgt.dat[k]);
                    snpldb.allele.push_back(bgt.allele[k]);

                }
            }
            else {
                snpldb.loc.push_back(gt.loc[j]);
                snpldb.chr.push_back(gt.chr[j]);
                snpldb.pos.push_back(gt.pos[j]);
                snpldb.dat.push_back(gt.dat[j]);
                snpldb.allele.push_back(gt.allele[j]);
            }
        }
    }

    if (nempty > 0)
        eprint("WARNING: skipped %td blocks without SNP\n", nempty);

    if (nsingle > 0)
        eprint("WARNING: skipped %td blocks with only one SNP\n", nsingle);

    snpldb.ind = gt.ind;
    snpldb.ploidy = gt.ploidy;

    auto allele = snpldb.allele;
    recode_haplotype(allele);
    write_block_allele(snpldb, allele);

    snpldb.allele.swap(allele);

    if (write_vcf(snpldb, par.out + ".vcf") != 0)
        return 1;

    eprint("INFO: SNPLDB has finished successfully\n");

    return 0;
}
