#include "rtmgwasgconv.h"

#include "types.h"
#include "print.h"
#include "cmdline.h"
#include "stringutil.h"
#include "vectorutil.h"
#include "vcf.h"
#include "ped.h"
#include "hmp.h"
#include "geno.h"

#ifndef RTM_GWAS_VERSION
#define RTM_GWAS_VERSION "unknown"
#endif // RTM_GWAS_VERSION

namespace {

    struct Parameter
    {
        std::string vcf;
        std::string ped;
        std::string hmp;
        std::string geno;
        std::string out;
        bool sort = false;
    } par;

    void show_help()
    {
        const char *help =
            "Usage: rtm-gwas-gconv [options...]\n"
            "Options:\n"
            "  --vcf <>   VCF genotype file (.vcf)\n"
            "  --ped <>   PLINK PED file (.ped)\n"
            "  --hmp <>   HapMap genotype file\n"
            "  --geno <>  General genotype file\n"
            "  --out <>   output file with format suffix (.vcf/.ped/.hmp/.geno)\n"
            "  --sort     sort locus in ascending chromosome/position order\n"
            "\n";
        eprint(help);
    }

    int parse_cmdline(int argc, char *argv[])
    {
        CmdLine cmd;

        cmd.add("--vcf", "");
        cmd.add("--ped", "");
        cmd.add("--hmp", "");
        cmd.add("--geno", "");
        cmd.add("--out", "");
        cmd.add("--sort");

        if (cmd.parse(argc, argv) != 0)
            return 1;

        par.vcf = cmd.get("--vcf");
        par.ped = cmd.get("--ped");
        par.hmp = cmd.get("--hmp");
        par.geno = cmd.get("--geno");
        par.out = cmd.get("--out");
        par.sort = cmd.has("--sort");

        return 0;
    }

    int load_genotype(Genotype &gt)
    {
        if (!par.vcf.empty())
            return read_vcf(par.vcf, gt);

        if (!par.ped.empty()) {
            std::string prefix = par.ped;
            if (ends_with(prefix, ".ped") || ends_with(prefix, ".map"))
                prefix.erase(length(prefix) - 4, 4);
            return read_ped(prefix, gt);
        }

        if (!par.hmp.empty())
            return read_hmp(par.hmp, gt);

        if (!par.geno.empty())
            return read_geno(par.geno, gt);

        eprint("ERROR: no input genotype file specified\n");

        return 1;
    }

    void sort_chrpos(Genotype &gt)
    {
        auto chr = stable_unique(gt.chr);

        isize_t m = length(gt.loc);
        std::vector<isize_t> z;
        z.reserve(m);

        for (auto &e : chr) {
            std::vector<isize_t> idx;
            for (isize_t i = 0; i < m; ++i) {
                if (gt.chr[i] == e)
                    idx.push_back(i);
            }

            auto pos = subset(gt.pos, idx);
            if (!std::is_sorted(pos.begin(), pos.end())) {
                auto ord = order_asc(pos);
                subset(idx, ord).swap(idx);
            }

            z.insert(z.end(), idx.begin(), idx.end());
        }

        subset(gt.loc, z).swap(gt.loc);
        subset(gt.chr, z).swap(gt.chr);
        subset(gt.pos, z).swap(gt.pos);
        subset(gt.dat, z).swap(gt.dat);
        subset(gt.allele, z).swap(gt.allele);
    }

} // namespace

int rtm_gwas_gconv(int argc, char *argv[])
{
    eprint("RTM-GWAS %s GCONV (Built on %s %s)\n", RTM_GWAS_VERSION, __DATE__, __TIME__);

    if (argc < 2) {
        show_help();
        return 1;
    }

    if (parse_cmdline(argc, argv) != 0)
        return 1;

    Genotype gt;

    eprint("INFO: reading genotype file...\n");

    if (load_genotype(gt) != 0)
        return 1;

    eprint("INFO: %td individuals, %td loci\n", length(gt.ind), length(gt.loc));

    if (par.sort)
        sort_chrpos(gt);

    int info = 0;

    if (ends_with(par.out, ".vcf"))
        info = write_vcf(gt, par.out);
    else if (ends_with(par.out, ".ped")) {
        std::string prefix = par.out.substr(0, length(par.out) - 4);
        info = write_ped(gt, prefix);
    }
    else if (ends_with(par.out, ".hmp"))
        info = write_hmp(gt, par.out);
    else if (ends_with(par.out, ".geno"))
        info = write_geno(gt, par.out);
    else {
        eprint("ERROR: unrecognized output format: %s\n", par.out);
        return 1;
    }

    if (info != 0)
        return 1;

    eprint("INFO: GCONV has finished successfully\n");

    return 0;
}
