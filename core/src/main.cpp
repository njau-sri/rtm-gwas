#include <memory>
#include <iostream>
#include <exception>
#include "appldb.h"
#include "appgsc.h"
#include "appassoc.h"
#include "appdata.h"
#include "appsum.h"

#ifdef _OPENMP
extern "C" {
    int omp_get_num_procs();
    int omp_get_max_threads();
    void omp_set_num_threads(int);
}
#endif

int main(int argc, char *argv[])
{
    std::cerr << "RTM-GWAS v1.0 (Built on " __DATE__ " at " __TIME__ ")\n"
                 "  https://github.com/njau-sri/rtm-gwas\n";

    if (argc < 2) {
        std::cerr << "Usage: rtm-gwas [command]\n"
                     "Command:\n"
                     "  ldb     SNPLDB marker construction\n"
                     "  gsc     Genetic similarity coefficient\n"
                     "  assoc   Association analysis\n"
                     "  data    Genotype management\n"
                     "  sum     Genotype summary\n";
        return 1;
    }

#ifdef _OPENMP
    int num_procs = omp_get_num_procs();
    if (num_procs < omp_get_max_threads() * 2)
        omp_set_num_threads(num_procs < 4 ? 1 : num_procs / 2);
#endif

    try {
        string command = argv[1];
        if (command == "ldb") {
            return std::make_shared<AppLDB>()->run(argc-1, argv+1);
        }
        else if (command == "gsc") {
            return std::make_shared<AppGSC>()->run(argc-1, argv+1);
        }
        else if (command == "assoc") {
            return std::make_shared<AppAssoc>()->run(argc-1, argv+1);
        }
        else if (command == "data") {
            return std::make_shared<AppData>()->run(argc-1, argv+1);
        }
        else if (command == "sum") {
            return std::make_shared<AppSum>()->run(argc-1, argv+1);
        }
        else {
            std::cerr << "ERROR: unrecognized command: " << command << "\n";
            return 1;
        }
    }
    catch (const std::exception &e) {
        std::cerr << "FATAL: exception caught: " << e.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cerr << "FATAL: unknown exception caught\n";
        return 1;
    }

    return 0;
}
