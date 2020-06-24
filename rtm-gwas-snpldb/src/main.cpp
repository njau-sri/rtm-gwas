#include <iostream>
#include <exception>

#include "rtmgwassnpldb.h"

int main(int argc, char *argv[])
{
    try {
        return rtm_gwas_snpldb(argc, argv);
    }
    catch (const std::exception &e) {
        std::cerr << "FATAL: exception caught: " << e.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cerr << "FATAL: unknown exception caught\n";
        return 1;
    }
}
