#include <iostream>
#include <exception>

#include "rtmgwasassoc.h"

int main(int argc, char *argv[])
{
    try {
        return rtm_gwas_assoc(argc, argv);
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
