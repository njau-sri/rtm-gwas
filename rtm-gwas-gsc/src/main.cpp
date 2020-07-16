#include <exception>
#include "print.h"
#include "rtmgwasgsc.h"

int main(int argc, char *argv[])
{
    try {
        return rtm_gwas_gsc(argc, argv);
    }
    catch (const std::exception &e) {
        eprint("FATAL: exception caught: %s\n", e.what());
        return 1;
    }
    catch (...) {
        eprint("FATAL: unknown exception caught\n");
        return 1;
    }
}
