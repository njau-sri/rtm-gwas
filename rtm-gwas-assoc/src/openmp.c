#include <stdio.h>
#include "openmp.h"

#ifdef _OPENMP

int omp_get_max_threads();

int omp_get_num_procs();

void omp_set_num_threads(int num_threads);

int call_omp_get_max_threads()
{
    return omp_get_max_threads();
}

int call_omp_get_num_procs()
{
    return omp_get_num_procs();
}

void call_omp_set_num_threads(int num_threads)
{
    omp_set_num_threads(num_threads);
}

#else

int call_omp_get_max_threads()
{
    fputs("WARNING: OpenMP not enabled: call_omp_get_max_threads\n", stderr);
    return 1;
}

int call_omp_get_num_procs()
{
    fputs("WARNING: OpenMP not enabled: call_omp_get_num_procs\n", stderr);
    return 1;
}

void call_omp_set_num_threads(int num_threads)
{
    fputs("WARNING: OpenMP not enabled: call_omp_set_num_threads\n", stderr);
    (void) num_threads;
}

#endif // _OPENMP
