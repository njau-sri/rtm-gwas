#ifndef OPENMP_H
#define OPENMP_H

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

    int call_omp_get_max_threads();

    int call_omp_get_num_procs();

    void call_omp_set_num_threads(int num_threads);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // OPENMP_H
