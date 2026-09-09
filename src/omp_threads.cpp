// Set the number of OpenMP threads for this process.
//
// The on-disk linear algebra is performed by BigDataStatMeth's C++ code, which
// reads the number of threads to use from the OpenMP runtime on every call
// (omp_get_max_threads()) rather than from an argument. Calling
// omp_set_num_threads() sets that maximum for the current process, so this is
// the supported way to size the algebra from R. Compiled without OpenMP the
// function is a no-op and the algebra runs single-threaded.

#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export(".mc_omp_set_threads")]]
int mc_omp_set_threads(int n) {
    if (n < 1) {
        n = 1;
    }
#ifdef _OPENMP
    omp_set_num_threads(n);
    return omp_get_max_threads();
#else
    return 1;
#endif
}
