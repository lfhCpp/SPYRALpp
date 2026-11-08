#include "NumpyArgsort.h"
#include <algorithm>
#include <numeric>
namespace rk4 {
#ifdef RK4_NUMPY_X86_ARGSORT
void NumpyArgsortAvx512(double*,size_t*,size_t);
void NumpyArgsortAvx2(double*,size_t*,size_t);
#endif
std::vector<size_t> NumpyArgsort(std::vector<double> values){
    std::vector<size_t> indices(values.size());std::iota(indices.begin(),indices.end(),0);
#ifdef RK4_NUMPY_X86_ARGSORT
    if(__builtin_cpu_supports("avx512f")&&__builtin_cpu_supports("avx512dq")&&__builtin_cpu_supports("avx512vl")){
        NumpyArgsortAvx512(values.data(),indices.data(),values.size());return indices;
    }
    if(__builtin_cpu_supports("avx2")){
        NumpyArgsortAvx2(values.data(),indices.data(),values.size());return indices;
    }
#endif
    // Portable fallback preserves a valid ordering. Border assignments at exact
    // ties may differ from SIMD NumPy; oracle tests make that visible.
    std::stable_sort(indices.begin(),indices.end(),[&](size_t a,size_t b){return values[a]<values[b];});
    return indices;
}
}
