#ifndef SPYRALPP_NUMPY_ARGSORT_H
#define SPYRALPP_NUMPY_ARGSORT_H
#include <vector>
#include <cstddef>
namespace rk4 {
// Equal-distance ordering is observable in sklearn's HDBSCAN border labels.
// Use the same pinned sorting kernel as NumPy 2.1.3 on x86 SIMD machines.
std::vector<size_t> NumpyArgsort(std::vector<double> values);
}
#endif
