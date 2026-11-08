#ifndef SPYRALPP_REINSCH_H
#define SPYRALPP_REINSCH_H
#include <vector>
#include <cmath>

namespace rk4 {
// Exact natural cubic smoothing: solve (R + lambda Q Q^T) u = Q y,
// then f = y - lambda Q^T u. The dual system is pentadiagonal;
// factorization, solve and storage are all O(n), without forming R^-1.
inline bool ReinschSmooth(const std::vector<double>& x,
                          const std::vector<double>& y, double lambda,
                          std::vector<double>& out,
                          std::vector<double>* second_derivatives = nullptr) {
    const size_t n = x.size();
    if (n < 5 || y.size() != n || !std::isfinite(lambda) || lambda < 0) return false;
    std::vector<long double> h(n-1);
    for (size_t i=0; i<n; ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return false;
        if (i) {
            h[i-1] = x[i]-x[i-1];
            if (!(h[i-1] > 0)) return false;
        }
    }
    if (lambda == 0 && !second_derivatives) { out=y; return true; }
    const size_t m=n-2;
    std::vector<long double> a(m), b(m), c(m), diag(m), l1(m), l2(m), rhs(m);
    for (size_t i=0; i<m; ++i) {
        a[i]=1/h[i]; b[i]=-1/h[i]-1/h[i+1]; c[i]=1/h[i+1];
        rhs[i]=a[i]*y[i]+b[i]*y[i+1]+c[i]*y[i+2];
    }
    for (size_t i=0; i<m; ++i) {
        if (i>=2) l2[i]=lambda*a[i]*c[i-2]/diag[i-2];
        if (i) {
            long double off=h[i]/6 + lambda*(a[i]*b[i-1]+b[i]*c[i-1]);
            if (i>=2) off-=l2[i]*diag[i-2]*l1[i-1];
            l1[i]=off/diag[i-1];
        }
        diag[i]=(h[i]+h[i+1])/3+lambda*(a[i]*a[i]+b[i]*b[i]+c[i]*c[i]);
        if (i) diag[i]-=l1[i]*l1[i]*diag[i-1];
        if (i>=2) diag[i]-=l2[i]*l2[i]*diag[i-2];
        if (!(diag[i]>0) || !std::isfinite(diag[i])) return false;
        if (i) rhs[i]-=l1[i]*rhs[i-1];
        if (i>=2) rhs[i]-=l2[i]*rhs[i-2];
    }
    std::vector<long double> u(m);
    for (size_t i=m; i-- > 0;) {
        u[i]=rhs[i]/diag[i];
        if (i+1<m) u[i]-=l1[i+1]*u[i+1];
        if (i+2<m) u[i]-=l2[i+2]*u[i+2];
    }
    std::vector<long double> result(y.begin(),y.end());
    for (size_t i=0; i<m; ++i) {
        result[i]-=lambda*a[i]*u[i];
        result[i+1]-=lambda*b[i]*u[i];
        result[i+2]-=lambda*c[i]*u[i];
    }
    for (double f:result) if (!std::isfinite(f)) return false;
    out.assign(result.begin(), result.end());
    if (second_derivatives) {
        second_derivatives->assign(n,0.0);
        for (size_t i=0;i<m;++i) (*second_derivatives)[i+1]=static_cast<double>(u[i]);
    }
    return true;
}
} // namespace rk4
#endif
