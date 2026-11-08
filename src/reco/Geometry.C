#include "Geometry.h"

#include <cmath>

namespace rk4 {

bool FitCircle(const std::vector<double>& xs, const std::vector<double>& ys,
               double& cx, double& cy, double& radius) {
    const size_t n = xs.size();
    if (n < 3 || ys.size() != n) return false;

    double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
    double sxz = 0.0, syz = 0.0, sz = 0.0;

    for (size_t i = 0; i < n; ++i) {
        const double x = xs[i];
        const double y = ys[i];
        const double z = x * x + y * y;
        sx += x;
        sy += y;
        sxx += x * x;
        syy += y * y;
        sxy += x * y;
        sxz += x * z;
        syz += y * z;
        sz += z;
    }

    // 求解
    // | sxx sxy sx | |A|   |sxz|
    // | sxy syy sy | |B| = |syz|
    // | sx  sy  n  | |C|   | sz|
    // 圆：x^2 + y^2 = A x + B y + C
    const double m00 = sxx, m01 = sxy, m02 = sx;
    const double m10 = sxy, m11 = syy, m12 = sy;
    const double m20 = sx,  m21 = sy,  m22 = static_cast<double>(n);
    const double b0 = sxz, b1 = syz, b2 = sz;

    const double det = m00 * (m11 * m22 - m12 * m21) -
                       m01 * (m10 * m22 - m12 * m20) +
                       m02 * (m10 * m21 - m11 * m20);
    if (std::fabs(det) < 1e-12) return false;

    const double detA = b0 * (m11 * m22 - m12 * m21) -
                        m01 * (b1 * m22 - m12 * b2) +
                        m02 * (b1 * m21 - m11 * b2);
    const double detB = m00 * (b1 * m22 - m12 * b2) -
                        b0 * (m10 * m22 - m12 * m20) +
                        m02 * (m10 * b2 - b1 * m20);
    const double detC = m00 * (m11 * b2 - b1 * m21) -
                        m01 * (m10 * b2 - b1 * m20) +
                        b0 * (m10 * m21 - m11 * m20);

    const double A = detA / det;
    const double B = detB / det;
    const double C = detC / det;

    cx = 0.5 * A;
    cy = 0.5 * B;
    const double r2 = C + cx * cx + cy * cy;
    if (r2 <= 0.0) return false;
    radius = std::sqrt(r2);
    return std::isfinite(radius) && radius > 0.0;
}

bool FitCircleSpyral(const std::vector<double>& xs,const std::vector<double>& ys,
                     double& cx,double& cy,double& radius) {
    const size_t n=xs.size();if(n<3||ys.size()!=n)return false;
    double mx=0,my=0;
    for(size_t i=0;i<n;++i){mx+=xs[i];my+=ys[i];}mx/=n;my/=n;
    double uu=0,vv=0,uv=0,uuu=0,uvv=0,uuv=0,vvv=0;
    for(size_t i=0;i<n;++i){
        const double u=xs[i]-mx,v=ys[i]-my;
        uu+=u*u;vv+=v*v;uv+=u*v;uuu+=u*u*u;uvv+=u*v*v;uuv+=u*u*v;vvv+=v*v*v;
    }
    const double det=uu*vv-uv*uv;if(det==0||!std::isfinite(det))return false;
    const double bx=(uuu+uvv)*.5,by=(uuv+vvv)*.5;
    cx=mx+(bx*vv-by*uv)/det;cy=my+(by*uu-bx*uv)/det;
    radius=0;for(size_t i=0;i<n;++i)radius+=std::hypot(xs[i]-cx,ys[i]-cy);radius/=n;
    return std::isfinite(cx)&&std::isfinite(cy)&&std::isfinite(radius)&&radius>0;
}

} // namespace rk4
