// 5000 事件：用 Reinsch 式自然三次平滑样条（节点=数据点，二阶差分惩罚 + 共轭梯度）
// 替代 P-spline，检验 80-150 点区间的 IQR 是否降到 Spyral 水平
R__ADD_INCLUDE_PATH(src)
R__ADD_INCLUDE_PATH(extern/eigen)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/Estimator.h"
#include "reco/Geometry.h"
#include "reco/ParticleId.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#ifndef GUARD_LIMIT
#define GUARD_LIMIT 0.30
#endif
using namespace rk4;

struct P { double x, y, z, q; };

// ---- 平滑节点值：min Σ(f-y)² + λΣ(Δ²f)² ----
static bool SmoothNodes(const std::vector<double>& y, double lam, std::vector<double>& out) {
    const int n = (int)y.size();
    if (n < 5) return false;
    auto apply = [&](const std::vector<double>& f, std::vector<double>& r) {
        r.assign(n, 0.0);
        for (int i = 1; i + 1 < n; ++i) {
            const double d2 = f[i + 1] - 2.0 * f[i] + f[i - 1];
            r[i - 1] += lam * d2;
            r[i] += -2.0 * lam * d2;
            r[i + 1] += lam * d2;
        }
        for (int i = 0; i < n; ++i) r[i] += f[i];
    };
    std::vector<double> x(n, 0.0), r = y, p = r, ap(n, 0.0);
    double rs = 0.0;
    for (int i = 0; i < n; ++i) rs += r[i] * r[i];
    for (int it = 0; it < 300; ++it) {
        apply(p, ap);
        double pap = 0.0;
        for (int i = 0; i < n; ++i) pap += p[i] * ap[i];
        if (!(pap > 1e-30)) break;
        const double a = rs / pap;
        for (int i = 0; i < n; ++i) { x[i] += a * p[i]; r[i] -= a * ap[i]; }
        double rn = 0.0;
        for (int i = 0; i < n; ++i) rn += r[i] * r[i];
        if (rn < 1e-20) break;
        const double b = rn / rs;
        for (int i = 0; i < n; ++i) p[i] = r[i] + b * p[i];
        rs = rn;
    }
    out = x;
    return true;
}

// ---- 自然三次样条（节点 = (x_i, f_i)），返回 1000 段积分 ----
static double CubicArclength(const std::vector<double>& x, const std::vector<double>& fx,
                             const std::vector<double>& fy, double za, double zb,
                             int samples = 1000) {
    const int n = (int)x.size();
    if (n < 3) return 0.0;
    auto secondDeriv = [&](const std::vector<double>& f, std::vector<double>& m) {
        m.assign(n, 0.0);
        std::vector<double> h(n - 1);
        for (int i = 0; i + 1 < n; ++i) h[i] = std::max(0.5, x[i + 1] - x[i]);  // h 下限 0.5mm（防病态）
        std::vector<double> a(n - 2), b(n - 2), c(n - 2), d(n - 2);
        for (int i = 1; i + 1 < n; ++i) {
            a[i - 1] = h[i - 1];
            b[i - 1] = 2.0 * (h[i - 1] + h[i]);
            c[i - 1] = h[i];
            d[i - 1] = 6.0 * ((f[i + 1] - f[i]) / h[i] - (f[i] - f[i - 1]) / h[i - 1]);
        }
        // 自然边界：m_0 = m_{n-1} = 0 ⇒ 首行去掉 m_0 的系数、末行去掉 m_{n-1} 的系数
        if (n - 2 > 0) {
            a[0] = 0.0;
            c[n - 3] = 0.0;
        }
        // Thomas
        const int k = n - 2;
        for (int i = 1; i < k; ++i) {
            const double w = a[i] / b[i - 1];
            b[i] -= w * c[i - 1];
            d[i] -= w * d[i - 1];
        }
        std::vector<double> s(k);
        if (k > 0) s[k - 1] = d[k - 1] / b[k - 1];
        for (int i = k - 2; i >= 0; --i) s[i] = (d[i] - c[i] * s[i + 1]) / b[i];
        for (int i = 0; i < k; ++i) m[i + 1] = s[i];
        // 钳制：|m| 不超过数据二阶差分的量级（抑制病态解）
        double scale = 0.0;
        for (int i = 1; i + 1 < n; ++i) {
            const double d2 = std::fabs(f[i + 1] - 2.0 * f[i] + f[i - 1]) / (h[i] * h[i]);
            if (d2 > scale) scale = d2;
        }
        const double mmax = 20.0 * std::max(scale, 1e-6);
        for (int i = 0; i < n; ++i) {
            if (m[i] > mmax) m[i] = mmax;
            if (m[i] < -mmax) m[i] = -mmax;
        }
    };
    std::vector<double> mx, my;
    secondDeriv(fx, mx);
    secondDeriv(fy, my);
#ifdef REINSCH_ZERO_M
    for (auto& v : mx) v = 0.0;
    for (auto& v : my) v = 0.0;
#endif
    auto eval = [&](const std::vector<double>& f, const std::vector<double>& m, double z) {
        if (z <= x[0]) return f[0];
        if (z >= x[n - 1]) return f[n - 1];
        int i = 0;
        while (i + 2 < n && x[i + 1] < z) ++i;
        const double h = x[i + 1] - x[i];
        const double A = (x[i + 1] - z) / h, B = (z - x[i]) / h;
        return A * f[i] + B * f[i + 1] +
               ((A * A * A - A) * m[i] + (B * B * B - B) * m[i + 1]) * h * h / 6.0;
    };
    double len = 0.0;
    double xp = eval(fx, mx, za), yp = eval(fy, my, za), zp = za;
    for (int i = 1; i <= samples; ++i) {
        const double z = za + (zb - za) * i / samples;
        const double xx = eval(fx, mx, z), yy = eval(fy, my, z);
        len += std::sqrt((xx - xp) * (xx - xp) + (yy - yp) * (yy - yp) + (z - zp) * (z - zp));
        xp = xx; yp = yy; zp = z;
    }
    return len;
}


// k=2 拆分：仅在中心距 > 3×簇内 RMS 且两半都 ≥25 点时返回 true
static bool SplitKMeans2(const std::vector<P>& pts, std::vector<P>& a, std::vector<P>& b) {
    const size_t n = pts.size();
    (void)a; (void)b;
    return false;  // 拆分实测无改善，关闭
    // 初值：最远点对（确定性）
    size_t i0 = 0, i1 = 0;
    double best = -1;
    for (size_t i = 0; i < n; i += 7) {
        for (size_t j = i + 1; j < n; j += 11) {
            const double dx = pts[i].x - pts[j].x, dy = pts[i].y - pts[j].y, dz = pts[i].z - pts[j].z;
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 > best) { best = d2; i0 = i; i1 = j; }
        }
    }
    double c1[3] = {pts[i0].x, pts[i0].y, pts[i0].z};
    double c2[3] = {pts[i1].x, pts[i1].y, pts[i1].z};
    std::vector<int> lab(n, 0);
    for (int it = 0; it < 20; ++it) {
        double s1[3] = {0, 0, 0}, s2[3] = {0, 0, 0};
        int n1 = 0, n2 = 0;
        for (size_t i = 0; i < n; ++i) {
            const double d1 = (pts[i].x - c1[0]) * (pts[i].x - c1[0]) + (pts[i].y - c1[1]) * (pts[i].y - c1[1]) + (pts[i].z - c1[2]) * (pts[i].z - c1[2]);
            const double d2 = (pts[i].x - c2[0]) * (pts[i].x - c2[0]) + (pts[i].y - c2[1]) * (pts[i].y - c2[1]) + (pts[i].z - c2[2]) * (pts[i].z - c2[2]);
            if (d1 <= d2) { lab[i] = 1; s1[0] += pts[i].x; s1[1] += pts[i].y; s1[2] += pts[i].z; ++n1; }
            else { lab[i] = 2; s2[0] += pts[i].x; s2[1] += pts[i].y; s2[2] += pts[i].z; ++n2; }
        }
        if (n1 == 0 || n2 == 0) return false;
        for (int k = 0; k < 3; ++k) { c1[k] = s1[k] / n1; c2[k] = s2[k] / n2; }
    }
    // 分离判据
    const double sep = std::sqrt((c1[0] - c2[0]) * (c1[0] - c2[0]) + (c1[1] - c2[1]) * (c1[1] - c2[1]) +
                                 (c1[2] - c2[2]) * (c1[2] - c2[2]));
    double w1 = 0.0, w2 = 0.0;
    int n1 = 0, n2 = 0;
    for (size_t i = 0; i < n; ++i) {
        if (lab[i] == 1) {
            w1 += (pts[i].x - c1[0]) * (pts[i].x - c1[0]) + (pts[i].y - c1[1]) * (pts[i].y - c1[1]) + (pts[i].z - c1[2]) * (pts[i].z - c1[2]);
            ++n1;
        } else {
            w2 += (pts[i].x - c2[0]) * (pts[i].x - c2[0]) + (pts[i].y - c2[1]) * (pts[i].y - c2[1]) + (pts[i].z - c2[2]) * (pts[i].z - c2[2]);
            ++n2;
        }
    }
    if (n1 < 25 || n2 < 25) return false;
    const double w = std::sqrt((w1 + w2) / (n1 + n2));
    if (sep < 2.0 * std::max(w, 1e-6)) return false;
    for (size_t i = 0; i < n; ++i) (lab[i] == 1 ? a : b).push_back(pts[i]);
    return true;
}


// GCV 选平滑参数：GCV(λ) = (RSS/n) / (1 - tr(H)/n)²
// tr(H) 用 Rademacher 向量做一次随机估计（H = (I+λDᵀD)^{-1}）
static double GcvLambda(const std::vector<double>& y) {
    const int n = (int)y.size();
    if (n < 12) return 1000.0;
    std::vector<double> u(n);
    unsigned int s = 12345u;
    for (int i = 0; i < n; ++i) {
        s = s * 1664525u + 1013904223u;
        u[i] = ((s >> 16) & 1u) ? 1.0 : -1.0;
    }
    const double grid[] = {10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0};
    double best = 1e300, bestlam = 1000.0;
    for (double lam : grid) {
        std::vector<double> f, hu;
        if (!SmoothNodes(y, lam, f)) continue;
        if (!SmoothNodes(u, lam, hu)) continue;
        double tr = 0.0, rss = 0.0;
        for (int i = 0; i < n; ++i) {
            tr += hu[i] * u[i];
            const double r = f[i] - y[i];
            rss += r * r;
        }
        double den = 1.0 - tr / n;
        if (den < 0.05) den = 0.05;
        const double gcv = (rss / n) / (den * den);
        if (gcv < best) { best = gcv; bestlam = lam; }
    }
    return bestlam;
}


#ifdef SPYRALPP_FAST_REINSCH
#include "reco/Reinsch.h"
#else
// ---- 精确三次平滑样条（Reinsch），用 Eigen 稀疏求解 ----
// min Σ(f_i-y_i)² + lam·∫(f'')²  ⇒  (I + lam P) f = y,  P = Qᵀ R⁻¹ Q（五对角）
static bool ReinschSmooth(const std::vector<double>& x, const std::vector<double>& y, double lam,
                          std::vector<double>& out) {
    const int n = (int)x.size();
    if (n < 5) return false;
    std::vector<double> h(n - 1);
    for (int i = 0; i + 1 < n; ++i) h[i] = std::max(1e-6, x[i + 1] - x[i]);
    const int m = n - 2;
    // R 三对角
    std::vector<double> ra(m, 0.0), rb(m, 0.0), rc(m, 0.0);
    for (int i = 0; i < m; ++i) {
        ra[i] = (i > 0) ? h[i] / 6.0 : 0.0;
        rb[i] = (h[i] + h[i + 1]) / 3.0;
        rc[i] = (i + 1 < m) ? h[i + 1] / 6.0 : 0.0;
    }
    {
        static int flag = 0;
        if (flag == 0) {   // DUMP_H
            FILE* fh = std::fopen("/tmp/rs_h.csv", "w");
            std::fprintf(fh, "h\n");
            for (int i = 0; i + 1 < n; ++i) std::fprintf(fh, "%.9f\n", h[i]);
            std::fclose(fh);
            ++flag;
        }
    }
    // Z = R⁻¹ Q，逐列用 Eigen 三对角 LDLT
    Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(m, n);
    {
        Eigen::SparseMatrix<double> R(m, m);
        std::vector<Eigen::Triplet<double>> trips;
        for (int i = 0; i < m; ++i) {
            trips.emplace_back(i, i, rb[i]);
            if (i > 0) trips.emplace_back(i, i - 1, ra[i]);
            if (i + 1 < m) trips.emplace_back(i, i + 1, rc[i]);
        }
        R.setFromTriplets(trips.begin(), trips.end());
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(R);
        if (solver.info() != Eigen::Success) return false;
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(m);
        for (int col = 0; col < n; ++col) {
            rhs.setZero();
            if (col < m) rhs[col] += 1.0 / h[col];
            if (col - 1 >= 0 && col - 1 < m) {
                const double s1 = (col < n - 1) ? 1.0 / h[col] : 0.0;
                rhs[col - 1] += -(s1 + 1.0 / h[col - 1]);
            }
            if (col - 2 >= 0 && col - 2 < m) rhs[col - 2] += 1.0 / h[col - 1];
            Z.col(col) = solver.solve(rhs);
        }
    }
    // P = Qᵀ Z（稀疏）
    Eigen::SparseMatrix<double> P(n, n);
    {
        std::vector<Eigen::Triplet<double>> trips;
        for (int r = 0; r < m; ++r) {
            const double q[3] = {1.0 / h[r], -(1.0 / h[r] + 1.0 / h[r + 1]), 1.0 / h[r + 1]};
            const int idx[3] = {r, r + 1, r + 2};
            for (int a = 0; a < 3; ++a)
                for (int j = 0; j < n; ++j)
                    if (std::fabs(Z(r, j)) > 0.0)
                        trips.emplace_back(idx[a], j, q[a] * Z(r, j));
        }
        P.setFromTriplets(trips.begin(), trips.end());
    }
    // (I + lam P) f = y
    Eigen::SparseMatrix<double> A = P;
    A *= lam;
    for (int i = 0; i < n; ++i) A.coeffRef(i, i) += 1.0;
    A.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    lu.compute(A);
    if (lu.info() != Eigen::Success) return false;
    Eigen::VectorXd yy(n);
    for (int i = 0; i < n; ++i) yy[i] = y[i];
    Eigen::VectorXd f = lu.solve(yy);
    if (lu.info() != Eigen::Success) return false;
    out.resize(n);
    for (int i = 0; i < n; ++i) out[i] = f[i];
    return true;
}

#endif

static bool Port(const std::vector<P>& raw, double B, EstimateResult& out) {
    std::vector<P> d = raw;
    if (d.size() < 10) return false;
    std::sort(d.begin(), d.end(), [](const P& a, const P& b) { return a.z < b.z; });
    const size_t n = d.size();
#ifdef SPYRALPP_PRESMOOTH
    {
        std::vector<double> z,x,y,q,fx,fy,fq;
        for (auto& point:d) {z.push_back(point.z);x.push_back(point.x);y.push_back(point.y);q.push_back(point.q);}
        if (!ReinschSmooth(z,x,100,fx) || !ReinschSmooth(z,y,100,fy) || !ReinschSmooth(z,q,100,fq)) return false;
        for (size_t i=0; i<n; ++i) { d[i].x=fx[i];d[i].y=fy[i];
#ifndef SPYRALPP_RAW_CHARGE
            d[i].q=fq[i];
#endif
 }
    }
#endif
    // 定向：离束流轴更近的一端作为顶点（Spyral 用 direction 翻转做同样的事）
    {
        const double r0 = std::hypot(d.front().x, d.front().y);
        const double r1 = std::hypot(d.back().x, d.back().y);
        if (r1 < 0.8 * r0) std::reverse(d.begin(), d.end());
    }
    double vx = d[0].x, vy = d[0].y;
    auto max_idx = [&](double cx0, double cy0) {
        double best = -1; size_t idx = 0;
        for (size_t i = 1; i < n; ++i) {
            const double r = std::hypot(d[i].x - cx0, d[i].y - cy0);
            if (r > best) { best = r; idx = i; }
        }
        return idx;
    };
    size_t maximum = max_idx(vx, vy);
    std::vector<double> xs, ys;
    for (size_t i = 0; i <= maximum; ++i) { xs.push_back(d[i].x); ys.push_back(d[i].y); }
    double cx = 0, cy = 0, R = 0;
    if (!FitCircle(xs, ys, cx, cy, R)) return false;
    // 圆拟合质量：点到圆的 RMSE / R（>0.3 说明曲率不可信）
    {
        double sse = 0.0;
        for (size_t i = 0; i <= maximum; ++i) {
            const double r = std::hypot(d[i].x - cx, d[i].y - cy);
            const double dd = r - R;
            sse += dd * dd;
        }
        const double rmse = std::sqrt(sse / std::max<size_t>(1, maximum + 1));
        if (R > 1.0 && rmse / R > 0.30) return false;
    }
    const double rc = std::hypot(cx, cy);
    if (rc > 1e-9) { vx = cx - R * cx / rc; vy = cy - R * cy / rc; }
#ifdef SPYRALPP_REFIT_VERTEX
    for (int iteration = 0; iteration < 2; ++iteration) {
        maximum = max_idx(vx, vy);
        xs.clear(); ys.clear();
        for (size_t i = 0; i <= maximum; ++i) { xs.push_back(d[i].x); ys.push_back(d[i].y); }
        if (!FitCircle(xs, ys, cx, cy, R)) return false;
        const double center_rho = std::hypot(cx, cy);
        if (center_rho > 1e-9) { vx = cx - R*cx/center_rho; vy = cy - R*cy/center_rho; }
    }
#endif
    maximum = max_idx(vx, vy);
    if (maximum < 10) return false;
    const size_t test_index = std::max<size_t>(10, maximum / 2);
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (size_t i = 0; i < test_index; ++i) {
        const double zz = d[i].z, rr = std::hypot(d[i].x - vx, d[i].y - vy);
        sx += zz; sy += rr; sxx += zz * zz; sxy += zz * rr;
    }
    const double den = test_index * sxx - sx * sx;
    if (std::fabs(den) < 1e-9) return false;
    const double slope = (test_index * sxy - sx * sy) / den;
    if (std::fabs(slope) < 1e-9) return false;
    double polar = std::atan(slope);
    if (polar < 0) polar += M_PI;
    // Spyral 的两道 cut（移植时漏了，会导致"接受太多低质量簇"→ 图上散点更多）
    const double vertex_rho = std::hypot(vx, vy);
    if (vertex_rho > 25.0) return false;                 // beam_region_radius
    if (polar <= 0.0 || polar >= M_PI) return false;
    if (polar >= 0.5 * M_PI && d.size() > 0) {
        // Forward/Backward 自洽：本移植版按 z 升序固定，等价于 Forward
        // 若 polar 落在 (90°,180°) 说明方向不自洽
        if (polar >= 0.9 * M_PI) return false;
    }
    double charge = d[0].q;
    size_t cutoff = maximum;
    for (size_t i = 0; i + 1 <= maximum; ++i) {
        if (std::hypot(d[i + 1].x, d[i + 1].y) > 152.0) { cutoff = i + 1; break; }
        charge += d[i + 1].q;
    }
    std::vector<double> zs, fx, fy;
    for (const auto& p : d) { zs.push_back(p.z); fx.push_back(p.x); fy.push_back(p.y); }
    // 去重：同一时间桶会产生相同 z，h=0 会让三阶样条出现 inf
    std::vector<double> uz, ufx, ufy;
    for (size_t i = 0; i < zs.size(); ++i) {
        if (!uz.empty() && zs[i] - uz.back() < 1e-6) {
            ufx.back() = 0.5 * (ufx.back() + fx[i]);
            ufy.back() = 0.5 * (ufy.back() + fy[i]);
            continue;
        }
        uz.push_back(zs[i]);
        ufx.push_back(fx[i]);
        ufy.push_back(fy[i]);
    }
    std::vector<double> gx, gy;
    // 精确 Reinsch 平滑样条，lam 与 Spyral 的 smoothing_factor 相同（100）
#ifdef SPYRALPP_PRESMOOTH
    gx=ufx; gy=ufy;
#else
    if (!ReinschSmooth(uz, ufx, 100.0, gx)) return false;
    if (!ReinschSmooth(uz, ufy, 100.0, gy)) return false;
#endif
    {
        static int dbg = 0;
        if (dbg == 0) {   // DUMP_RS：第一个簇存 CSV
            FILE* fd = std::fopen("/tmp/rs_dump.csv", "w");
            std::fprintf(fd, "z,xraw,xsmooth\n");
            for (size_t i = 0; i < uz.size(); ++i)
                std::fprintf(fd, "%.6f,%.6f,%.6f\n", uz[i], ufx[i], gx[i]);
            std::fclose(fd);
        }
        if (dbg < 5) {
            double mx = 0.0;
            for (size_t i = 0; i < uz.size(); ++i) mx = std::max(mx, std::fabs(gx[i] - ufx[i]));
            double chord = 0.0;
            for (size_t i = 0; i + 1 < uz.size(); ++i)
                chord += std::sqrt((ufx[i+1]-ufx[i])*(ufx[i+1]-ufx[i]) + (ufy[i+1]-ufy[i])*(ufy[i+1]-ufy[i]) + (uz[i+1]-uz[i])*(uz[i+1]-uz[i]));
            std::printf("DBG_RS n=%zu maxf-y=%.2f chord=%.1f zspan=%.1f\n", uz.size(), mx,
                        chord, std::fabs(uz.back() - uz.front()));
            ++dbg;
        }
    }
    const double len = CubicArclength(uz, gx, gy, d[0].z, d[cutoff].z);
    {
        // DBG_LEN：与折线长度对照
        double lin = 0.0;
        for (size_t i = 0; i + 1 <= cutoff; ++i) {
            const double dx = d[i + 1].x - d[i].x, dy = d[i + 1].y - d[i].y;
            const double dz = d[i + 1].z - d[i].z;
            lin += std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        static int dbg = 0;
        if (dbg < 5) {
            std::printf("DBG_LEN n=%zu charge=%.1f cubic=%.1f linear=%.1f zspan=%.1f\n", n, charge,
                        len, lin, std::fabs(d[cutoff].z - d[0].z));
            ++dbg;
        }
    }
    // 全局弧长标定：与 Spyral 的 dE/dx 比中位 0.871 → 弧长乘 1/0.871
    const double len_cal = len;   // 正确样条后不再需要 0.871 标定
    if (!(len_cal > 1.0)) return false;
    // 用平滑后的路径重拟合圆（抑制抖动 → Bρ 更稳）
    const double R_raw = R;
#ifndef SPYRALPP_PRESMOOTH
    {
        std::vector<double> px, py;
        for (size_t i = 0; i < uz.size() && i <= cutoff; ++i) { px.push_back(gx[i]); py.push_back(gy[i]); }
        double cx3 = 0, cy3 = 0, R3 = 0;
        if (px.size() >= 10 && FitCircle(px, py, cx3, cy3, R3) && R3 > 0) {
            double sse = 0.0;
            for (size_t i = 0; i < px.size(); ++i) {
                const double r = std::hypot(px[i] - cx3, py[i] - cy3);
                sse += (r - R3) * (r - R3);
            }
            const double rmse = std::sqrt(sse / px.size());
            if (rmse / R3 > 0.10) return false;   // 曲率不可信
#ifdef SPYRALPP_HELIX_ANGLE
            double sz=0, sp=0, szz=0, szp=0, phase=0, previous=0;
            for (size_t i=0; i<px.size(); ++i) {
                const double current=std::atan2(py[i]-cy3,px[i]-cx3);
                if (!i) phase=current;
                else phase+=std::remainder(current-previous,2*M_PI);
                previous=current;
                const double zz=uz[i];
                sz+=zz; sp+=phase; szz+=zz*zz; szp+=zz*phase;
            }
            const double denom=px.size()*szz-sz*sz;
            if (std::fabs(denom)>1e-9) {
                const double rate=(px.size()*szp-sz*sp)/denom;
                if (std::fabs(rate)>1e-9) polar=std::atan(R3*std::fabs(rate));
            }
#endif

                        R = std::sqrt(std::max(1e-6, R_raw) * R3);   // 两个独立半径的几何平均
        }
    }
#endif
    out.brho = B * (R * 0.001) / std::max(1e-6, std::fabs(std::sin(polar)));
    out.polar = polar;

    // 截断均值 dE/dx：沿首弧（到 152mm 截断）按 ~20 mm 分段，
    // 每段 dE/dx = 段电荷/段长，丢掉最高 30% 后平均（抗 Landau 尾部）。
    double dedx_trunc = 0.0;
    {
        double stot = 0.0;
        for (size_t i = 1; i <= cutoff; ++i) {
            const double dx = d[i].x - d[i - 1].x, dy = d[i].y - d[i - 1].y,
                         dz = d[i].z - d[i - 1].z;
            stot += std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        if (stot > 20.0) {
            // 用平滑后的弧长做标尺（折线被 pad 抖动拉长约 5 倍）
            const double scale_s = len_cal / stot;
            const double Ls = len_cal;
            const int nb = std::max(1, (int)std::ceil(Ls / std::max(8.0, std::min(25.0, Ls / 10.0))));
            const double binlen = Ls / nb;
            std::vector<double> qb(nb, 0.0);
            std::vector<int> cb(nb, 0);
            double acc = 0.0;
            for (size_t i = 0; i <= cutoff; ++i) {
                const int b = std::min(nb - 1, (int)((acc * scale_s) / binlen));
                qb[b] += d[i].q;
                cb[b] += 1;
                if (i < cutoff) {
                    const double dx = d[i + 1].x - d[i].x, dy = d[i + 1].y - d[i].y,
                                 dz = d[i + 1].z - d[i].z;
                    acc += std::sqrt(dx * dx + dy * dy + dz * dz);
                }
            }
            std::vector<double> v;
            for (int b = 0; b < nb; ++b)
                if (cb[b] > 0) v.push_back(qb[b] / binlen);
            if (v.size() >= 2) {
                std::sort(v.begin(), v.end());
                const int keep = std::max(1, (int)(v.size() * 0.7));
                double ssum = 0.0;
                for (int i = 0; i < keep; ++i) ssum += v[i];
                dedx_trunc = ssum / keep;
            }
        }
    }
    out.arclength = len_cal;
    {
        const double d_simple = charge / len_cal;
        // 截断均值要按全局刻度归一到同量级（0.821 是实测中位比）
        const double d_tr = (dedx_trunc > 0.0) ? dedx_trunc / 0.821 : d_simple;
        out.dEdx = d_simple;   // 纯 charge/arclength（与 Spyral 一致）
    }
    return true;
}

void port_reinsch(const char* input = "../results/our_clusters5k.csv",
                  const char* output = "/tmp/port_reinsch5k.csv") {
    std::map<std::pair<int, int>, std::vector<P>> groups;
    std::ifstream in(input);
    if (!in) throw std::runtime_error("Cannot open cluster input");
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() < 6) continue;
        groups[{(int)v[0], (int)v[1]}].push_back(P{v[2], v[3], v[4], v[5]});
    }
    FILE* fc = std::fopen(output, "w");
    if (!fc) throw std::runtime_error("Cannot open estimate output");
    std::fprintf(fc, "event,cluster,brho,dedx,polar_deg\n");
    size_t nok = 0;
    std::vector<double> brho, dedx;
    for (auto& kv : groups) {
        if (kv.second.size() < 30) continue;
        // 递归拆分：400+ 点的簇可能不止两条径迹
        std::vector<std::vector<P>> queue{kv.second};
        bool any = false;
        while (!queue.empty()) {
            std::vector<P> cur = queue.back();
            queue.pop_back();
            std::vector<P> a, b;
            if (cur.size() >= 400 && SplitKMeans2(cur, a, b)) {
                any = true;
                queue.push_back(a);
                queue.push_back(b);
                continue;
            }
            if (cur.size() < 30) continue;
            EstimateResult r;
            if (!Port(cur, 2.85, r)) continue;
            std::fprintf(fc, "%d,%d,%.6f,%.6f,%.3f\n", kv.first.first, kv.first.second,
                         r.brho, r.dEdx, r.polar * 180 / M_PI);
            brho.push_back(r.brho); dedx.push_back(r.dEdx); ++nok;
        }
        (void)any;
        continue;
        EstimateResult r;
        if (!Port(kv.second, 2.85, r)) continue;
        std::fprintf(fc, "%d,%d,%.6f,%.6f,%.3f\n", kv.first.first, kv.first.second, r.brho,
                     r.dEdx, r.polar * 180 / M_PI);
        brho.push_back(r.brho); dedx.push_back(r.dEdx); ++nok;
    }
    std::fclose(fc);
    std::sort(brho.begin(), brho.end());
    std::sort(dedx.begin(), dedx.end());
    std::printf("Reinsch 版: N=%zu Bρ中位=%.3f dEdx中位=%.0f\n", nok,
                brho.empty() ? 0.0 : brho[brho.size() / 2],
                dedx.empty() ? 0.0 : dedx[dedx.size() / 2]);
}
