#include "Estimator.h"

#include "catima/catima.h"
#include <utility>
#include "../core/Units.h"
#include "Geometry.h"
#include "Spline.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace rk4 {


namespace {

// 离散平滑样条：min Σ(f_i - y_i)^2 + λ Σ (f_{i+1} - 2 f_i + f_{i-1})^2
// 系统 (I + λ DᵀD) f = y 对称正定、半带宽 2，用共轭梯度求解。
bool SmoothingFitPentadiag(const std::vector<double>& y, double lambda,
                           std::vector<double>& out) {
    const int n = static_cast<int>(y.size());
    if (n < 5) return false;

    auto apply = [&](const std::vector<double>& f, std::vector<double>& r) {
        r.assign(n, 0.0);
        for (int i = 1; i + 1 < n; ++i) {
            const double d2 = f[i + 1] - 2.0 * f[i] + f[i - 1];
            // D^T 作用回 i-1, i, i+1
            r[i - 1] += lambda * d2;
            r[i]     += -2.0 * lambda * d2;
            r[i + 1] += lambda * d2;
        }
        for (int i = 0; i < n; ++i) r[i] += f[i];
    };

    std::vector<double> x(n, 0.0), r = y, pdir = r, ap(n, 0.0);
    double rsold = 0.0;
    for (int i = 0; i < n; ++i) rsold += r[i] * r[i];
    const double tol = 1e-12 * std::max(1.0, rsold);
    for (int it = 0; it < 400; ++it) {
        apply(pdir, ap);
        double pap = 0.0;
        for (int i = 0; i < n; ++i) pap += pdir[i] * ap[i];
        if (!(pap > 1e-30)) break;
        const double alpha = rsold / pap;
        for (int i = 0; i < n; ++i) { x[i] += alpha * pdir[i]; r[i] -= alpha * ap[i]; }
        double rsnew = 0.0;
        for (int i = 0; i < n; ++i) rsnew += r[i] * r[i];
        if (rsnew < tol) break;
        const double beta = rsnew / rsold;
        for (int i = 0; i < n; ++i) pdir[i] = r[i] + beta * pdir[i];
        rsold = rsnew;
    }
    out = x;
    return true;
}

} // namespace

namespace {

// 沿 z 逐点局部二次拟合（Savitzky-Golay 的不均匀间距版）：
// 对每个点取窗口内的 (z, x) 做二次最小二乘，再在中心点求值。
void LocalQuadraticSmooth(const std::vector<double>& zs, std::vector<double>& xs, int half_win) {
    const int n = static_cast<int>(zs.size());
    if (n < 5) return;
    std::vector<double> out(xs.size(), 0.0);
    for (int i = 0; i < n; ++i) {
        const int lo = std::max(0, i - half_win);
        const int hi = std::min(n - 1, i + half_win);
        double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
        double t0 = 0, t1 = 0, t2 = 0;
        for (int j = lo; j <= hi; ++j) {
            const double z = zs[static_cast<size_t>(j)] - zs[static_cast<size_t>(i)];
            const double y = xs[static_cast<size_t>(j)];
            const double w = 1.0;
            s0 += w;
            s1 += w * z;
            s2 += w * z * z;
            s3 += w * z * z * z;
            s4 += w * z * z * z * z;
            t0 += w * y;
            t1 += w * y * z;
            t2 += w * y * z * z;
        }
        // 3x3 正规方程解 a0（中心点值）
        const double m00 = s0, m01 = s1, m02 = s2;
        const double m11 = s2, m12 = s3, m22 = s4;
        const double det = m00 * (m11 * m22 - m12 * m12) - m01 * (m01 * m22 - m12 * m02) +
                           m02 * (m01 * m12 - m11 * m02);
        if (std::fabs(det) < 1e-12) {
            out[static_cast<size_t>(i)] = xs[static_cast<size_t>(i)];
            continue;
        }
        const double inv = 1.0 / det;
        const double a0 = ((m11 * m22 - m12 * m12) * t0 + (m02 * m12 - m01 * m22) * t1 +
                           (m01 * m12 - m02 * m11) * t2) * inv;
        out[static_cast<size_t>(i)] = a0;
    }
    xs.swap(out);
}

} // namespace

bool Estimator::Estimate(const Cluster& cluster, double bfield_t,
                         EstimateResult& out) const {
    return Estimate(cluster, PointCloud(), bfield_t, out);
}

bool Estimator::Estimate(const Cluster& cluster, const PointCloud& cloud, double bfield_t,
                         EstimateResult& out) const {
    if (use_reinsch_) return EstimateReinsch(cluster, bfield_t, out);
    out = EstimateResult{};

    std::vector<CloudPoint> pts = cluster.points;
    if (pts.size() < 10) return false;

    std::sort(pts.begin(), pts.end(),
              [](const CloudPoint& a, const CloudPoint& b) { return a.z < b.z; });

    // 可选：对 x(z)、y(z) 做平滑样条降噪（对应 Spyral 的 make_smoothing_spline）
    if (smooth_ && pts.size() >= 15) {
        std::vector<double> zs, xs, ys;
        zs.reserve(pts.size());
        xs.reserve(pts.size());
        ys.reserve(pts.size());
        for (const auto& q : pts) {
            zs.push_back(q.z);
            xs.push_back(q.x);
            ys.push_back(q.y);
        }
        if (SmoothClusterXY(zs, xs, ys, lambda_)) {
            for (size_t i = 0; i < pts.size(); ++i) {
                pts[i].z = zs[i];
                pts[i].x = xs[i];
                pts[i].y = ys[i];
            }
        }
    }

    // 让 index 0 是顶点：Backward 时反向
    if (cluster.direction == Direction::Backward) {
        std::reverse(pts.begin(), pts.end());
    }

    out.direction = static_cast<int>(cluster.direction);
    out.vertex_x = pts.front().x;
    out.vertex_y = pts.front().y;
    out.vertex_z = pts.front().z;

    const double vx = out.vertex_x;
    const double vy = out.vertex_y;

    // 第一段弧：取到顶点横向距离最大的位置。
    // 不用"第一次下降"来判定，因为 pad 量化会让 rho 在起始处抖动，
    // 直接把弧长截到十几个点，圆拟合会完全失真。
    size_t arc_end = 0;
    double rho_max = -1.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        const double rho = std::hypot(pts[i].x - vx, pts[i].y - vy);
        if (rho > rho_max) {
            rho_max = rho;
            arc_end = i;
        }
    }
    if (arc_end < 10) arc_end = pts.size() - 1;
    if (arc_end < 10) return false;

    // 首弧全长的终点（未被 152 mm 截断）：圆拟合与 polar 用它；
    // Spyral 的 152 mm 只用于"电荷与弧长"，不用来限制圆拟合。
    const size_t arc_full = arc_end;

    // 只积到小 pad 区边界（Spyral 的做法）：|(x,y)| > small_pad_radius 就截断，
    // 否则大小 pad 电荷刻度混在一起，dE/dx 会有系统偏差。
    if (small_pad_radius_mm_ > 0.0) {
        size_t cut = 0;
        while (cut <= arc_end &&
               std::hypot(pts[cut].x, pts[cut].y) <= small_pad_radius_mm_) {
            ++cut;
        }
        if (cut > 10) arc_end = cut - 1;
    }

    // 圆拟合
    std::vector<double> xs, ys;
    xs.reserve(arc_end + 1);
    ys.reserve(arc_end + 1);
    double charge = 0.0;
    for (size_t i = 0; i <= arc_end; ++i) {
        charge += pts[i].integral;
    }
    // 圆拟合用整条首弧：之前只用起始 40%，曲率被高估 → 半径偏小 → Bρ 偏小 3 倍
    const size_t fit_end = arc_full; // 整条首弧（整簇拟合实测 Bρ 更小：0.138 vs 0.154）
    for (size_t i = 0; i <= fit_end; ++i) {
        xs.push_back(pts[i].x);
        ys.push_back(pts[i].y);
    }
    out.charge = charge;

    double cx = 0.0, cy = 0.0, radius = 0.0;
    if (!FitCircle(xs, ys, cx, cy, radius)) return false;



    // 矢高法半径：R = L^2/(8h) + h/2。最小二乘圆在短弧上会低估半径
    // （实测 Bρ 偏低 ~23%），矢高法对短弧更稳，取两者的几何中位。
    if (xs.size() >= 10) {
        const double x0 = xs.front(), y0 = ys.front();
        const double x1 = xs.back(), y1 = ys.back();
        const double L = std::hypot(x1 - x0, y1 - y0);
        if (L > 5.0) {
            const double ux = (x1 - x0) / L, uy = (y1 - y0) / L;
            double h = 0.0;
            for (size_t i = 1; i + 1 < xs.size(); ++i) {
                const double dx = xs[i] - x0, dy = ys[i] - y0;
                const double perp = std::fabs(dx * uy - dy * ux);
                if (perp > h) h = perp;
            }
            if (h > 1e-3) {
                const double r_sag = L * L / (8.0 * h) + 0.5 * h;
                // 实验：不混入矢高（Spyral 只用最小二乘半径），看 Bρ 是否回到 0.255
                (void)r_sag;

                // 短弧/大曲率：最小二乘与矢高都可能偏，另算"三点外接圆半径中位"
                // （确定性的伪随机取样，避免每次结果不同）
                double r_trip = 0.0;
                {
                    const size_t m = xs.size();
                    std::vector<double> rs;
                    rs.reserve(200);
                    unsigned int seed = 12345u;
                    for (int it = 0; it < 200; ++it) {
                        seed = seed * 1664525u + 1013904223u;
                        const size_t i0 = seed % m;
                        seed = seed * 1664525u + 1013904223u;
                        const size_t i1 = seed % m;
                        seed = seed * 1664525u + 1013904223u;
                        const size_t i2 = seed % m;
                        if (i0 == i1 || i1 == i2 || i0 == i2) continue;
                        const double ax = xs[i0], ay = ys[i0];
                        const double bx = xs[i1], by = ys[i1];
                        const double cx2 = xs[i2], cy2 = ys[i2];
                        const double d = 2.0 * (ax * (by - cy2) + bx * (cy2 - ay) + cx2 * (ay - by));
                        if (std::fabs(d) < 1e-6) continue;
                        const double a2 = ax * ax + ay * ay;
                        const double b2 = bx * bx + by * by;
                        const double c2 = cx2 * cx2 + cy2 * cy2;
                        const double ux = (a2 * (by - cy2) + b2 * (cy2 - ay) + c2 * (ay - by)) / d;
                        const double uy = (a2 * (cx2 - bx) + b2 * (ax - cx2) + c2 * (bx - ax)) / d;
                        const double rr = std::hypot(ax - ux, ay - uy);
                        if (std::isfinite(rr) && rr > 0.0) rs.push_back(rr);
                    }
                    if (rs.size() > 20) {
                        std::sort(rs.begin(), rs.end());
                        r_trip = rs[rs.size() / 2];
                    }
                }
                if (r_trip > 0.0) {
                    // 三个估计取中位，压制单个估计的病态
                    std::vector<double> cand{radius, r_sag, r_trip};
                    std::sort(cand.begin(), cand.end());
                    radius = cand[1];
                }
            }
        }
    }
    out.center_x = cx;
    out.center_y = cy;
    out.radius = radius;

    // 极角：rho-顶点 对 z 的线性回归
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    // Spyral: 只用起始半段做 rho-vs-z 回归（test_index = max(10, maximum/2)）
    const size_t fit_points = std::max<size_t>(10, arc_full / 2);
    const size_t n = std::min(arc_full + 1, fit_points);
    for (size_t i = 0; i < n; ++i) {
        const double x = pts[i].z;
        const double y = std::hypot(pts[i].x - vx, pts[i].y - vy);
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    const double denom = n * sxx - sx * sx;
    if (std::fabs(denom) < 1e-9) return false;
    const double slope = (n * sxy - sx * sy) / denom;
    out.polar = std::atan(slope);
    // Backward 径迹沿 -z 走，极角应落在 (90°, 180°)
    if (out.polar < 0.0) out.polar += Units::PI;
    if (out.polar <= 0.0 || out.polar >= Units::PI) return false;

    //------------------------------------------------------------------
    // 顶点精修：对"到束流轴的距离 rho"关于 z 做线性回归，外推到 rho = 0。
    // 取簇首点当顶点误差一个 pad 量级（10-20 mm），会把 beam_region 切割全杀掉。
    //------------------------------------------------------------------
    double vx2 = vx;
    double vy2 = vy;
    {
        double sx2 = 0.0, sy2 = 0.0, sxx2 = 0.0, sxy2 = 0.0;
        for (size_t i = 0; i <= arc_end; ++i) {
            const double z = pts[i].z;
            const double rho = std::hypot(pts[i].x, pts[i].y); // 到束流轴
            sx2 += z;
            sy2 += rho;
            sxx2 += z * z;
            sxy2 += z * rho;
        }
        const double den2 = n * sxx2 - sx2 * sx2;
        if (std::fabs(den2) > 1e-9) {
            const double a2 = (n * sxy2 - sx2 * sy2) / den2;
            const double b2 = (sy2 - a2 * sx2) / n;
            if (std::fabs(a2) > 1e-6) {
                double zv = -b2 / a2;
                const double z_lo = std::min(pts[0].z, pts[arc_end].z) - 100.0;
                const double z_hi = std::max(pts[0].z, pts[arc_end].z) + 100.0;
                if (zv < z_lo) zv = z_lo;
                if (zv > z_hi) zv = z_hi;
                // 在 z = zv 处线性插值 x、y（沿轨迹）
                double xv = pts[0].x;
                double yv = pts[0].y;
                for (size_t i = 0; i + 1 <= arc_end; ++i) {
                    const double z0 = pts[i].z;
                    const double z1 = pts[i + 1].z;
                    if ((zv >= std::min(z0, z1)) && (zv <= std::max(z0, z1))) {
                        const double dz = z1 - z0;
                        const double t = (std::fabs(dz) > 1e-9) ? (zv - z0) / dz : 0.0;
                        xv = pts[i].x + t * (pts[i + 1].x - pts[i].x);
                        yv = pts[i].y + t * (pts[i + 1].y - pts[i].y);
                        break;
                    }
                }
                out.vertex_x = xv;
                out.vertex_y = yv;
                out.vertex_z = zv;
                vx2 = xv;
                vy2 = yv;
            }
        }
    }

    // 顶点必须在束流轴附近（Spyral: vertex_rho > beam_region_radius → 丢弃）
    const double vertex_rho = std::hypot(vx2, vy2);
    if (beam_region_radius_mm_ > 0.0 && vertex_rho > beam_region_radius_mm_) return false;

    // 方向自洽：Forward 应 polar < 90°，Backward 应 > 90°
    if (cluster.direction == Direction::Forward && out.polar >= 0.5 * 3.14159265358979323846)
        return false;
    if (cluster.direction == Direction::Backward && out.polar < 0.5 * 3.14159265358979323846)
        return false;
    // direction == None 时不做方向切（否则会把大半个点云砍掉）

    // 方位角：从圆心指向顶点，再转 90 度（与 Spyral 同口径）
    double az = std::atan2(vy2 - cy, vx2 - cx) + 0.5 * 3.14159265358979323846;
    while (az < 0.0) az += 2.0 * 3.14159265358979323846;
    while (az > 2.0 * 3.14159265358979323846) az -= 2.0 * 3.14159265358979323846;
    out.azimuthal = az;

    // Rigidity is a magnitude. Preserve the signed field at propagation
    // boundaries, but use |B| here so a negative coordinate convention does
    // not produce a negative brho and reject an otherwise valid estimate.
    const double sin_theta = std::sin(out.polar);
    if (sin_theta > 1e-6) {
        out.brho = std::fabs(bfield_t) * (radius * 0.001) / sin_theta; // Tm
    }

    // 弧长：先把轨迹平滑（论文强调的关键一步），再沿平滑后的轨迹细采样积分。
    // 直接用折线求和会被 pad 量化噪声放大长度 → dE/dx 偏小且散。
    double length = 0.0;
    {
        const int kSamples = 1000;
        const double z0 = pts[0].z;
        const double z1 = pts[arc_end].z;

        // 样条用整条簇拟合（Spyral: apply_smoothing_splines 作用于完整 cluster），
        // 积分区间再截到 first_arc —— 短弧时样条被整条轨迹约束，不会乱翘。
        std::vector<double> zs, xs, ys;
        zs.reserve(pts.size());
        xs.reserve(pts.size());
        ys.reserve(pts.size());
        for (size_t i = 0; i < pts.size(); ++i) {
            zs.push_back(pts[i].z);
            xs.push_back(pts[i].x);
            ys.push_back(pts[i].y);
        }
        {
            std::vector<size_t> ord(zs.size());
            for (size_t i = 0; i < ord.size(); ++i) ord[i] = i;
            std::sort(ord.begin(), ord.end(),
                      [&](size_t a, size_t b) { return zs[a] < zs[b]; });
            std::vector<double> z2(zs.size()), x2(zs.size()), y2(zs.size());
            for (size_t i = 0; i < ord.size(); ++i) {
                z2[i] = zs[ord[i]];
                x2[i] = xs[ord[i]];
                y2[i] = ys[ord[i]];
            }
            zs.swap(z2);
            xs.swap(x2);
            ys.swap(y2);
        }
        // Spyral 用 make_smoothing_spline（penalized least squares，节点在数据点上）。
        // 普通移动窗口平滑会改变曲线形状；这里用 P-spline（B 样条 + 二阶差分惩罚），
        // 是"最小二乘 + 曲率惩罚"，不会把径迹拉弯。λ = smooth_lambda（Spyral 用 100）。
        const int npts = static_cast<int>(zs.size());
        // 离散平滑样条（节点 = 数据点）：min Σ(f-y)² + λΣ(Δ²f)²
        std::vector<double> xs_s, ys_s;
        const bool spline_ok = (npts >= 12) && SmoothingFitPentadiag(xs, lambda_, xs_s) &&
                               SmoothingFitPentadiag(ys, lambda_, ys_s);
        if (spline_ok) {
            xs = xs_s;
            ys = ys_s;
        }
        if (!spline_ok) {
            const int half_win = (npts >= 20) ? 6 : (npts >= 10 ? 3 : 0);
            LocalQuadraticSmooth(zs, xs, half_win);
            LocalQuadraticSmooth(zs, ys, half_win);
        }

        // 采样：一律在（平滑后的）数据点上线性插值
        auto at = [&](double z, double& x, double& y) {
            size_t lo = 0;
            while (lo + 1 <= arc_end &&
                   ((z1 > z0) ? (zs[lo + 1] < z) : (zs[lo + 1] > z))) {
                ++lo;
            }
            const size_t hi = std::min(lo + 1, arc_end);
            const double dz = zs[hi] - zs[lo];
            const double t = (std::fabs(dz) > 1e-9) ? (z - zs[lo]) / dz : 0.0;
            x = xs[lo] + t * (xs[hi] - xs[lo]);
            y = ys[lo] + t * (ys[hi] - ys[lo]);
        };

        double xp = 0.0, yp = 0.0, zp = z0;
        at(z0, xp, yp);
        for (int s = 1; s <= kSamples; ++s) {
            const double z = z0 + (z1 - z0) * s / kSamples;
            double x = 0.0, y = 0.0;
            at(z, x, y);
            const double dx = x - xp, dy = y - yp, dz = z - zp;
            length += std::sqrt(dx * dx + dy * dy + dz * dz);
            xp = x;
            yp = y;
            zp = z;
        }
    }
    // 护栏：样条弧长与折线弧长差距过大时（>25%），说明过冲或混入杂点，退回折线长度
    {
        double chord = 0.0;
        for (size_t i = 0; i < arc_end; ++i) {
            const double dx = pts[i + 1].x - pts[i].x;
            const double dy = pts[i + 1].y - pts[i].y;
            const double dz = pts[i + 1].z - pts[i].z;
            chord += std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        if (chord > 1.0 && length > 1.25 * chord) length = chord;
    }

    out.arclength = length;
    if (length > 1.0) {
        out.dEdx = charge / length;

        // polar 保留上面的 atan(drho/dz) 回归（Spyral 口径）；几何定义先停用

    }

    //------------------------------------------------------------------
    // 沿路径的 Bethe 平均：用 catima 从顶点动量出发逐步降能，
    // 得到 ∫(dE/dx)ds / L（MeV/mm）。这是 PID 的正确参考量：
    // 与 dE/dx 一样是"路径平均"，自动吸收 β 沿程演化。
    //------------------------------------------------------------------
    if (target_ready_ && out.brho > 0.0 && length > 1.0) {
        const double mass = MassData::GetNuclearMassNuclear(1, 1) * Units::u_MeV;
        catima::Material mat({static_cast<double>(target_A_), target_Z_,
                              static_cast<double>(target_stoich_)});
        mat.density(target_density_);
        catima::Projectile proj(1, 1);

        double p = 299.792458 * out.brho; // MeV/c（Z = 1）
        const int n = 40;
        const double ds = length / n;
        double sum = 0.0;
        int used = 0;
        for (int i = 0; i < n; ++i) {
            const double E = std::sqrt(p * p + mass * mass) - mass;
            if (!(E > 0.5)) break;
            proj.T = E / mass;
            const double dd = catima::dedx(proj, mat) * target_density_ / 10.0; // MeV/mm
            if (!(dd > 0.0)) break;
            sum += dd;
            ++used;
            double Enew = E - dd * pred_strength_ * ds;
            if (Enew < 0.5) Enew = 0.5;
            p = std::sqrt(Enew * (Enew + 2.0 * mass));
        }
        if (used > 4) out.dedx_pred = sum / used;
    }

    //------------------------------------------------------------------
    // dE/dx：pad 足迹并集 + 电荷收集（charge_radius_mm_ > 0 时启用）
    //   每个命中点在弧长上覆盖 [s-p/2, s+p/2]（p = pad 间距），分母取并集长度；
    //   分子再加簇外邻 pad 的电荷（补回阈值以下丢掉的边缘电荷）。
    //   这样"没有读出的空段"既不进分子也不进分母，不同长度的径迹才可比。
    //------------------------------------------------------------------
    if (charge_radius_mm_ > 0.0) {
        std::vector<double> s(arc_end + 1, 0.0);
        for (size_t i = 1; i <= arc_end; ++i) {
            const double dx = pts[i].x - pts[i - 1].x;
            const double dy = pts[i].y - pts[i - 1].y;
            const double dz = pts[i].z - pts[i - 1].z;
            s[i] = s[i - 1] + std::sqrt(dx * dx + dy * dy + dz * dz);
        }

        const double half_pitch = 0.5 * pad_pitch_mm_;
        std::vector<std::pair<double, double>> spans;
        double q_total = 0.0;
        for (size_t i = 0; i <= arc_end; ++i) {
            spans.emplace_back(s[i] - half_pitch, s[i] + half_pitch);
            q_total += pts[i].integral;
        }

        if (!cloud.points.empty()) {
            for (const auto& cp : cloud.points) {
                bool in_cluster = false;
                double best_d2 = 1e300;
                size_t best_i = 0;
                for (size_t i = 0; i <= arc_end; ++i) {
                    const double dx = cp.x - pts[i].x;
                    const double dy = cp.y - pts[i].y;
                    const double dz = cp.z - pts[i].z;
                    const double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 < best_d2) {
                        best_d2 = d2;
                        best_i = i;
                    }
                    if (d2 < 1e-6) {
                        in_cluster = true;
                        break;
                    }
                }
                if (in_cluster) continue;
                if (std::sqrt(best_d2) <= charge_radius_mm_) {
                    q_total += cp.integral;
                    spans.emplace_back(s[best_i] - half_pitch, s[best_i] + half_pitch);
                }
            }
        }

        std::sort(spans.begin(), spans.end());
        double covered = 0.0;
        double cur_lo = spans.front().first;
        double cur_hi = spans.front().second;
        for (size_t i = 1; i < spans.size(); ++i) {
            if (spans[i].first > cur_hi) {
                covered += (cur_hi - cur_lo);
                cur_lo = spans[i].first;
                cur_hi = spans[i].second;
            } else {
                cur_hi = std::max(cur_hi, spans[i].second);
            }
        }
        covered += (cur_hi - cur_lo);
        if (covered > 1e-6) {
            out.charge = q_total;
            out.arclength = covered;
            out.dEdx = q_total / covered;
        }
    }

    out.ok = std::isfinite(out.brho) && out.brho > 0.0;
    return out.ok;
}

} // namespace rk4
