#include "Spline.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace rk4 {
namespace {

// Cox-de Boor 递归：第 i 个、order 阶 B 样条基在 x 处的值
double BSplineRec(int i, int order, const std::vector<double>& t, double x) {
    if (order == 1) {
        const int L = static_cast<int>(t.size()) - 1; // 最后区间下标
        if (i >= L) {
            // x 恰为最大值时，最后一个基函数取 1
            return (x >= t[static_cast<size_t>(L)] && i == L - 1) ? 1.0 : 0.0;
        }
        return (x >= t[static_cast<size_t>(i)] && x < t[static_cast<size_t>(i + 1)]) ? 1.0 : 0.0;
    }

    double sum = 0.0;
    const double d1 = t[static_cast<size_t>(i + order - 1)] - t[static_cast<size_t>(i)];
    if (d1 > 1e-12) {
        sum += (x - t[static_cast<size_t>(i)]) / d1 * BSplineRec(i, order - 1, t, x);
    }
    const double d2 = t[static_cast<size_t>(i + order)] - t[static_cast<size_t>(i + 1)];
    if (d2 > 1e-12) {
        sum += (t[static_cast<size_t>(i + order)] - x) / d2 *
               BSplineRec(i + 1, order - 1, t, x);
    }
    return sum;
}

// 稠密线性方程组求解（部分选主元高斯消元）
bool SolveLinear(std::vector<std::vector<double>>& a, std::vector<double>& b,
                 std::vector<double>& x) {
    const int n = static_cast<int>(b.size());
    if (n == 0 || static_cast<int>(a.size()) != n) return false;

    for (int col = 0; col < n; ++col) {
        // 选主元
        int piv = col;
        double best = std::fabs(a[static_cast<size_t>(col)][static_cast<size_t>(col)]);
        for (int r = col + 1; r < n; ++r) {
            const double v = std::fabs(a[static_cast<size_t>(r)][static_cast<size_t>(col)]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-14) return false;
        if (piv != col) {
            std::swap(a[static_cast<size_t>(piv)], a[static_cast<size_t>(col)]);
            std::swap(b[static_cast<size_t>(piv)], b[static_cast<size_t>(col)]);
        }

        const double diag = a[static_cast<size_t>(col)][static_cast<size_t>(col)];
        for (int r = col + 1; r < n; ++r) {
            const double factor = a[static_cast<size_t>(r)][static_cast<size_t>(col)] / diag;
            if (factor == 0.0) continue;
            for (int c = col; c < n; ++c) {
                a[static_cast<size_t>(r)][static_cast<size_t>(c)] -=
                    factor * a[static_cast<size_t>(col)][static_cast<size_t>(c)];
            }
            b[static_cast<size_t>(r)] -= factor * b[static_cast<size_t>(col)];
        }
    }

    x.assign(static_cast<size_t>(n), 0.0);
    for (int r = n - 1; r >= 0; --r) {
        double sum = b[static_cast<size_t>(r)];
        for (int c = r + 1; c < n; ++c) {
            sum -= a[static_cast<size_t>(r)][static_cast<size_t>(c)] * x[static_cast<size_t>(c)];
        }
        x[static_cast<size_t>(r)] = sum / a[static_cast<size_t>(r)][static_cast<size_t>(r)];
    }
    return true;
}

} // namespace

void SmoothingSpline::BuildKnots(int n_knots) {
    // 内部节点等距分布在 (xmin, xmax)
    const int n_int = std::max(1, n_knots);
    m_ = n_int + order_;
    knots_.assign(static_cast<size_t>(2 * order_ + n_int), xmin_);

    for (int i = 0; i < n_int; ++i) {
        const double f = static_cast<double>(i + 1) / static_cast<double>(n_int + 1);
        knots_[static_cast<size_t>(order_ + i)] = xmin_ + (xmax_ - xmin_) * f;
    }
    for (int i = 0; i < order_; ++i) {
        knots_[static_cast<size_t>(order_ + n_int + i)] = xmax_;
    }
}

double SmoothingSpline::Basis(int i, int order, double x) const {
    return BSplineRec(i, order, knots_, x);
}

double SmoothingSpline::Eval(double x) const {
    if (!ready_) return 0.0;
    if (x <= xmin_) {
        // 在左端只有第一个基非零且为 1
        double sum = 0.0;
        for (int i = 0; i < m_; ++i) {
            sum += coef_[static_cast<size_t>(i)] * Basis(i, order_, xmin_);
        }
        return sum;
    }
    if (x >= xmax_) {
        return coef_[static_cast<size_t>(m_ - 1)];
    }
    double sum = 0.0;
    for (int i = 0; i < m_; ++i) {
        sum += coef_[static_cast<size_t>(i)] * Basis(i, order_, x);
    }
    return sum;
}

std::vector<double> SmoothingSpline::Smooth(const std::vector<double>& xs) const {
    std::vector<double> out;
    out.reserve(xs.size());
    for (double x : xs) out.push_back(Eval(x));
    return out;
}

void SmoothingSpline::AddPenalty(std::vector<std::vector<double>>& a, double lambda) const {
    // 二阶差分算子作用在系数上：D2 c 的第 j 项 = c[j] - 2c[j+1] + c[j+2]
    // 惩罚项 lambda * D2^T D2（m x m）
    if (lambda <= 0.0) return;
    for (int j = 0; j + 2 < m_; ++j) {
        const int idx[3] = {j, j + 1, j + 2};
        const double w[3] = {1.0, -2.0, 1.0};
        for (int p = 0; p < 3; ++p) {
            for (int q = 0; q < 3; ++q) {
                a[static_cast<size_t>(idx[p])][static_cast<size_t>(idx[q])] +=
                    lambda * w[p] * w[q];
            }
        }
    }
}

bool SmoothingSpline::Solve(std::vector<std::vector<double>>& a, std::vector<double>& b,
                            std::vector<double>& out) const {
    return SolveLinear(a, b, out);
}

bool SmoothingSpline::Fit(const std::vector<double>& xs, const std::vector<double>& ys,
                          double lambda, int n_knots) {
    ready_ = false;
    coef_.clear();

    if (xs.size() < 5 || xs.size() != ys.size()) return false;

    xmin_ = xs.front();
    xmax_ = xs.back();
    if (!(xmax_ > xmin_)) return false;

    BuildKnots(std::min(n_knots, static_cast<int>(xs.size()) - 2));

    const size_t n = xs.size();

    // 设计矩阵 B (n x m)
    std::vector<std::vector<double>> bmat(n, std::vector<double>(static_cast<size_t>(m_), 0.0));
    for (size_t r = 0; r < n; ++r) {
        for (int c = 0; c < m_; ++c) {
            bmat[r][static_cast<size_t>(c)] = Basis(c, order_, xs[r]);
        }
    }

    // 正规方程 A = B^T B + lambda * D2^T D2, rhs = B^T y
    std::vector<std::vector<double>> a(static_cast<size_t>(m_),
                                       std::vector<double>(static_cast<size_t>(m_), 0.0));
    std::vector<double> rhs(static_cast<size_t>(m_), 0.0);

    for (int i = 0; i < m_; ++i) {
        for (int j = 0; j < m_; ++j) {
            double s = 0.0;
            for (size_t r = 0; r < n; ++r) {
                s += bmat[r][static_cast<size_t>(i)] * bmat[r][static_cast<size_t>(j)];
            }
            a[static_cast<size_t>(i)][static_cast<size_t>(j)] = s;
        }
        double sy = 0.0;
        for (size_t r = 0; r < n; ++r) {
            sy += bmat[r][static_cast<size_t>(i)] * ys[r];
        }
        rhs[static_cast<size_t>(i)] = sy;
    }

    AddPenalty(a, lambda);

    if (!Solve(a, rhs, coef_)) return false;

    ready_ = true;
    return true;
}

bool SmoothClusterXY(std::vector<double>& zs, std::vector<double>& xs,
                     std::vector<double>& ys, double lambda) {
    if (zs.size() < 5 || zs.size() != xs.size() || zs.size() != ys.size()) return false;

    // 按 z 升序（调用方通常已排好，这里保证一下）
    std::vector<size_t> order(zs.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&zs](size_t a, size_t b) { return zs[a] < zs[b]; });

    std::vector<double> sz(zs.size()), sx(xs.size()), sy(ys.size());
    for (size_t i = 0; i < order.size(); ++i) {
        sz[i] = zs[order[i]];
        sx[i] = xs[order[i]];
        sy[i] = ys[order[i]];
    }

    SmoothingSpline spx, spy;
    if (!spx.Fit(sz, sx, lambda)) return false;
    if (!spy.Fit(sz, sy, lambda)) return false;

    for (size_t i = 0; i < sz.size(); ++i) {
        xs[i] = spx.Eval(sz[i]);
        ys[i] = spy.Eval(sz[i]);
    }
    zs.swap(sz);
    return true;
}

} // namespace rk4
