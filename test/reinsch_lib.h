#ifndef RK4_REINSCH_LIB_H
#define RK4_REINSCH_LIB_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

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
                for (int j = 0; j < n; ++j)   // j 遍历全部列（关键修复）
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
    if (n <= 20) {   // UNIT_A_DUMP
        FILE* fa = std::fopen("/tmp/unit_A.csv", "w");
        std::fprintf(fa, "i,j,a\n");
        for (int k = 0; k < A.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
                std::fprintf(fa, "%d,%d,%.12g\n", (int)it.row(), (int)it.col(), it.value());
        std::fclose(fa);
    }
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
