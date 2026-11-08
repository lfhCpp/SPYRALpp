#ifndef RK4_RECO_SPLINE_H
#define RK4_RECO_SPLINE_H

//----------------------------------------------------------------------
// 平滑样条（P-spline）：三次 B 样条基 + 二阶差分惩罚
//
// 对应 Spyral 里的 make_smoothing_spline 用途：拟合前对簇内 x(z)、y(z) 降噪。
//
// 求解 (B^T B + lambda * D2^T D2) c = B^T y，其中 B 是 B 样条设计矩阵，
// D2 是系数上的二阶差分算子。lambda 越大越平滑；基函数个数固定（约 12+4 个），
// 所以正规方程只有十几阶，直接稠密高斯消元即可。
//----------------------------------------------------------------------

#include <vector>

namespace rk4 {

class SmoothingSpline {
public:
    // xs 需升序；lambda 为平滑强度（>=0）；n_knots 为内部节点数
    bool Fit(const std::vector<double>& xs, const std::vector<double>& ys,
             double lambda = 1.0, int n_knots = 12);

    double Eval(double x) const;

    // 在给定 xs 上求值（用于替换原始数据）
    std::vector<double> Smooth(const std::vector<double>& xs) const;

    bool Ready() const { return ready_; }

private:
    void BuildKnots(int n_knots);
    double Basis(int i, int order, double x) const;
    void AddPenalty(std::vector<std::vector<double>>& a, double lambda) const;
    bool Solve(std::vector<std::vector<double>>& a, std::vector<double>& b,
               std::vector<double>& out) const;

    bool ready_ = false;
    int order_ = 4;          // 三次（degree 3）
    double xmin_ = 0.0, xmax_ = 0.0;
    std::vector<double> knots_;
    std::vector<double> coef_;
    int m_ = 0;              // 基函数个数
};

// 对簇内坐标做平滑：把 (z -> x) 和 (z -> y) 分别拟合后替换。
// 输入点的 z 需升序。返回是否成功。
bool SmoothClusterXY(std::vector<double>& zs, std::vector<double>& xs,
                     std::vector<double>& ys, double lambda = 1.0);

} // namespace rk4

#endif // RK4_RECO_SPLINE_H
