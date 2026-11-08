#ifndef RK4_RECO_GEOMETRY_H
#define RK4_RECO_GEOMETRY_H

//----------------------------------------------------------------------
// 几何拟合工具
//----------------------------------------------------------------------

#include <vector>

namespace rk4 {

// 最小二乘圆拟合（Kåsa 代数解）
// 输入 N 个点，输出圆心与半径；点数不足或退化时返回 false
bool FitCircle(const std::vector<double>& xs, const std::vector<double>& ys,
               double& cx, double& cy, double& radius);

// Spyral reference geometry: centered fit and mean radial distance.
// Reject fewer than three points: a two-point circle is underdetermined.
bool FitCircleSpyral(const std::vector<double>& xs,const std::vector<double>& ys,
                     double& cx,double& cy,double& radius);

} // namespace rk4

#endif // RK4_RECO_GEOMETRY_H
