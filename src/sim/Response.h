#ifndef RK4_SIM_RESPONSE_H
#define RK4_SIM_RESPONSE_H

//----------------------------------------------------------------------
// GET 电子学响应
//
// 响应函数取自 attpc_engine（原文献 NIM A 论文）：
//   resp(t) = c1 * exp(-3 c2) * c2^3 * sin(c2)
// 一个点上沉积的电子数乘响应，取最大值作幅度、求和作积分。
//----------------------------------------------------------------------

#include "DetectorConfig.h"

#include <vector>

namespace rk4 {

class GetResponse {
public:
    void Init(const DetectorConfig& cfg);

    // electrons 为该点沉积的电子数（已含增益）
    void Apply(double electrons, double& amplitude, double& integral) const;

    int Size() const { return static_cast<int>(response_.size()); }

private:
    std::vector<double> response_;
};

} // namespace rk4

#endif // RK4_SIM_RESPONSE_H
