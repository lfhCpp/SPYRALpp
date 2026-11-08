#ifndef RK4_SIM_PADPLANE_H
#define RK4_SIM_PADPLANE_H

//----------------------------------------------------------------------
// pad 平面：几何查表
//
// 数据来自 padxy.csv（每个 pad 的中心，mm）与 pad_scale.csv（大小 pad 标记）。
// 空间查找用 2 mm 的桶网格 + 邻域最近中心搜索，避免依赖 attpc_engine 的
// pad_grid.npz（二进制、格式耦合）。
//----------------------------------------------------------------------

#include "BeamPads.h"

#include <string>
#include <unordered_set>
#include <vector>

namespace rk4 {

class PadPlane {
public:
    bool Load(const std::string& dir, const std::string& xy_file = "padxy.csv",
              const std::string& scale_file = "pad_scale.csv");

    // 返回最近的 pad 编号；超出灵敏区返回 -1
    int FindPad(double x_mm, double y_mm) const;

    // 可选：加载 pad 时间修正（Spyral 的 pad_time_correction.csv）
    bool LoadTimeOffsets(const std::string& dir,
                         const std::string& file = "pad_time_correction.csv");

    // pad 的时间偏移（单位 TB）；未加载时为 0
    double TimeOffset(int pad) const {
        return time_offset_.empty() ? 0.0 : time_offset_[pad];
    }

    bool IsBeamPad(int pad) const { return beam_pads_.count(pad) > 0; }

    int NumPads() const { return static_cast<int>(x_.size()); }
    double X(int pad) const { return x_[pad]; }
    double Y(int pad) const { return y_[pad]; }
    double Scale(int pad) const { return scale_[pad]; }

private:
    void BuildGrid();

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> scale_;
    std::vector<double> time_offset_;
    std::unordered_set<int> beam_pads_;

    // 桶网格
    double cell_size_ = 2.0; // mm
    double x_min_ = -300.0, y_min_ = -300.0;
    int nx_ = 0, ny_ = 0;
    std::vector<std::vector<int>> cells_;
};

} // namespace rk4

#endif // RK4_SIM_PADPLANE_H
