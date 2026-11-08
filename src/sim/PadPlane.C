#include "PadPlane.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>

namespace rk4 {

namespace {
bool ReadCsvColumn(const std::string& path, std::vector<double>& out, int column = 0) {
    std::ifstream file(path.c_str());
    if (!file.is_open()) {
        std::fprintf(stderr, "[PadPlane] 无法打开 %s\n", path.c_str());
        return false;
    }

    // 整个文件读进来，按 \n 或 \r 切分。
    // Spyral 的 pad_time_correction.csv 用的是老式 CR 行尾，只按 \n 会读成一行。
    std::string text((std::istreambuf_iterator<char>(file)),
                     std::istreambuf_iterator<char>());

    out.clear();
    size_t pos = 0;
    bool first_line = true;

    while (pos < text.size()) {
        size_t end = text.find_first_of("\r\n", pos);
        if (end == std::string::npos) end = text.size();

        std::string line = text.substr(pos, end - pos);
        pos = end + 1;

        if (line.empty()) continue;
        if (first_line) { // 跳过表头
            first_line = false;
            continue;
        }

        // 取第 column 个逗号分隔字段
        size_t field_start = 0;
        int col = 0;
        while (col < column) {
            const size_t comma = line.find(',', field_start);
            if (comma == std::string::npos) break;
            field_start = comma + 1;
            ++col;
        }
        if (col == column) {
            out.push_back(std::atof(line.c_str() + field_start));
        }
    }

    return !out.empty();
}

} // namespace

bool PadPlane::Load(const std::string& dir, const std::string& xy_file,
                    const std::string& scale_file) {
    std::vector<double> xs, ys;
    if (!ReadCsvColumn(dir + "/" + xy_file, xs, 0)) return false;
    if (!ReadCsvColumn(dir + "/" + xy_file, ys, 1)) return false;
    if (xs.size() != ys.size()) {
        std::fprintf(stderr, "[PadPlane] x/y 数量不一致\n");
        return false;
    }

    std::vector<double> scales;
    if (!ReadCsvColumn(dir + "/" + scale_file, scales, 0)) return false;
    if (scales.size() != xs.size()) {
        std::fprintf(stderr, "[PadPlane] pad 数与 scale 数不一致 (%zu vs %zu)\n",
                     xs.size(), scales.size());
        return false;
    }

    x_ = xs;
    y_ = ys;
    scale_ = scales;

    for (int i = 0; i < kNumBeamPads; ++i) {
        beam_pads_.insert(kBeamPads[i]);
    }

    BuildGrid();
    return true;
}

bool PadPlane::LoadTimeOffsets(const std::string& dir, const std::string& file) {
    std::vector<double> offsets;
    if (!ReadCsvColumn(dir + "/" + file, offsets, 0)) return false;
    if (offsets.size() != x_.size()) {
        std::fprintf(stderr, "[PadPlane] 时间修正数与 pad 数不一致 (%zu vs %zu)\n",
                     offsets.size(), x_.size());
        return false;
    }
    time_offset_ = offsets;
    return true;
}

void PadPlane::BuildGrid() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    double xmin = nan, xmax = nan, ymin = nan, ymax = nan;
    for (size_t i = 0; i < x_.size(); ++i) {
        xmin = std::isnan(xmin) ? x_[i] : std::min(xmin, x_[i]);
        xmax = std::isnan(xmax) ? x_[i] : std::max(xmax, x_[i]);
        ymin = std::isnan(ymin) ? y_[i] : std::min(ymin, y_[i]);
        ymax = std::isnan(ymax) ? y_[i] : std::max(ymax, y_[i]);
    }
    x_min_ = xmin;
    y_min_ = ymin;
    nx_ = static_cast<int>((xmax - xmin) / cell_size_) + 1;
    ny_ = static_cast<int>((ymax - ymin) / cell_size_) + 1;

    cells_.assign(static_cast<size_t>(nx_ * ny_), std::vector<int>());
    for (size_t i = 0; i < x_.size(); ++i) {
        const int ix = static_cast<int>((x_[i] - x_min_) / cell_size_);
        const int iy = static_cast<int>((y_[i] - y_min_) / cell_size_);
        if (ix >= 0 && ix < nx_ && iy >= 0 && iy < ny_) {
            cells_[static_cast<size_t>(iy * nx_ + ix)].push_back(static_cast<int>(i));
        }
    }
}

int PadPlane::FindPad(double x_mm, double y_mm) const {
    const int ix = static_cast<int>((x_mm - x_min_) / cell_size_);
    const int iy = static_cast<int>((y_mm - y_min_) / cell_size_);

    int best = -1;
    double best_d2 = 9.0; // 3 mm 以内才算命中

    // 从 0 环开始向外扩，找到第一个有候选的环后最多再看一环
    for (int r = 0; r <= 4; ++r) {
        bool found = false;
        for (int dy = -r; dy <= r; ++dy) {
            for (int dx = -r; dx <= r; ++dx) {
                if (r > 0 && std::abs(dx) != r && std::abs(dy) != r) continue;
                const int cx = ix + dx;
                const int cy = iy + dy;
                if (cx < 0 || cx >= nx_ || cy < 0 || cy >= ny_) continue;
                for (int pad : cells_[static_cast<size_t>(cy * nx_ + cx)]) {
                    const double ddx = x_[pad] - x_mm;
                    const double ddy = y_[pad] - y_mm;
                    const double d2 = ddx * ddx + ddy * ddy;
                    if (d2 < best_d2) {
                        best_d2 = d2;
                        best = pad;
                        found = true;
                    }
                }
            }
        }
        if (found && r > 1) break;
    }

    return best;
}

} // namespace rk4
