#include "Response.h"

#include "../core/Units.h"

#include <algorithm>
#include <cmath>

namespace rk4 {

void GetResponse::Init(const DetectorConfig& cfg) {
    const int n_tb = cfg.num_tb > 0 ? cfg.num_tb : 512;
    response_.assign(static_cast<size_t>(n_tb), 0.0);

    const double c1 = 4095.0 * Units::E_CHARGE / cfg.amp_gain / 1.0e-15;
    const double denom = cfg.shaping_ns * cfg.clock_mhz * 1.0e-3;

    for (int i = 0; i < n_tb; ++i) {
        const double c2 = static_cast<double>(i) / denom;
        double value = c1 * std::exp(-3.0 * c2) * c2 * c2 * c2 * std::sin(c2);
        if (value < 0.0) value = 0.0;
        response_[static_cast<size_t>(i)] = value;
    }
}

void GetResponse::Apply(double electrons, double& amplitude, double& integral) const {
    double max_value = 0.0;
    double sum_value = 0.0;
    for (double r : response_) {
        double signal = r * electrons;
        if (signal > 4095.0) signal = 4095.0; // 11 bit 饱和
        if (signal > max_value) max_value = signal;
        sum_value += signal;
    }
    amplitude = max_value;
    integral = sum_value;
}

} // namespace rk4
