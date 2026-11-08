#include "PeakFinder.h"
#include "NumpyArgsort.h"

// SciPy find_peaks order: local maxima, distance, prominence, continuous width.
// Hardware threshold is applied by GetTrace only after these selections.
std::vector<int> PeakFinder::findLocalMaxima(const std::vector<double>& data) {
    std::vector<int> maxima;
    const int n = static_cast<int>(data.size());
    for (int i = 1; i < n - 1; ++i) {
        if (data[i-1] >= data[i]) continue;
        int right = i + 1;
        while (right < n - 1 && data[right] == data[i]) ++right;
        if (data[right] < data[i]) maxima.push_back((i + right - 1) / 2);
        i = right - 1;
    }
    return maxima;
}

double PeakFinder::calculateProminence(const std::vector<double>& data, int p,
                                       int& left_base, int& right_base) {
    double left = data[p], right = data[p];
    left_base = right_base = p;
    for (int i = p; i >= 0 && data[i] <= data[p]; --i) {
        if (data[i] < left) { left = data[i]; left_base = i; }
    }
    for (int i = p; i < static_cast<int>(data.size()) && data[i] <= data[p]; ++i) {
        if (data[i] < right) { right = data[i]; right_base = i; }
    }
    return data[p] - std::max(left, right);
}

void PeakFinder::calculateWidth(const std::vector<double>& data, int p, double prominence,
                                int left_base, int right_base, double& left, double& right) {
    const double height = data[p] - params_.rel_height * prominence;
    int i = p;
    while (left_base < i && height < data[i]) --i;
    left = i;
    if (data[i] < height) left += (height - data[i]) / (data[i+1] - data[i]);
    i = p;
    while (i < right_base && height < data[i]) ++i;
    right = i;
    if (data[i] < height) right -= (height - data[i]) / (data[i-1] - data[i]);
}

std::vector<Peak> PeakFinder::findPeaks(std::mt19937& rng) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    return findPeaksWithDither([&] { return uniform(rng); });
}

std::vector<Peak> PeakFinder::findPeaksWithDither(const std::function<double()>& draw) {
    peaks_.clear();
    dither_draws_ = 0;
    if (trace_.size() < 3) return peaks_;
    if (!std::isfinite(params_.peak_separation) || params_.peak_separation < 1 || params_.rel_height < 0)
        throw std::invalid_argument("Peak distance must be >= 1 and rel_height >= 0");
    const auto candidates = findLocalMaxima(trace_);
    std::vector<double> heights;
    heights.reserve(candidates.size());
    for (int p : candidates) heights.push_back(trace_[p]);
    auto order = rk4::NumpyArgsort(std::move(heights));
    std::vector<bool> keep(candidates.size(), true);
    // SciPy walks ascending argsort in reverse; ties are observable.
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const size_t k = *it;
        if (!keep[k]) continue;
        for (size_t j = k; j > 0;) {
            --j;
            if (candidates[k] - candidates[j] >= params_.peak_separation) break;
            keep[j] = false;
        }
        for (size_t j = k + 1; j < candidates.size(); ++j) {
            if (candidates[j] - candidates[k] >= params_.peak_separation) break;
            keep[j] = false;
        }
    }
    for (size_t k = 0; k < candidates.size(); ++k) {
        if (!keep[k]) continue;
        const int p = candidates[k];
        int lb, rb;
        const double prominence = calculateProminence(trace_, p, lb, rb);
        if (prominence < params_.peak_prominence) continue;
        double left, right;
        calculateWidth(trace_, p, prominence, lb, rb, left, right);
        const double width = right - left;
        if (width < params_.peak_min_width || width > params_.peak_max_width) continue;
        Peak peak;
        peak.centroid = p + draw();
        peak.dither_index = dither_draws_++;
        peak.amplitude = trace_[p];
        peak.prominence = prominence;
        peak.width = width;
        peak.left_ip = static_cast<int>(std::floor(left));
        peak.right_ip = static_cast<int>(std::ceil(right));
        for (int i = peak.left_ip; i < peak.right_ip; ++i)
            peak.integral += std::abs(trace_[i]);
        if (peak.amplitude > params_.peak_threshold) peaks_.push_back(peak);
    }
    return peaks_;
}
