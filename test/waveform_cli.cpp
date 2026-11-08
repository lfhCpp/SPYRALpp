#include "reco/WaveProcessor.h"
#include "reco/PeakFinder.h"
#include <iostream>
#include <iomanip>
int main() {
    int n, integer_input;
    double scale;
    PeakFinderParameters p;
    std::cout << std::setprecision(17);
    while (std::cin >> n >> scale >> integer_input >> p.peak_separation >> p.peak_prominence
                    >> p.peak_max_width >> p.peak_threshold >> p.rel_height) {
        std::vector<double> raw(n), filtered(n);
        for (double& x : raw) std::cin >> x;
        auto wave = WaveProcessor::createUnique(n, scale);
        if (!wave) return 2;
        wave->setIntegerInput(integer_input);
        if (!wave->process(raw.data(), filtered.data())) return 3;
        for (double x : filtered) std::cout << x << ' ';
        std::cout << '\n';
        if (integer_input) for (double& x : filtered) x = std::trunc(x);
        PeakFinder finder(p);
        finder.setTrace(filtered);
        auto peaks = finder.findPeaksWithDither([] { return 0.0; });
        std::cout << finder.DitherDraws() << ' ' << peaks.size() << '\n';
        for (const auto& q : peaks)
            std::cout << q.centroid << ' ' << q.amplitude << ' ' << q.integral << ' '
                      << q.left_ip << ' ' << q.right_ip << ' ' << q.prominence << ' '
                      << q.width << ' ' << q.dither_index << '\n';
    }
}
