#ifndef SPYRALPP_FRIB_CONFIG_H
#define SPYRALPP_FRIB_CONFIG_H
#include "PeakFinder.h"
namespace rk4 {
struct FribConfig {
    double baseline_scale=100, delay=1100;
    int max_multiplicity=1;
    PeakFinderParameters peaks;
    FribConfig(){
        peaks.peak_separation=50;peaks.peak_prominence=20;
        peaks.peak_max_width=500;peaks.peak_threshold=100;
        peaks.rel_height=.85;peaks.peak_min_width=0;
    }
};
}
#endif
