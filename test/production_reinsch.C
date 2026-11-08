R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/Estimator.h"
#include <map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
void production_reinsch(const char* input="../results/our_clusters5k.csv", const char* output="../results/production_reinsch5k.csv") {
    std::ifstream in(input);
    if (!in) throw std::runtime_error("Cannot open clusters");
    std::map<std::pair<int,int>,rk4::Cluster> groups;
    std::string line; std::getline(in,line);
    while (std::getline(in,line)) {
        std::replace(line.begin(),line.end(),',',' ');
        std::istringstream row(line); int event,label; rk4::CloudPoint p;
        if (!(row>>event>>label>>p.x>>p.y>>p.z>>p.integral)) throw std::runtime_error("Malformed cluster CSV");
        groups[{event,label}].points.push_back(p);
    }
    std::ofstream out(output);
    if (!out) throw std::runtime_error("Cannot open output");
    out<<"event,cluster,brho,dedx,polar_deg\n"<<std::setprecision(12);
    rk4::Estimator estimator; estimator.SetUseReinsch(true);
    size_t accepted=0;
    for (auto& entry:groups) {
        rk4::EstimateResult r;
        if (!estimator.Estimate(entry.second,2.85,r)) continue;
        if (!r.ok || !std::isfinite(r.vertex_z) || r.radius<=0 || r.arclength<=0)
            throw std::runtime_error("Invalid accepted estimate");
        out<<entry.first.first<<','<<entry.first.second<<','<<r.brho<<','<<r.dEdx<<','<<r.polar*180/M_PI<<'\n';
        ++accepted;
    }
    std::printf("Production Reinsch: accepted=%zu / %zu\n",accepted,groups.size());
}
