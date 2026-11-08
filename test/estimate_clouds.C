// PID-only replay of the default local cluster/estimate configuration.
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
void estimate_clouds(const char* input,const char* output) {
    std::vector<rk4::PointCloud> clouds;
    if(!rk4::ReadPointClouds(input,clouds)) throw std::runtime_error("Cannot read clouds");
    rk4::Clusterer clusterer;clusterer.SetParams(rk4::ClusterParams{});
    rk4::Estimator estimator;estimator.SetUseReinsch(true,100);
    estimator.SetSmallPadRadius(152);estimator.SetBeamRegionRadius(25);
    std::ofstream out(output);
    if(!out) throw std::runtime_error("Cannot open PID CSV");
    out<<"event,cluster,ok,brho,dedx,polar_deg\n"<<std::setprecision(14);
    size_t total=0,accepted=0;
    for(const auto& cloud:clouds) {
        std::vector<rk4::Cluster> clusters;clusterer.Run(cloud,clusters);
        for(const auto& cluster:clusters) {
            rk4::EstimateResult result;bool ok=estimator.Estimate(cluster,2.85,result);
            out<<cloud.event<<','<<cluster.label<<','<<(ok?1:0)<<','<<result.brho<<','<<result.dEdx<<','<<result.polar*180/std::acos(-1.0)<<'\n';
            ++total;if(ok)++accepted;
        }
    }
    out.close();if(!out) throw std::runtime_error("Incomplete PID CSV");
    std::printf("PID replay: clouds=%zu, clusters=%zu, valid=%zu\n",clouds.size(),total,accepted);
}
