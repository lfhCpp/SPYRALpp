#include "RecoConfig.h"
// Raw pointclouds -> native clusters -> PID estimates, with separate compute timers.
#include "io/CloudIO.h"
#include "io/Config.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
using namespace rk4;
using Clock=std::chrono::steady_clock;
static double Seconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }
int main(int argc,char** argv) {
    if(argc!=4&&argc!=5) { std::cerr<<"Usage: rk4_pid clouds.root estimates.csv config [--export-clusters]\n";return 1; }
    try {
    const bool export_clusters=argc==5;
    if(export_clusters&&std::string(argv[4])!="--export-clusters")throw std::runtime_error("Unknown PID argument");
    const auto start=Clock::now();
    Config cfg;if(!cfg.Load(argv[3]))throw std::runtime_error("Cannot load configuration");
    std::vector<PointCloud> clouds;
    const auto loading=Clock::now();
    if(!ReadPointClouds(argv[1],clouds))throw std::runtime_error("Cannot read pointclouds");
    const double load_seconds=Seconds(loading);
    const ClusterParams cluster_params=ReadClusterParams(cfg);

    Clusterer clusterer;clusterer.SetParams(cluster_params);
    Estimator estimator;ConfigureEstimator(estimator,cfg);
    const double field=cfg.GetDouble("detector.bfield_t",2.85);
    std::ofstream out(argv[2]);if(!out)throw std::runtime_error("Cannot open output");
    out<<"event,cluster,ok,brho,dedx,polar_deg\n"<<std::setprecision(17);
    std::ofstream cluster_file;
    if(export_clusters){
        cluster_file.open(std::string(argv[2])+".clusters.csv");
        if(!cluster_file)throw std::runtime_error("Cannot open cluster CSV");
        cluster_file<<"event,cluster,x,y,z,integral\n"<<std::setprecision(17);
    }
    double cluster_seconds=0,estimate_seconds=0;
    size_t total=0,accepted=0,processed=0,nonempty=0;std::set<long long> accepted_events;
    for(const auto& cloud:clouds) {
        if(!cloud.points.empty())++nonempty;
        std::vector<Cluster> clusters;
        const auto cluster_start=Clock::now();clusterer.Run(cloud,clusters);
        cluster_seconds+=Seconds(cluster_start);
        for(const auto& cluster:clusters) {
            if(export_clusters)for(const auto& point:cluster.points)
                cluster_file<<cloud.event<<','<<cluster.label<<','<<point.x<<','<<point.y<<','<<point.z<<','<<point.integral<<'\n';
            EstimateResult result;const auto estimate_start=Clock::now();
            const bool ok=estimator.Estimate(cluster,field,result);
            estimate_seconds+=Seconds(estimate_start);
            out<<cloud.event<<','<<cluster.label<<','<<(ok?1:0)<<','<<result.brho<<','<<result.dEdx<<','<<result.polar*180/std::acos(-1.0)<<'\n';
            ++total;if(ok){++accepted;accepted_events.insert(cloud.event);}
        }
        if(++processed%1000==0)std::cout<<"PROGRESS clouds "<<processed<<'/'<<clouds.size()<<std::endl;
    }
    if(export_clusters){cluster_file.close();if(!cluster_file)throw std::runtime_error("Incomplete cluster CSV");}
    out.close();if(!out)throw std::runtime_error("Incomplete CSV output");
    const double elapsed=Seconds(start);
    std::ofstream timing(std::string(argv[2])+".timing.json");
    timing<<std::setprecision(12)<<"{\n  \"stored_clouds\": "<<clouds.size()
          <<",\n  \"nonempty_clouds\": "<<nonempty
          <<",\n  \"clusters\": "<<total<<",\n  \"valid_estimates\": "<<accepted
          <<",\n  \"events_with_valid_estimates\": "<<accepted_events.size()
          <<",\n  \"cloud_load_seconds\": "<<load_seconds
          <<",\n  \"cluster_compute_seconds\": "<<cluster_seconds
          <<",\n  \"estimate_compute_seconds\": "<<estimate_seconds
          <<",\n  \"total_seconds\": "<<elapsed<<"\n}\n";
    if(!timing)throw std::runtime_error("Cannot write timing metadata");
    std::cout<<"PID: clouds="<<clouds.size()<<", clusters="<<total<<", valid="<<accepted<<", seconds="<<elapsed<<std::endl;
    return 0;
    } catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
