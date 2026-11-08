#include "app/RecoConfig.h"
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <unistd.h>
int main(){
    try{
        rk4::Config empty;auto defaults=rk4::ReadClusterParams(empty);
        if(!defaults.use_hdbscan||defaults.hdbscan_allow_single||defaults.min_cloud_size!=50||defaults.min_cluster_size!=5)
            throw std::runtime_error("Production defaults differ from the reference profile");
        char name[]="/tmp/spyralpp-config-XXXXXX";int descriptor=mkstemp(name);
        if(descriptor<0)throw std::runtime_error("Cannot create temporary config");close(descriptor);
        {std::ofstream f(name);f<<"[cluster]\nmin_cloud_size = 50\nmin_cluster_size = 10\nuse_hdbscan = false # explicit DBSCAN\nhdbscan_allow_single = true ; explicit option\n[estimate]\nsmoothing_factor = 300\nmethod = \"reinsch\"\n";}
        rk4::Config cfg;if(!cfg.Load(name))throw std::runtime_error("Load failed");auto p=rk4::ReadClusterParams(cfg);
        if(p.use_hdbscan||!p.hdbscan_allow_single||p.min_cloud_size!=50||p.min_cluster_size!=10)
            throw std::runtime_error("Cloud threshold was aliased to cluster size, or boolean parsing failed");
        if(cfg.GetDouble("estimate.reinsch_lambda")!=300||cfg.GetString("estimate.method")!="reinsch")
            throw std::runtime_error("Smoothing alias or quoted token failed");
        {std::ofstream f(name);f<<"Outer:\n  Inner:\n    value: 7\n    label: \"literal # text\"\n";}
        if(!cfg.Load(name)||cfg.GetInt("outer.inner.value")!=7||cfg.GetString("outer.inner.label")!="literal # text")
            throw std::runtime_error("Nested YAML or quoted comment failed");
        unlink(name);std::cout<<"PASS: HDBSCAN defaults, distinct thresholds, booleans, lambda alias and nested YAML\n";return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
