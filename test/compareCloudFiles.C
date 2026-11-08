R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include "TSystem.h"
#include <cstdio>
void compareCloudFiles(const char* first="../results/local_raw5000_1t.root",const char* second="../results/local_raw5000_8t.root") {
    std::vector<rk4::PointCloud>a,b;
    bool ok=rk4::ReadPointClouds(first,a)&&rk4::ReadPointClouds(second,b)&&a.size()==b.size();
    size_t points=0;
    for(size_t i=0;ok&&i<a.size();++i){
        ok=a[i].event==b[i].event&&a[i].points.size()==b[i].points.size();
        for(size_t j=0;ok&&j<a[i].points.size();++j){
            const auto&x=a[i].points[j];const auto&y=b[i].points[j];
            ok=x.pad==y.pad&&x.tb==y.tb&&x.x==y.x&&x.y==y.y&&x.z==y.z&&
                x.amplitude==y.amplitude&&x.integral==y.integral&&x.scale==y.scale;
            ++points;
        }
    }
    std::printf("Pointcloud equality: %s, events=%zu, points=%zu\n",ok?"PASS":"FAIL",a.size(),points);
    gSystem->Exit(ok?0:1);
}
