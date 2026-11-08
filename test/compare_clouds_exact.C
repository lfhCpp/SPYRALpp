R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include <fstream>
#include <stdexcept>
void compare_clouds_exact(const char* first,const char* second,const char* output) {
    std::vector<rk4::PointCloud> a,b;
    if(!rk4::ReadPointClouds(first,a)||!rk4::ReadPointClouds(second,b))throw std::runtime_error("Cannot read clouds");
    size_t differences=a.size()!=b.size(),points=0;
    for(size_t i=0;i<std::min(a.size(),b.size());++i){
        if(a[i].event!=b[i].event||a[i].points.size()!=b[i].points.size()){++differences;continue;}
        for(size_t j=0;j<a[i].points.size();++j){
            const auto& x=a[i].points[j];const auto& y=b[i].points[j];++points;
            if(x.pad!=y.pad||x.tb!=y.tb||x.x!=y.x||x.y!=y.y||x.z!=y.z||x.amplitude!=y.amplitude||x.integral!=y.integral||x.scale!=y.scale)++differences;
        }
    }
    std::ofstream out(output);
    out<<"{\"clouds\": "<<a.size()<<", \"points\": "<<points<<", \"exact_differences\": "<<differences<<"}\n";
    if(differences)throw std::runtime_error("Clouds differ");
}
