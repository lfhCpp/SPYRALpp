R__ADD_INCLUDE_PATH(src)
R__ADD_INCLUDE_PATH(extern/hdf5/include)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include <hdf5.h>
#include <fstream>
#include <map>
#include <array>
#include <stdexcept>
#include <iomanip>
void compare_raw_clouds(const char* reference, const char* native, const char* output, int limit=12432) {
    std::vector<rk4::PointCloud> clouds;
    if(!rk4::ReadPointClouds(native,clouds))throw std::runtime_error("Cannot read native clouds");
    std::map<long long,const rk4::PointCloud*> by_event;
    for(const auto& c:clouds)by_event[c.event]=&c;
    hid_t file=H5Fopen(reference,H5F_ACC_RDONLY,H5P_DEFAULT),group=H5Gopen2(file,"cloud",H5P_DEFAULT);
    size_t events=0,points=0,missing=0,extra=0,shape=0,differing=0;
    std::array<double,8> maxerr{};
    std::vector<int> bad;
    for(int e=0;e<limit;++e) {
        const std::string key="cloud_"+std::to_string(e);
        const bool exists=H5Lexists(group,key.c_str(),H5P_DEFAULT)>0;
        const auto it=by_event.find(e);
        if(!exists){if(it!=by_event.end())++extra;continue;}
        ++events;
        if(it==by_event.end()){++missing;continue;}
        hid_t ds=H5Dopen2(group,key.c_str(),H5P_DEFAULT),sp=H5Dget_space(ds);
        hsize_t dims[2];H5Sget_simple_extent_dims(sp,dims,nullptr);
        std::vector<double> data(dims[0]*8);
        if(dims[0])H5Dread(ds,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data.data());
        H5Sclose(sp);H5Dclose(ds);points+=dims[0];
        const auto& cloud=*it->second;
        if(cloud.points.size()!=dims[0]){++shape;if(bad.size()<20)bad.push_back(e);continue;}
        bool diff=false;
        for(size_t i=0;i<dims[0];++i){
            const auto& p=cloud.points[i];
            const double values[8]={p.x,p.y,p.z,p.amplitude,p.integral,double(p.pad),p.tb,p.scale};
            for(size_t j=0;j<8;++j){
                double error=std::abs(values[j]-data[i*8+j]);maxerr[j]=std::max(maxerr[j],error);
                if(error>1e-9)diff=true;
            }
        }
        if(diff){++differing;if(bad.size()<20)bad.push_back(e);}
    }
    H5Gclose(group);H5Fclose(file);
    std::ofstream out(output);out<<std::setprecision(17);
    out<<"{\n  \"reference_events\": "<<events<<", \"reference_points\": "<<points
       <<", \"native_events\": "<<clouds.size()<<",\n  \"missing\": "<<missing<<", \"extra\": "<<extra
       <<", \"shape_differences\": "<<shape<<", \"value_differences\": "<<differing<<",\n  \"max_absolute_errors_xyz_amp_integral_pad_tb_scale\": [";
    for(int i=0;i<8;++i)out<<(i?", ":"")<<maxerr[i];out<<"],\n  \"first_differing_events\": [";
    for(size_t i=0;i<bad.size();++i)out<<(i?", ":"")<<bad[i];out<<"]\n}\n";
    std::printf("Compared %zu events %zu points: missing=%zu extra=%zu shapes=%zu values=%zu\n",events,points,missing,extra,shape,differing);
}
