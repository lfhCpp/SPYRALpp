R__ADD_INCLUDE_PATH(src)
R__ADD_INCLUDE_PATH(extern/hdf5/include)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include <hdf5.h>
#include <stdexcept>
#include <cstdio>
void import_reference_clouds(const char* input,const char* output){
    hid_t file=H5Fopen(input,H5F_ACC_RDONLY,H5P_DEFAULT);
    if(file<0)throw std::runtime_error("Cannot open reference pointclouds");
    hid_t group=H5Gopen2(file,"cloud",H5P_DEFAULT);
    long long first=0,last=0;
    for(const auto& name:{"min_event","max_event"}){
        hid_t attr=H5Aopen(group,name,H5P_DEFAULT);
        H5Aread(attr,H5T_NATIVE_LLONG,std::string(name)=="min_event"?&first:&last);H5Aclose(attr);
    }
    std::vector<rk4::PointCloud> clouds;size_t points=0;
    for(long long event=first;event<=last;++event){
        char name[64];std::snprintf(name,sizeof(name),"cloud_%lld",event);
        if(H5Lexists(group,name,H5P_DEFAULT)<=0)continue;
        hid_t dataset=H5Dopen2(group,name,H5P_DEFAULT),space=H5Dget_space(dataset);
        hsize_t dims[2];H5Sget_simple_extent_dims(space,dims,nullptr);
        if(dims[1]!=8)throw std::runtime_error("Expected Nx8 pointcloud");
        std::vector<double> data(dims[0]*8);
        if(dims[0]&&H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data.data())<0)
            throw std::runtime_error("Cannot read pointcloud dataset");
        H5Sclose(space);H5Dclose(dataset);
        if(!dims[0])continue;
        rk4::PointCloud cloud;cloud.event=event;
        for(size_t i=0;i<dims[0];++i){
            const double* row=&data[8*i];rk4::CloudPoint p;
            p.x=row[0];p.y=row[1];p.z=row[2];p.amplitude=row[3];p.integral=row[4];
            p.pad=static_cast<int>(row[5]);p.tb=row[6];p.scale=row[7];cloud.points.push_back(p);
        }
        points+=cloud.points.size();clouds.push_back(std::move(cloud));
    }
    H5Gclose(group);H5Fclose(file);
    if(!rk4::WritePointClouds(output,clouds))throw std::runtime_error("Cannot write ROOT clouds");
    std::printf("Reference clouds imported: %zu events, %zu points\n",clouds.size(),points);
}
