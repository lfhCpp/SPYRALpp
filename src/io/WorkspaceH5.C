#include "WorkspaceIO.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <limits>
#ifdef RK4_WITH_HDF5
#include <hdf5.h>
namespace rk4 {
namespace {
struct H5 {
    hid_t id=-1;herr_t(*close)(hid_t)=nullptr;
    H5(hid_t handle,herr_t(*closer)(hid_t)):id(handle),close(closer){if(id<0)throw std::runtime_error("HDF5 stage operation failed");}
    ~H5(){if(id>=0)close(id);}
    H5(const H5&)=delete;H5&operator=(const H5&)=delete;
    operator hid_t()const{return id;}
};
void Require(bool good,const std::string& message){if(!good)throw std::runtime_error(message);}
bool HasAttr(hid_t object,const char* name){return H5Aexists(object,name)>0;}
double RealAttr(hid_t object,const char* name,double fallback){
    if(!HasAttr(object,name))return fallback;H5 attr(H5Aopen(object,name,H5P_DEFAULT),H5Aclose);double value;
    Require(H5Aread(attr,H5T_NATIVE_DOUBLE,&value)>=0,"Cannot read HDF5 attribute "+std::string(name));return value;
}
long long IntAttr(hid_t object,const char* name,long long fallback){
    if(!HasAttr(object,name))return fallback;H5 attr(H5Aopen(object,name,H5P_DEFAULT),H5Aclose);long long value;
    Require(H5Aread(attr,H5T_NATIVE_LLONG,&value)>=0,"Cannot read HDF5 attribute "+std::string(name));return value;
}
void PutInt(hid_t object,const char* name,long long value){
    H5 space(H5Screate(H5S_SCALAR),H5Sclose);H5 attr(H5Acreate2(object,name,H5T_STD_I64LE,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    Require(H5Awrite(attr,H5T_NATIVE_LLONG,&value)>=0,"Cannot write HDF5 integer attribute");
}
void PutReal(hid_t object,const char* name,double value){
    H5 space(H5Screate(H5S_SCALAR),H5Sclose);H5 attr(H5Acreate2(object,name,H5T_IEEE_F64LE,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    Require(H5Awrite(attr,H5T_NATIVE_DOUBLE,&value)>=0,"Cannot write HDF5 real attribute");
}
void PutText(hid_t object,const char* name,const std::string& value){
    H5 type(H5Tcopy(H5T_C_S1),H5Tclose);H5Tset_size(type,value.size()+1);H5Tset_strpad(type,H5T_STR_NULLTERM);
    H5 space(H5Screate(H5S_SCALAR),H5Sclose);H5 attr(H5Acreate2(object,name,type,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    Require(H5Awrite(attr,type,value.c_str())>=0,"Cannot write HDF5 text attribute");
}
std::string TextAttr(hid_t object,const char* name,const std::string& fallback){
    if(!HasAttr(object,name))return fallback;H5 attr(H5Aopen(object,name,H5P_DEFAULT),H5Aclose);H5 type(H5Aget_type(attr),H5Tclose);
    if(H5Tis_variable_str(type)){char* text=nullptr;Require(H5Aread(attr,type,&text)>=0,"Cannot read HDF5 text");std::string value=text?text:"";if(text)H5free_memory(text);return value;}
    std::vector<char> text(H5Tget_size(type)+1,0);Require(H5Aread(attr,type,text.data())>=0,"Cannot read HDF5 text");return std::string(text.data());
}
EventMetadata Meta(hid_t object,long long event){
    EventMetadata m;m.orig_run=IntAttr(object,"orig_run",-1);m.orig_event=IntAttr(object,"orig_event",event);
    m.ic_amplitude=RealAttr(object,"ic_amplitude",-1);m.ic_centroid=RealAttr(object,"ic_centroid",-1);
    m.ic_integral=RealAttr(object,"ic_integral",-1);m.ic_multiplicity=RealAttr(object,"ic_multiplicity",-1);return m;
}
void PutMeta(hid_t object,const EventMetadata& m,long long event){
    PutInt(object,"orig_run",m.orig_run);PutInt(object,"orig_event",m.orig_event<0?event:m.orig_event);
    PutReal(object,"ic_amplitude",m.ic_amplitude);
    if(m.ic_multiplicity>=0){PutInt(object,"ic_centroid",static_cast<long long>(m.ic_centroid));PutInt(object,"ic_integral",static_cast<long long>(m.ic_integral));PutInt(object,"ic_multiplicity",static_cast<long long>(m.ic_multiplicity));}
    else{PutReal(object,"ic_centroid",m.ic_centroid);PutReal(object,"ic_integral",m.ic_integral);PutReal(object,"ic_multiplicity",m.ic_multiplicity);}
}
StageFileInfo Info(hid_t group){
    StageFileInfo info;info.first_event=IntAttr(group,"min_event",0);info.last_event=IntAttr(group,"max_event",-1);
    info.run=IntAttr(group,"run_number",-1);info.raw_events=IntAttr(group,"raw_events",0);
    info.excluded_ic=IntAttr(group,"excluded_ic",0);info.good_downscale=IntAttr(group,"good_ic_downscale",0);
    if(HasAttr(group,"dither_seed")){H5 a(H5Aopen(group,"dither_seed",H5P_DEFAULT),H5Aclose);Require(H5Aread(a,H5T_NATIVE_UINT64,&info.seed)>=0,"Cannot read random seed");}
    info.waveform_method=TextAttr(group,"waveform_method","baseline");return info;
}
void PutInfo(hid_t group,const StageFileInfo& info){
    PutInt(group,"min_event",info.first_event);PutInt(group,"max_event",info.last_event);PutInt(group,"run_number",info.run);
    PutInt(group,"raw_events",info.raw_events);PutInt(group,"excluded_ic",info.excluded_ic);PutInt(group,"good_ic_downscale",info.good_downscale);
    H5 space(H5Screate(H5S_SCALAR),H5Sclose);H5 seed(H5Acreate2(group,"dither_seed",H5T_STD_U64LE,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    Require(H5Awrite(seed,H5T_NATIVE_UINT64,&info.seed)>=0,"Cannot write random seed");
    PutText(group,"waveform_method",info.waveform_method);PutText(group,"spyralpp_schema","workspace-1");
}
std::vector<std::pair<long long,std::string>> Children(hid_t group,const std::string& prefix){
    std::vector<std::string> names;
    auto callback=[](hid_t,const char* name,const H5L_info2_t*,void* data)->herr_t{static_cast<std::vector<std::string>*>(data)->push_back(name);return 0;};
    Require(H5Literate(group,H5_INDEX_NAME,H5_ITER_INC,nullptr,callback,&names)>=0,"Cannot list stage events");
    std::vector<std::pair<long long,std::string>> events;
    for(const auto& name:names)if(name.compare(0,prefix.size(),prefix)==0){
        const auto suffix=name.substr(prefix.size());if(suffix.empty()||suffix.find_first_not_of("0123456789")!=std::string::npos)continue;
        const auto event=std::stoll(suffix);events.emplace_back(event,name);
    }
    std::sort(events.begin(),events.end());return events;
}
std::vector<double> Matrix(hid_t dataset,size_t columns,size_t& rows){
    H5 space(H5Dget_space(dataset),H5Sclose);Require(H5Sget_simple_extent_ndims(space)==2,"Stage cloud must be a matrix");
    hsize_t shape[2]={0,0};H5Sget_simple_extent_dims(space,shape,nullptr);Require(shape[1]==columns,"Unexpected stage cloud column count");
    Require(shape[0]<=std::numeric_limits<size_t>::max()/columns,"Stage cloud size overflows");rows=static_cast<size_t>(shape[0]);std::vector<double> data(rows*columns);
    if(rows)Require(H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data.data())>=0,"Cannot read stage cloud");
    Require(std::all_of(data.begin(),data.end(),[](double x){return std::isfinite(x);}),"Nonfinite stage cloud value");return data;
}
void PutMatrix(hid_t parent,const std::string& name,const std::vector<double>& data,size_t columns){
    Require(data.size()%columns==0,"Invalid output cloud shape");hsize_t shape[2]={data.size()/columns,columns};H5 space(H5Screate_simple(2,shape,nullptr),H5Sclose);
    H5 dataset(H5Dcreate2(parent,name.c_str(),H5T_IEEE_F64LE,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Dclose);
    if(!data.empty())Require(H5Dwrite(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data.data())>=0,"Cannot write stage cloud");
}
void ExpandRange(StageFileInfo& info,const std::vector<std::pair<long long,std::string>>& events){
    if(events.empty())return;
    info.first_event=std::min(info.first_event,events.front().first);info.last_event=std::max(info.last_event,events.back().first);
}
}
PointcloudData ReadPointcloudH5(const std::string& path){
    H5 file(H5Fopen(path.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);
    Require(H5Lexists(file,"cloud",H5P_DEFAULT)>0,"Pointcloud file lacks /cloud: "+path);H5 group(H5Gopen2(file,"cloud",H5P_DEFAULT),H5Gclose);
    PointcloudData output;output.info=Info(group);const auto names=Children(group,"cloud_");ExpandRange(output.info,names);output.events.reserve(names.size());
    for(const auto& entry:names){
        H5 dataset(H5Dopen2(group,entry.second.c_str(),H5P_DEFAULT),H5Dclose);size_t rows=0;const auto data=Matrix(dataset,8,rows);
        PointCloud cloud;cloud.event=entry.first;cloud.metadata=Meta(dataset,cloud.event);cloud.points.reserve(rows);
        for(size_t i=0;i<rows;++i){const double* r=&data[i*8];CloudPoint point;
            Require(r[5]>=std::numeric_limits<int>::min()&&r[5]<=std::numeric_limits<int>::max()&&r[5]==std::trunc(r[5]),"Invalid pad identifier");
            point.x=r[0];point.y=r[1];point.z=r[2];point.amplitude=r[3];point.integral=r[4];point.pad=static_cast<int>(r[5]);point.tb=r[6];point.scale=r[7];cloud.points.push_back(point);
        }
        output.events.push_back(std::move(cloud));
    }
    return output;
}
void WritePointcloudH5(const std::string& path,const PointcloudData& input){
    H5 file(H5Fcreate(path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);H5 group(H5Gcreate2(file,"cloud",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);PutInfo(group,input.info);
    for(const auto& cloud:input.events){
        Require(!cloud.pending_dither,"Cannot write unresolved dithering");std::vector<double> data;data.reserve(cloud.points.size()*8);
        for(const auto& p:cloud.points)data.insert(data.end(),{p.x,p.y,p.z,p.amplitude,p.integral,static_cast<double>(p.pad),p.tb,p.scale});
        const auto name="cloud_"+std::to_string(cloud.event);PutMatrix(group,name,data,8);H5 dataset(H5Dopen2(group,name.c_str(),H5P_DEFAULT),H5Dclose);PutMeta(dataset,cloud.metadata,cloud.event);
    }
}
ClusterData ReadClusterH5(const std::string& path){
    H5 file(H5Fopen(path.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);Require(H5Lexists(file,"cluster",H5P_DEFAULT)>0,"Cluster file lacks /cluster: "+path);
    H5 group(H5Gopen2(file,"cluster",H5P_DEFAULT),H5Gclose);ClusterData output;output.info=Info(group);
    const auto names=Children(group,"event_");ExpandRange(output.info,names);output.events.reserve(names.size());
    for(const auto& entry:names){
        H5 event(H5Gopen2(group,entry.second.c_str(),H5P_DEFAULT),H5Gclose);ClusterEvent item;item.event=entry.first;item.metadata=Meta(event,item.event);
        const long long count=IntAttr(event,"nclusters",-1);Require(count>=0&&count<=std::numeric_limits<int>::max(),"Invalid nclusters");item.clusters.reserve(static_cast<size_t>(count));
        for(int index=0;index<count;++index){
            const auto name="cluster_"+std::to_string(index);Require(H5Lexists(event,name.c_str(),H5P_DEFAULT)>0,"Missing declared cluster "+name);
            H5 cluster(H5Gopen2(event,name.c_str(),H5P_DEFAULT),H5Gclose);H5 dataset(H5Dopen2(cluster,"cloud",H5P_DEFAULT),H5Dclose);size_t rows=0;const auto data=Matrix(dataset,5,rows);
            Cluster c;c.label=index;c.original_label=static_cast<int>(IntAttr(cluster,"label",index));
            const auto direction=IntAttr(cluster,"direction",-1);Require(direction>=-1&&direction<=1,"Invalid cluster direction");c.direction=static_cast<Direction>(direction);c.points.reserve(rows);
            for(size_t i=0;i<rows;++i){const auto* r=&data[i*5];CloudPoint p;p.x=r[0];p.y=r[1];p.z=r[2];p.integral=r[3];p.scale=r[4];c.points.push_back(p);}
            item.clusters.push_back(std::move(c));
        }
        output.events.push_back(std::move(item));
    }
    return output;
}
void WriteClusterH5(const std::string& path,const ClusterData& input){
    H5 file(H5Fcreate(path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);H5 group(H5Gcreate2(file,"cluster",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);PutInfo(group,input.info);
    for(const auto& event:input.events){
        const auto name="event_"+std::to_string(event.event);H5 item(H5Gcreate2(group,name.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);PutMeta(item,event.metadata,event.event);PutInt(item,"nclusters",event.clusters.size());
        for(size_t index=0;index<event.clusters.size();++index){
            const auto& c=event.clusters[index];const auto key="cluster_"+std::to_string(index);H5 cluster(H5Gcreate2(item,key.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
            PutInt(cluster,"label",c.original_label<0?c.label:c.original_label);PutInt(cluster,"direction",static_cast<int>(c.direction));
            std::vector<double> data;data.reserve(c.points.size()*5);for(const auto& p:c.points)data.insert(data.end(),{p.x,p.y,p.z,p.integral,p.scale});PutMatrix(cluster,"cloud",data,5);
        }
    }
}
}
#else
namespace rk4 {
PointcloudData ReadPointcloudH5(const std::string&){throw std::runtime_error("HDF5 support is disabled");}
void WritePointcloudH5(const std::string&,const PointcloudData&){throw std::runtime_error("HDF5 support is disabled");}
ClusterData ReadClusterH5(const std::string&){throw std::runtime_error("HDF5 support is disabled");}
void WriteClusterH5(const std::string&,const ClusterData&){throw std::runtime_error("HDF5 support is disabled");}
}
#endif
