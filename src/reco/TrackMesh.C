#include "TrackMesh.h"
#include "../core/MassData.h"
#include "../core/Units.h"
#include "../sim/ParticleDynamics.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <regex>
#include <sstream>
#include <streambuf>
#include <cstdio>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace rk4 {
namespace {
// Increment when the propagator/material physics used to build a mesh changes.
constexpr int PhysicsRevision=2;
bool CountFor(const MeshConfig& c,size_t& count) {
    if(c.ke_bins<2||c.polar_bins<2||c.n_time<2||!(c.ke_min>0)||!(c.ke_max>c.ke_min)||
       !(c.polar_min>=0)||!(c.polar_max<=180)||!(c.polar_max>c.polar_min)||!(c.dt>0)||
       !std::isfinite(c.ke_min)||!std::isfinite(c.ke_max)||!std::isfinite(c.polar_min)||
       !std::isfinite(c.polar_max)||!std::isfinite(c.dt)||!std::isfinite(c.z_ref))return false;
    count=3;
    const size_t max=static_cast<size_t>(std::numeric_limits<std::streamsize>::max())/sizeof(double);
    for(int n:{c.n_time,c.ke_bins,c.polar_bins}){
        if(count>max/static_cast<size_t>(n))return false;
        count*=static_cast<size_t>(n);
    }
    return true;
}
std::string ByteOrder(){const unsigned int value=1;return *reinterpret_cast<const unsigned char*>(&value)?"little":"big";}
class MemoryBuffer : public std::streambuf {
    char* begin_;
public:
    MemoryBuffer(void* begin,size_t size):begin_(static_cast<char*>(begin)){setg(begin_,begin_,begin_+size);}
    size_t Position()const{return static_cast<size_t>(gptr()-begin_);}
};
}
size_t TrackMesh::Index(int t,int c,int k,int p)const {
    return ((static_cast<size_t>(t)*3+c)*cfg_.ke_bins+k)*cfg_.polar_bins+p;
}
void TrackMesh::ResetStorage(){
    if(map_base_)::munmap(map_base_,map_length_);
    map_base_=nullptr;map_length_=0;shared_view_=nullptr;shared_count_=0;
    data_.clear();has_context_=false;python_mesh_=false;composition_.clear();
}
TrackMesh::~TrackMesh(){ResetStorage();}
double TrackMesh::EnergyAt(int index)const {
    if(index==0)return cfg_.ke_min;
    if(index==cfg_.ke_bins-1)return cfg_.ke_max;
    const double fraction=static_cast<double>(index)/(cfg_.ke_bins-1);
    return cfg_.log_ke?std::exp(std::log(cfg_.ke_min)+fraction*std::log(cfg_.ke_max/cfg_.ke_min)):
                       cfg_.ke_min+fraction*(cfg_.ke_max-cfg_.ke_min);
}
bool TrackMesh::Build(const MeshConfig& cfg,int A,int Z,const Target& gas,double b,double e){
    ResetStorage();cfg_=cfg;particle_A_=A;particle_Z_=Z;
    size_t count=0;
    if(!CountFor(cfg_,count)||!std::isfinite(b)||!std::isfinite(e)||
       !(gas.density()>0)||!std::isfinite(gas.density())||A<1||Z<1||Z>A||!MassData::IsLoaded())return false;
    particle_mass_=MassData::GetNuclearMassNuclear(A,Z);
    if(!(particle_mass_>0))return false;
    composition_=gas.compound();
    if(!Utils::validCompound(composition_))return false;
    ParticleDynamics dynamics;
    dynamics.ConfigureTarget(composition_,gas.density());
    // Match Spyral 1.0's mesh convention: TrackGenerator stores positive field
    // magnitudes in metadata, then passes -B and -E to the equation of motion.
    // The detector config continues to use positive magnitudes at this API.
    dynamics.ConfigureFields(-b,-e);
    dynamics.SetStopZ(-.5,1.5);dynamics.SetStopRho(.5);dynamics.SetParticle(A,Z);
    try {data_.assign(count,0.0);}catch(const std::exception&){return false;}
    std::vector<TrackPoint> track;const double vertex[3]={0,0,cfg_.z_ref};
    for(int k=0;k<cfg_.ke_bins;++k)for(int p=0;p<cfg_.polar_bins;++p){
        const double angle=cfg_.polar_min+(cfg_.polar_max-cfg_.polar_min)*p/(cfg_.polar_bins-1);
        dynamics.Propagate(EnergyAt(k),vertex,angle*Units::DEG2RAD,0,track,cfg_.dt,cfg_.n_time);
        for(int t=0;t<cfg_.n_time;++t){
            const TrackPoint point=track.empty()?TrackPoint{}:track[std::min<size_t>(t,track.size()-1)];
            data_[Index(t,0,k,p)]=point.x*1000;data_[Index(t,1,k,p)]=point.y*1000;data_[Index(t,2,k,p)]=point.z*1000;
        }
    }
    bfield_=b;efield_=e;density_=gas.density();has_context_=true;return true;
}
bool TrackMesh::ReadHeader(std::istream& in){
    std::string magic;int version=0;
    if(!(in>>magic>>version)||magic!="RK4MESH"||version<1||version>3)return false;
    particle_A_=particle_Z_=0;
    if(version>=2&&!(in>>particle_A_>>particle_Z_))return false;
    if(!(in>>cfg_.n_time>>cfg_.ke_bins>>cfg_.polar_bins>>cfg_.ke_min>>cfg_.ke_max>>
         cfg_.polar_min>>cfg_.polar_max>>cfg_.dt>>cfg_.z_ref))return false;
    cfg_.log_ke=false;
    if(version==3){
        int revision=0,log=0,atoms=0;std::string order;
        if(!(in>>revision>>order>>log>>bfield_>>efield_>>density_>>particle_mass_>>atoms)||
           revision!=PhysicsRevision||order!=ByteOrder()||(log!=0&&log!=1)||atoms<1||atoms>10000||
           !std::isfinite(bfield_)||!std::isfinite(efield_)||!(density_>0)||!std::isfinite(density_)||
           !(particle_mass_>0)||!std::isfinite(particle_mass_))return false;
        cfg_.log_ke=log!=0;
        composition_.resize(static_cast<size_t>(atoms));
        for(auto& atom:composition_)if(!(in>>atom.a>>atom.z>>atom.s))return false;
        if(!Utils::validCompound(composition_))return false;
        has_context_=true;
    }
    size_t count;
    if(!CountFor(cfg_,count))return false;
    // The first binary byte can itself be whitespace. Consume only the header newline.
    char newline=0;if(!in.get(newline))return false;
    if(newline=='\r'&&!in.get(newline))return false;
    return newline=='\n';
}
bool TrackMesh::Load(const std::string& path){
    ResetStorage();std::ifstream in(path,std::ios::binary);
    if(!in||!ReadHeader(in))return false;
    size_t count=0;if(!CountFor(cfg_,count))return false;
    const auto offset=in.tellg();in.seekg(0,std::ios::end);const auto length=in.tellg();
    if(offset<0||length<offset||static_cast<uint64_t>(length-offset)<count*sizeof(double))return false;
    in.seekg(offset);
    std::vector<double> loaded;
    try {loaded.resize(count);}catch(const std::exception&){return false;}
    if(!in.read(reinterpret_cast<char*>(loaded.data()),static_cast<std::streamsize>(count*sizeof(double))))return false;
    if(!std::all_of(loaded.begin(),loaded.end(),[](double x){return std::isfinite(x);}))return false;
    data_.swap(loaded);return true;
}
bool TrackMesh::LoadPythonNpy(const std::string& path,const MeshConfig& cfg,int A,int Z,
                              double bfield,double efield){
    ResetStorage();cfg_=cfg;cfg_.log_ke=false;particle_A_=A;particle_Z_=Z;
    if(A<1||Z<1||Z>A||cfg.ke_bins<2||cfg.polar_bins<2||cfg.n_time<2||
       !std::isfinite(bfield)||!std::isfinite(efield))return false;
    // Python's cache filename is not a physics identity. Require its sidecar
    // and reject a same-shaped mesh generated at a different field/config.
    std::string meta_path=path;
    const auto extension=meta_path.rfind('.');
    if(extension!=std::string::npos)meta_path.resize(extension);
    meta_path+=".json";
    std::ifstream metadata(meta_path);
    if(!metadata)return false;
    const std::string json((std::istreambuf_iterator<char>(metadata)),{});
    const auto number=[&](const char* key,double& value){
        const std::regex pattern(std::string("\\\"")+key+"\\\"\\s*:\\s*(-?[0-9]+(?:\\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)");
        std::smatch match;if(!std::regex_search(json,match,pattern))return false;
        try{value=std::stod(match[1].str());return std::isfinite(value);}catch(...){return false;}
    };
    const auto matches=[&](const char* key,double expected){double value=0;
        return number(key,value)&&std::abs(value-expected)<=1e-12*std::max(1.0,std::abs(expected));};
    if(!matches("bfield",bfield)||!matches("efield",efield)||
       !matches("time_steps",cfg.n_time)||!matches("ke_min",cfg.ke_min)||
       !matches("ke_max",cfg.ke_max)||!matches("ke_bins",cfg.ke_bins)||
       !matches("polar_min",cfg.polar_min)||!matches("polar_max",cfg.polar_max)||
       !matches("polar_bins",cfg.polar_bins))return false;
    const int fd=::open(path.c_str(),O_RDONLY);if(fd<0)return false;struct stat st{};
    if(::fstat(fd,&st)!=0||st.st_size<16){::close(fd);return false;}
    const size_t length=static_cast<size_t>(st.st_size);void* base=::mmap(nullptr,length,PROT_READ,MAP_SHARED,fd,0);::close(fd);
    if(base==MAP_FAILED)return false;
    const auto* bytes=static_cast<const unsigned char*>(base);
    const unsigned char magic[6]={0x93,'N','U','M','P','Y'};
    if(std::memcmp(bytes,magic,6)!=0){::munmap(base,length);return false;}
    const unsigned major=bytes[6],minor=bytes[7];size_t offset=0,header=0;
    if(major==1){header=static_cast<size_t>(bytes[8])|(static_cast<size_t>(bytes[9])<<8);offset=10;}
    else if(major==2||major==3){if(length<12){::munmap(base,length);return false;}
        header=static_cast<size_t>(bytes[8])|(static_cast<size_t>(bytes[9])<<8)|(static_cast<size_t>(bytes[10])<<16)|(static_cast<size_t>(bytes[11])<<24);offset=12;}
    else{::munmap(base,length);return false;}
    if(header>length-offset){::munmap(base,length);return false;}
    const std::string descriptor(reinterpret_cast<const char*>(bytes+offset),header);
    if(descriptor.find("'descr': '<f8'")==std::string::npos||
       descriptor.find("'fortran_order': False")==std::string::npos){::munmap(base,length);return false;}
    const auto shape_key=descriptor.find("'shape': (");
    if(shape_key==std::string::npos){::munmap(base,length);return false;}
    const auto shape_end=descriptor.find(')',shape_key);
    if(shape_end==std::string::npos){::munmap(base,length);return false;}
    std::vector<size_t> dims;std::stringstream shape(descriptor.substr(shape_key+10,shape_end-shape_key-10));std::string token;
    while(std::getline(shape,token,',')){try{auto value=std::stoull(token);if(value)dims.push_back(value);}catch(...){}}
    if(dims.size()!=4||dims[0]!=static_cast<size_t>(cfg.n_time)||dims[1]!=3||
       dims[2]!=static_cast<size_t>(cfg.ke_bins)||dims[3]!=static_cast<size_t>(cfg.polar_bins)){
        ::munmap(base,length);return false;}
    size_t count=0;if(!CountFor(cfg_,count)||count>(length-offset-header)/sizeof(double)){
        ::munmap(base,length);return false;}
    map_base_=base;map_length_=length;shared_view_=reinterpret_cast<const char*>(bytes+offset+header);
    shared_count_=count;python_mesh_=true;return true;
}
bool TrackMesh::LoadShared(const std::string& path){
    ResetStorage();const int fd=::open(path.c_str(),O_RDONLY);if(fd<0)return false;
    struct stat st{};
    if(::fstat(fd,&st)!=0||st.st_size<=0){::close(fd);return false;}
    const size_t length=static_cast<size_t>(st.st_size);
    void* base=::mmap(nullptr,length,PROT_READ,MAP_SHARED,fd,0);::close(fd);
    if(base==MAP_FAILED)return false;
    MemoryBuffer buffer(base,length);std::istream in(&buffer);
    size_t count=0;
    if(!ReadHeader(in)||!CountFor(cfg_,count)||count>(length-buffer.Position())/sizeof(double)){
        ::munmap(base,length);return false;
    }
    // Keep a byte view: ASCII header length is not necessarily double aligned.
    const char* view=static_cast<const char*>(base)+buffer.Position();
    for(size_t i=0;i<count;++i){double value;std::memcpy(&value,view+i*sizeof(double),sizeof(double));
        if(!std::isfinite(value)){::munmap(base,length);return false;}}
    map_base_=base;map_length_=length;shared_count_=count;shared_view_=view;return true;
}
bool TrackMesh::Save(const std::string& path)const {
    if(Empty())return false;
    std::string temporary=path+".tmp.XXXXXX";
    std::vector<char> name(temporary.begin(),temporary.end());name.push_back('\0');
    const int fd=::mkstemp(name.data());if(fd<0)return false;::close(fd);
    temporary=name.data();
    std::ofstream out(temporary,std::ios::binary);
    if(!out){std::remove(temporary.c_str());return false;}
    out<<std::setprecision(17)<<"RK4MESH "<<(has_context_?3:2)<<'\n'<<particle_A_<<' '<<particle_Z_<<'\n'
       <<cfg_.n_time<<' '<<cfg_.ke_bins<<' '<<cfg_.polar_bins<<'\n'
       <<cfg_.ke_min<<' '<<cfg_.ke_max<<' '<<cfg_.polar_min<<' '<<cfg_.polar_max<<'\n'
       <<cfg_.dt<<' '<<cfg_.z_ref<<'\n';
    if(has_context_){
        out<<PhysicsRevision<<' '<<ByteOrder()<<' '<<cfg_.log_ke<<' '<<bfield_<<' '<<efield_<<' '
           <<density_<<' '<<particle_mass_<<' '<<composition_.size()<<'\n';
        for(const auto& atom:composition_)out<<atom.a<<' '<<atom.z<<' '<<atom.s<<'\n';
    }
    const char* bytes=shared_view_?shared_view_:reinterpret_cast<const char*>(data_.data());
    out.write(bytes,static_cast<std::streamsize>(Count()*sizeof(double)));out.close();
    if(!out||std::rename(temporary.c_str(),path.c_str())!=0){std::remove(temporary.c_str());return false;}
    return true;
}
bool TrackMesh::MatchesConfig(const MeshConfig& c,int A,int Z)const {
    return !Empty()&&particle_A_==A&&particle_Z_==Z&&cfg_.ke_bins==c.ke_bins&&
      cfg_.polar_bins==c.polar_bins&&cfg_.n_time==c.n_time&&cfg_.ke_min==c.ke_min&&
      cfg_.ke_max==c.ke_max&&cfg_.polar_min==c.polar_min&&cfg_.polar_max==c.polar_max&&
      cfg_.dt==c.dt&&cfg_.z_ref==c.z_ref&&cfg_.log_ke==c.log_ke;
}
bool TrackMesh::MatchesConfig(const MeshConfig& c,int A,int Z,const Target& gas,double b,double e)const {
    if(!has_context_||!MatchesConfig(c,A,Z)||bfield_!=b||efield_!=e||density_!=gas.density()||
       particle_mass_!=MassData::GetNuclearMassNuclear(A,Z))return false;
    const auto composition=gas.compound();if(composition.size()!=composition_.size())return false;
    for(size_t i=0;i<composition.size();++i)if(composition[i].a!=composition_[i].a||
        composition[i].z!=composition_[i].z||composition[i].s!=composition_[i].s)return false;
    return true;
}
bool TrackMesh::IsInRange(double ke,double angle)const {
    const double deg=angle/Units::DEG2RAD;
    return ke>=cfg_.ke_min&&ke<=cfg_.ke_max&&deg>=cfg_.polar_min&&deg<=cfg_.polar_max;
}
bool TrackMesh::Interpolate(double ke,double angle,std::vector<double>& xyz)const {
    if(Empty()||!IsInRange(ke,angle))return false;
    const double fk=(cfg_.log_ke?std::log(ke/cfg_.ke_min)/std::log(cfg_.ke_max/cfg_.ke_min):
                     (ke-cfg_.ke_min)/(cfg_.ke_max-cfg_.ke_min))*(cfg_.ke_bins-1);
    const double fp=(angle/Units::DEG2RAD-cfg_.polar_min)/(cfg_.polar_max-cfg_.polar_min)*(cfg_.polar_bins-1);
    const int k0=std::clamp(static_cast<int>(std::floor(fk)),0,cfg_.ke_bins-2);
    const int p0=std::clamp(static_cast<int>(std::floor(fp)),0,cfg_.polar_bins-2);
    const double wk=(ke-EnergyAt(k0))/(EnergyAt(k0+1)-EnergyAt(k0)),wp=fp-p0;
    xyz.resize(static_cast<size_t>(cfg_.n_time)*3);
    for(int t=0;t<cfg_.n_time;++t)for(int c=0;c<3;++c){
        const double lo=At(Index(t,c,k0,p0))*(1-wp)+At(Index(t,c,k0,p0+1))*wp;
        const double hi=At(Index(t,c,k0+1,p0))*(1-wp)+At(Index(t,c,k0+1,p0+1))*wp;
        xyz[static_cast<size_t>(t)*3+c]=lo*(1-wk)+hi*wk;
    }
    return true;
}
bool TrackMesh::InterpolatePython(double ke,double angle,std::vector<double>& xyz)const{
    if(!python_mesh_||Empty()||!IsInRange(ke,angle))return false;
    const double polar=angle;
    const double pwidth=(cfg_.polar_max*Units::DEG2RAD-cfg_.polar_min*Units::DEG2RAD)/cfg_.polar_bins;
    const double kwidth=(cfg_.ke_max-cfg_.ke_min)/cfg_.ke_bins;
    const double pf=(polar-cfg_.polar_min*Units::DEG2RAD)/pwidth;
    const double kf=(ke-cfg_.ke_min)/kwidth;
    int p0=static_cast<int>(std::floor(pf)),k0=static_cast<int>(std::floor(kf));
    // Python's bins-based cell width leaves an extrapolating fringe at the top edge.
    p0=std::clamp(p0,0,cfg_.polar_bins-2);k0=std::clamp(k0,0,cfg_.ke_bins-2);
    const double p1=std::min(p0+1,cfg_.polar_bins-1),k1=std::min(k0+1,cfg_.ke_bins-1);
    const double p_low=cfg_.polar_min*Units::DEG2RAD+p0*pwidth;
    const double k_low=cfg_.ke_min+k0*kwidth;
    const double wp=(polar-p_low)/pwidth,wk=(ke-k_low)/kwidth;
    xyz.resize(static_cast<size_t>(cfg_.n_time)*3);
    for(int t=0;t<cfg_.n_time;++t)for(int c=0;c<3;++c){
        const double q11=At(Index(t,c,k0,p0)),q12=At(Index(t,c,k0,static_cast<int>(p1)));
        const double q21=At(Index(t,c,static_cast<int>(k1),p0)),q22=At(Index(t,c,static_cast<int>(k1),static_cast<int>(p1)));
        xyz[static_cast<size_t>(t)*3+c]=q11*(1-wk)*(1-wp)+q21*wk*(1-wp)+q12*(1-wk)*wp+q22*wk*wp;
    }
    return true;
}
}
