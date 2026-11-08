#include "reco/TrackMesh.h"
#include "core/Units.h"
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <cstring>
#include <cstdint>
#include <unistd.h>
static void Check(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
struct DeuteriumAtDensity:Target {
    double rho;explicit DeuteriumAtDensity(double value):rho(value){}
    double density()const override{return rho;}
    double numberDensity()const override{return 0;}
    double arealDensity()const override{return 0;}
    double molarMass()const override{return 4;}
    std::vector<Element> compound()const override{return {{2,1,2}};}
    std::string str()const override{return "deuterium";}
};
int main(){
    namespace fs=std::filesystem;
    const auto dir=fs::temp_directory_path()/("spyral_mesh_context_"+std::to_string(getpid()));fs::create_directories(dir);
    const std::string path=(dir/"mesh.bin").string();
    GasTarget gas({{1,1,2}},600),other_pressure({{1,1,2}},700);DeuteriumAtDensity other_species(gas.density());
    rk4::MeshConfig cfg;cfg.ke_bins=4;cfg.polar_bins=4;cfg.n_time=30;
    rk4::TrackMesh mesh;
    Check(mesh.Build(cfg,1,1,gas,2.85,60000),"Build physics mesh");
    Check(mesh.Save(path),"Save mesh");
    rk4::TrackMesh loaded,mapped;
    Check(loaded.Load(path)&&mapped.LoadShared(path),"Read owned/mapped mesh");
    for(const auto* m:{&loaded,&mapped}){
        Check(m->MatchesConfig(cfg,1,1,gas,2.85,60000),"Valid context rejected");
        Check(!m->MatchesConfig(cfg,1,1,gas,2.8,60000),"Changed magnetic field reused cache");
        Check(!m->MatchesConfig(cfg,1,1,gas,2.85,55000),"Changed electric field reused cache");
        Check(!m->MatchesConfig(cfg,1,1,other_pressure,2.85,60000),"Changed density reused cache");
        Check(!m->MatchesConfig(cfg,1,1,other_species,2.85,60000),"Changed gas species reused cache");
    }
    std::vector<double> a,b;
    Check(loaded.Interpolate(20,70*Units::DEG2RAD,a)&&mapped.Interpolate(20,70*Units::DEG2RAD,b)&&a==b,"Owned/mapped interpolation differs");
    Check(mapped.Save(path),"Atomic overwrite of mapped source failed");
    Check(mapped.Interpolate(20,70*Units::DEG2RAD,b)&&a==b,"Mapped source invalidated by overwrite");
    auto changed=cfg;changed.ke_bins=5;
    Check(mapped.Build(changed,1,1,gas,2.85,60000)&&!mapped.IsShared()&&mapped.MatchesConfig(changed,1,1,gas,2.85,60000),"Rebuild retained old mmap");
    changed.ke_max=changed.ke_min;
    Check(!mapped.Build(changed,1,1,gas,2.85,60000)&&mapped.Empty(),"Failed build retained stale mesh");
    const auto truncated=(dir/"truncated.bin").string();fs::copy_file(path,truncated);fs::resize_file(truncated,fs::file_size(truncated)-8);
    Check(!loaded.Load(truncated)&&loaded.Empty(),"Truncated owned cache accepted");
    Check(!mapped.LoadShared(truncated)&&mapped.Empty(),"Truncated mapped cache accepted");
    const auto legacy=(dir/"legacy.bin").string();
    {
        std::ofstream out(legacy,std::ios::binary);out<<"RK4MESH 2\n1 1\n2 2 2\n2 10 20 160\n1e-10 0.5\n";
        std::vector<double> payload(24,1.0);uint64_t bits=0x3ff000000000000aULL;std::memcpy(&payload[0],&bits,8);
        out.write(reinterpret_cast<const char*>(payload.data()),payload.size()*8);
    }
    Check(loaded.Load(legacy)&&mapped.LoadShared(legacy),"Legacy cache cannot be read");
    Check(loaded.Interpolate(2,20*Units::DEG2RAD,a)&&mapped.Interpolate(2,20*Units::DEG2RAD,b)&&a==b,"Binary whitespace byte skipped");
    rk4::MeshConfig old;old.ke_bins=old.polar_bins=old.n_time=2;old.ke_min=2;old.ke_max=10;old.polar_min=20;old.polar_max=160;old.log_ke=false;
    Check(!mapped.MatchesConfig(old,1,1,gas,2.85,60000),"Legacy cache claimed known physics context");
    {std::ofstream out(legacy,std::ios::binary);out<<"RK4MESH 2\n1 1\n-1 2 2\n2 10 20 160\n1e-10 0.5\n";}
    Check(!loaded.Load(legacy)&&loaded.Empty(),"Negative cache dimensions accepted");
    fs::remove_all(dir);
}
