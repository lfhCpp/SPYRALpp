#include "testPipeline.C"
void diagnose_pipeline(){
    const std::string dir="/tmp/spyral_fitter_audit";
    gSystem->mkdir(dir.c_str(),true);WriteConfig(dir+"/config.cfg");
    FILE* f=std::fopen((dir+"/config.cfg").c_str(),"a");std::fprintf(f,"\n[run]\nthreads=1\n[output]\nstages=1\n");std::fclose(f);
    const std::string command="build/rk4_simulate "+dir+"/config.cfg "+dir+"/clouds.root 5 > "+dir+"/simulate.log 2>&1";
    if(gSystem->Exec(command.c_str()))throw std::runtime_error("Simulation failed");
    if(gSystem->Exec(("build/rk4_analyze "+dir+"/clouds.root "+dir+"/fit.root "+dir+"/config.cfg > "+dir+"/analyze.log 2>&1").c_str()))throw std::runtime_error("Analysis failed");
    std::vector<Event> events;ReadKinematics(dir+"/clouds_kin.root",events);
    std::vector<FitRecord> records;ReadFitResults(dir+"/fit.root",records);
    for(const auto& e:events){
        std::printf("Event %lld: KE=%.8g, true brho=%.8g, polar=%.8g, vertex z=%.8g mm\n",e.id,e.particles[2].Kinetic(),e.particles[2].P()/299.8,std::acos(e.particles[2].pz/e.particles[2].P())*180/std::acos(-1.),e.vz*1000);
        for(const auto& r:records)if(r.event==e.id)std::printf("  fitted brho=%.8g, deviation=%.5g%%, cost=%.5g, z=%.8g, polar=%.8g, ok=%d\n",r.brho,100*(r.brho*299.8/e.particles[2].P()-1),r.cost,r.vz,r.polar*180/std::acos(-1.),r.ok);
    }
}
