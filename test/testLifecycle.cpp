#include "sim/Stopping.h"
#include "sim/KinematicsGenerator.h"
#include "sim/Simulator.h"
#include "reco/PointCloudBuilder.h"
#include <stdexcept>
static void Check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
struct InvalidDensity : Target {
    double density()const override{return 0;}
    double numberDensity()const override{return 0;}
    double arealDensity()const override{return 0;}
    double molarMass()const override{return 2;}
    std::vector<Element> compound()const override{return {{1,1,2}};}
    std::string str()const override{return "invalid";}
};
int main(){
    GasTarget gas({{1,1,2}},600);InvalidDensity invalid;
    StoppingPower stop;Check(stop.Init(12,4,gas),"Valid stopping init");
    Check(!stop.Init(12,4,invalid) && !stop.IsValid(),"Failed stopping init retained valid state");
    rk4::GeneratorConfig cfg;cfg.reaction={12,4,1,1,1,1};cfg.beam_ke_mev=150;
    rk4::KinematicsGenerator generator;rk4::Event event;
    Check(generator.Init(cfg,nullptr) && generator.Generate(event),"Valid generator init");
    cfg.reaction.beam_A=0;Check(!generator.Init(cfg,nullptr),"Invalid generator accepted");
    Check(!generator.Generate(event) && generator.Trials()==0,"Failed generator init still generates/trials");
    rk4::PadPlane pads;Check(pads.Load("data"),"Pads load");
    rk4::PointCloudBuilder builder;rk4::DetectorConfig detector;PeakFinderParameters peaks;
    Check(builder.Init(detector,&pads,peaks),"Valid builder init");
    Check(!builder.Init(detector,nullptr,peaks),"Invalid builder accepted");
    rk4::PointCloud cloud;rk4::TraceEvent raw;
    Check(!builder.Build(raw,cloud),"Failed builder init still builds");
    rk4::DetectorSimulator sim;Check(sim.Init(detector,&gas),"Valid detector init");
    detector.pad_dir="/nonexistent/spyral_test_data";
    Check(!sim.Init(detector,&gas),"Invalid detector accepted");
    Check(!sim.Simulate(event,cloud),"Failed detector init still simulates");
}
