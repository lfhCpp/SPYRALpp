#include "io/CloudIO.h"
#include "TFile.h"
#include "TTree.h"
#include <filesystem>
#include <vector>
#include <limits>
#include <stdexcept>
#include <string>
#include <unistd.h>
static void Check(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
static void Malformed(const std::string& path,int kind){
    TFile file(path.c_str(),"RECREATE");TTree tree("cloud","cloud");
    Long64_t event=1;Int_t npoints=1;
    std::vector<double> x{1},y{2},z{3},amplitude{4},integral{5},tb{6},scale{1};
    std::vector<int> pad{7};
    if(kind==1)y.clear();
    if(kind==2)z[0]=std::numeric_limits<double>::quiet_NaN();
    tree.Branch("event",&event);tree.Branch("npoints",&npoints);
    tree.Branch("x",&x);tree.Branch("y",&y);tree.Branch("z",&z);
    tree.Branch("amplitude",&amplitude);tree.Branch("integral",&integral);
    tree.Branch("pad",&pad);tree.Branch("tb",&tb);
    if(kind!=0)tree.Branch("scale",&scale);
    tree.Fill();tree.Write();
}
int main(){
    const auto dir=std::filesystem::temp_directory_path()/("spyral_cloud_validation_"+std::to_string(getpid()));
    std::filesystem::create_directories(dir);
    std::vector<rk4::PointCloud> clouds(1);clouds[0].event=9;
    clouds[0].points.push_back(rk4::CloudPoint{});
    const auto good=(dir/"good.root").string();
    Check(rk4::WritePointClouds(good,clouds),"Cannot write valid cloud");
    Check(rk4::ReadPointClouds(good,clouds) && clouds.size()==1,"Cannot read valid cloud");
    for(int kind=0;kind<3;++kind){
        const auto bad=(dir/("bad"+std::to_string(kind)+".root")).string();Malformed(bad,kind);
        Check(!rk4::ReadPointClouds(bad,clouds),"Malformed cloud accepted");
        Check(clouds.empty(),"Failed read retained stale or partial clouds");
    }
    rk4::PointCloud pending;pending.pending_dither=true;
    const auto output=(dir/"pending.root").string();
    Check(!rk4::WritePointClouds(output,{pending}),"Pending dither serialized");
    Check(!std::filesystem::exists(output),"Invalid cloud created an output file");
    std::filesystem::remove_all(dir);
}
