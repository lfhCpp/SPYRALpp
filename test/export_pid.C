#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
void export_pid(const char* input,const char* output,const char* tree_name="estimate") {
    TFile file(input,"READ");
    if (file.IsZombie()) throw std::runtime_error("Cannot open ROOT input");
    auto* tree=dynamic_cast<TTree*>(file.Get(tree_name));
    if (!tree) throw std::runtime_error("Missing PID tree");
    std::ofstream out(output);
    if (!out) throw std::runtime_error("Cannot open CSV output");
    out<<"event,cluster,ok,brho,dedx,polar_deg\n"<<std::setprecision(14);
    TTreeReader reader(tree);
    TTreeReaderValue<Long64_t> event(reader,"event");
    TTreeReaderValue<Int_t> cluster(reader,"cluster"),ok(reader,"ok");
    TTreeReaderValue<Double_t> brho(reader,"brho"),dedx(reader,"dedx"),polar(reader,"polar");
    Long64_t entries=0;
    while(reader.Next()) {
        out<<*event<<','<<*cluster<<','<<*ok<<','<<*brho<<','<<*dedx<<','<<*polar*180/std::acos(-1.0)<<'\n';
        ++entries;
    }
    out.close();
    if(!out || entries!=tree->GetEntries()) throw std::runtime_error("Incomplete PID export");
    std::printf("Exported %lld rows from %s to %s\n",entries,tree_name,output);
}
