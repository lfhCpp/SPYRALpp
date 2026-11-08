#include "TreeIO.h"

#include "TFile.h"
#include "TTree.h"

#include <cstdio>

namespace rk4 {
namespace {

struct BranchBuffer {
    Long64_t event = 0;
    Double_t vx = 0.0, vy = 0.0, vz = 0.0;
    Double_t beam_energy = 0.0;
    Int_t npart = 0;
    Int_t A[kMaxParticles] = {0};
    Int_t Z[kMaxParticles] = {0};
    Int_t label[kMaxParticles] = {0};
    Double_t px[kMaxParticles] = {0.0};
    Double_t py[kMaxParticles] = {0.0};
    Double_t pz[kMaxParticles] = {0.0};
    Double_t E[kMaxParticles] = {0.0};
};

void BindBranches(TTree& tree, BranchBuffer& b) {
    tree.Branch("event", &b.event, "event/L");
    tree.Branch("vx", &b.vx, "vx/D");
    tree.Branch("vy", &b.vy, "vy/D");
    tree.Branch("vz", &b.vz, "vz/D");
    tree.Branch("beam_energy", &b.beam_energy, "beam_energy/D");
    tree.Branch("npart", &b.npart, "npart/I");
    tree.Branch("A", b.A, "A[npart]/I");
    tree.Branch("Z", b.Z, "Z[npart]/I");
    tree.Branch("label", b.label, "label[npart]/I");
    tree.Branch("px", b.px, "px[npart]/D");
    tree.Branch("py", b.py, "py[npart]/D");
    tree.Branch("pz", b.pz, "pz[npart]/D");
    tree.Branch("E", b.E, "E[npart]/D");
}

// 读取时不能再用 Branch()（会在已有树上重复建分支），必须 SetBranchAddress
void BindReadBranches(TTree& tree, BranchBuffer& b) {
    tree.SetBranchAddress("event", &b.event);
    tree.SetBranchAddress("vx", &b.vx);
    tree.SetBranchAddress("vy", &b.vy);
    tree.SetBranchAddress("vz", &b.vz);
    tree.SetBranchAddress("beam_energy", &b.beam_energy);
    tree.SetBranchAddress("npart", &b.npart);
    tree.SetBranchAddress("A", b.A);
    tree.SetBranchAddress("Z", b.Z);
    tree.SetBranchAddress("label", b.label);
    tree.SetBranchAddress("px", b.px);
    tree.SetBranchAddress("py", b.py);
    tree.SetBranchAddress("pz", b.pz);
    tree.SetBranchAddress("E", b.E);
}

void FillBuffer(const Event& ev, BranchBuffer& b) {
    b.event = ev.id;
    b.vx = ev.vx;
    b.vy = ev.vy;
    b.vz = ev.vz;
    b.beam_energy = ev.beam_energy;

    int n = static_cast<int>(ev.particles.size());
    if (n > kMaxParticles) {
        std::fprintf(stderr, "[TreeIO] 事件 %lld 粒子数 %d 超过上限 %d，截断\n",
                     ev.id, n, kMaxParticles);
        n = kMaxParticles;
    }
    b.npart = n;
    for (int i = 0; i < n; ++i) {
        const Particle& p = ev.particles[i];
        b.A[i] = p.A;
        b.Z[i] = p.Z;
        b.label[i] = p.label;
        b.px[i] = p.px;
        b.py[i] = p.py;
        b.pz[i] = p.pz;
        b.E[i] = p.E;
    }
}

} // namespace

bool WriteKinematics(const std::string& path, const std::vector<Event>& events,
                     const std::string& tree_name) {
    TFile file(path.c_str(), "RECREATE");
    if (file.IsZombie()) {
        std::fprintf(stderr, "[TreeIO] 无法创建文件 %s\n", path.c_str());
        return false;
    }

    TTree tree(tree_name.c_str(), "AT-TPC kinematics events");
    BranchBuffer buf;
    BindBranches(tree, buf);

    for (const auto& ev : events) {
        FillBuffer(ev, buf);
        tree.Fill();
    }

    tree.Write();
    file.Close();
    return true;
}

bool ReadKinematics(const std::string& path, std::vector<Event>& events,
                    const std::string& tree_name) {
    TFile file(path.c_str(), "READ");
    if (file.IsZombie()) {
        std::fprintf(stderr, "[TreeIO] 无法打开文件 %s\n", path.c_str());
        return false;
    }

    TTree* tree = dynamic_cast<TTree*>(file.Get(tree_name.c_str()));
    if (tree == nullptr) {
        std::fprintf(stderr, "[TreeIO] 文件 %s 中没有树 %s\n", path.c_str(),
                     tree_name.c_str());
        return false;
    }

    BranchBuffer buf;
    BindReadBranches(*tree, buf);

    const Long64_t n_entries = tree->GetEntries();
    events.clear();
    events.reserve(static_cast<size_t>(n_entries));

    for (Long64_t i = 0; i < n_entries; ++i) {
        tree->GetEntry(i);
        Event ev;
        ev.id = buf.event;
        ev.vx = buf.vx;
        ev.vy = buf.vy;
        ev.vz = buf.vz;
        ev.beam_energy = buf.beam_energy;
        for (int j = 0; j < buf.npart; ++j) {
            Particle p;
            p.A = buf.A[j];
            p.Z = buf.Z[j];
            p.label = buf.label[j];
            p.px = buf.px[j];
            p.py = buf.py[j];
            p.pz = buf.pz[j];
            p.E = buf.E[j];
            ev.particles.push_back(p);
        }
        events.push_back(ev);
    }

    file.Close();
    return true;
}

} // namespace rk4
