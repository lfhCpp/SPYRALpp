#include "ResultIO.h"

#include "TFile.h"
#include "TTree.h"

#include <cstdio>

namespace rk4 {
namespace {

void BindWrite(TTree& tree, FitRecord& r) {
    tree.Branch("event", &r.event, "event/L");
    tree.Branch("cluster", &r.cluster, "cluster/I");
    tree.Branch("brho", &r.brho, "brho/D");
    tree.Branch("sigma_brho", &r.sigma_brho, "sigma_brho/D");
    tree.Branch("polar", &r.polar, "polar/D");
    tree.Branch("azimuthal", &r.azimuthal, "azimuthal/D");
    tree.Branch("vx", &r.vx, "vx/D");
    tree.Branch("vy", &r.vy, "vy/D");
    tree.Branch("vz", &r.vz, "vz/D");
    tree.Branch("cost", &r.cost, "cost/D");
    tree.Branch("chi2ndf", &r.chi2ndf, "chi2ndf/D");
    tree.Branch("ok", &r.ok, "ok/I");
    tree.Branch("pid_index", &r.pid_index, "pid_index/I");
    tree.Branch("pid_sep", &r.pid_sep, "pid_sep/D");
    tree.Branch("dedx", &r.dedx, "dedx/D");
    tree.Branch("fit_state", &r.fit_state, "fit_state/I");
    tree.Branch("pid_state", &r.pid_state, "pid_state/I");
    tree.Branch("pid_unknown_reason", &r.pid_unknown_reason, "pid_unknown_reason/I");
    tree.Branch("ic_observation_state", &r.ic_observation_state, "ic_observation_state/I");
    tree.Branch("pid_state_basis", &r.pid_state_basis, "pid_state_basis/I");
    tree.Branch("ic_amplitude", &r.ic_amplitude, "ic_amplitude/D");
    tree.Branch("migrad_status", &r.migrad_status, "migrad_status/I");
    tree.Branch("migrad_edm", &r.migrad_edm, "migrad_edm/D");
}

void BindRead(TTree& tree, FitRecord& r) {
    tree.SetBranchAddress("event", &r.event);
    tree.SetBranchAddress("cluster", &r.cluster);
    tree.SetBranchAddress("brho", &r.brho);
    tree.SetBranchAddress("sigma_brho", &r.sigma_brho);
    tree.SetBranchAddress("polar", &r.polar);
    tree.SetBranchAddress("azimuthal", &r.azimuthal);
    tree.SetBranchAddress("vx", &r.vx);
    tree.SetBranchAddress("vy", &r.vy);
    tree.SetBranchAddress("vz", &r.vz);
    tree.SetBranchAddress("cost", &r.cost);
    tree.SetBranchAddress("chi2ndf", &r.chi2ndf);
    tree.SetBranchAddress("ok", &r.ok);
    tree.SetBranchAddress("pid_index", &r.pid_index);
    tree.SetBranchAddress("pid_sep", &r.pid_sep);
    tree.SetBranchAddress("dedx", &r.dedx);
    if (tree.GetBranch("fit_state")) tree.SetBranchAddress("fit_state", &r.fit_state);
    if (tree.GetBranch("pid_state")) tree.SetBranchAddress("pid_state", &r.pid_state);
    if (tree.GetBranch("pid_unknown_reason")) tree.SetBranchAddress("pid_unknown_reason", &r.pid_unknown_reason);
    if (tree.GetBranch("ic_observation_state")) tree.SetBranchAddress("ic_observation_state", &r.ic_observation_state);
    if (tree.GetBranch("pid_state_basis")) tree.SetBranchAddress("pid_state_basis", &r.pid_state_basis);
    if (tree.GetBranch("ic_amplitude")) tree.SetBranchAddress("ic_amplitude", &r.ic_amplitude);
    if (tree.GetBranch("migrad_status")) tree.SetBranchAddress("migrad_status", &r.migrad_status);
    if (tree.GetBranch("migrad_edm")) tree.SetBranchAddress("migrad_edm", &r.migrad_edm);
}

} // namespace

bool WriteFitResults(const std::string& path, const std::vector<FitRecord>& records,
                     const std::string& tree_name) {
    TFile file(path.c_str(), "RECREATE");
    if (file.IsZombie()) return false;

    TTree tree(tree_name.c_str(), "fit results");
    FitRecord r;
    BindWrite(tree, r);

    for (const auto& rec : records) {
        r = rec;
        tree.Fill();
    }
    tree.Write();
    file.Close();
    return true;
}

bool ReadFitResults(const std::string& path, std::vector<FitRecord>& records,
                    const std::string& tree_name) {
    TFile file(path.c_str(), "READ");
    if (file.IsZombie()) return false;

    TTree* tree = dynamic_cast<TTree*>(file.Get(tree_name.c_str()));
    if (tree == nullptr) return false;

    FitRecord r;
    BindRead(*tree, r);

    records.clear();
    const Long64_t n = tree->GetEntries();
    records.reserve(static_cast<size_t>(n));
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        records.push_back(r);
    }
    file.Close();
    return true;
}

} // namespace rk4
