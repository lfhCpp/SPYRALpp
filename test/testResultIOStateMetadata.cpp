#include "io/ResultIO.h"

#include "TFile.h"
#include "TTree.h"

#include <cstdio>
#include <stdexcept>
#include <unistd.h>

namespace {
void Check(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

std::string TempPath() {
    char name[] = "/tmp/spyralpp-result-state-XXXXXX";
    const int fd = mkstemp(name);
    if (fd < 0) throw std::runtime_error("mkstemp failed");
    close(fd);
    unlink(name);
    return name;
}

void WriteLegacyTree(const std::string& path) {
    TFile file(path.c_str(), "RECREATE");
    TTree tree("fit", "legacy fit results");
    long long event = 17;
    int cluster = 2, ok = 1, pid_index = -1;
    double brho = .2, sigma_brho = .01, polar = 1., azimuthal = 2.;
    double vx = 1., vy = 2., vz = 3., cost = 12., chi2ndf = 0., pid_sep = 0., dedx = 4.;
    tree.Branch("event", &event, "event/L");
    tree.Branch("cluster", &cluster, "cluster/I");
    tree.Branch("brho", &brho, "brho/D");
    tree.Branch("sigma_brho", &sigma_brho, "sigma_brho/D");
    tree.Branch("polar", &polar, "polar/D");
    tree.Branch("azimuthal", &azimuthal, "azimuthal/D");
    tree.Branch("vx", &vx, "vx/D");tree.Branch("vy", &vy, "vy/D");tree.Branch("vz", &vz, "vz/D");
    tree.Branch("cost", &cost, "cost/D");tree.Branch("chi2ndf", &chi2ndf, "chi2ndf/D");
    tree.Branch("ok", &ok, "ok/I");tree.Branch("pid_index", &pid_index, "pid_index/I");
    tree.Branch("pid_sep", &pid_sep, "pid_sep/D");tree.Branch("dedx", &dedx, "dedx/D");
    tree.Fill();tree.Write();file.Close();
}
}

int main() {
    const std::string legacy = TempPath(), current = TempPath();
    try {
        WriteLegacyTree(legacy);
        std::vector<rk4::FitRecord> records;
        Check(rk4::ReadFitResults(legacy, records), "legacy fit tree did not load");
        Check(records.size() == 1 && records[0].event == 17 && records[0].ok == 1,
              "legacy fit values changed while reading");
        Check(records[0].fit_state == static_cast<int>(rk4::FitState::NOT_RECORDED) &&
                  records[0].pid_state == static_cast<int>(rk4::PidState::NOT_RECORDED) &&
                  records[0].migrad_status == -1,
              "legacy tree must not invent PID or fit-state metadata");

        rk4::FitRecord row;
        row.event = 18;row.cluster = 1;row.ok = 1;row.brho = .3;row.cost = 19.;
        row.fit_state = static_cast<int>(rk4::FitState::FIT_ACCEPTED);
        row.pid_state = static_cast<int>(rk4::PidState::PID_UNKNOWN);
        row.pid_unknown_reason = static_cast<int>(rk4::PidUnknownReason::CALIBRATION_UNAVAILABLE);
        row.ic_observation_state = static_cast<int>(rk4::IcObservationState::MISSING_OBSERVATION);
        row.ic_amplitude = -1.;row.migrad_status = 1;row.migrad_edm = .01;
        Check(rk4::WriteFitResults(current, {row}), "new fit tree write failed");
        Check(rk4::ReadFitResults(current, records) && records.size() == 1,
              "new fit tree did not load");
        Check(records[0].fit_state == row.fit_state && records[0].pid_state == row.pid_state &&
                  records[0].pid_unknown_reason == row.pid_unknown_reason &&
                  records[0].ic_observation_state == row.ic_observation_state &&
                  records[0].migrad_status == 1,
              "new state metadata did not round-trip");
        unlink(legacy.c_str());unlink(current.c_str());
        return 0;
    } catch (...) {
        unlink(legacy.c_str());unlink(current.c_str());throw;
    }
}
