#include "StageIO.h"

#include "TFile.h"
#include "TTree.h"

#include <cstdio>

namespace rk4 {

bool WriteStages(const std::string& path, const std::vector<ClusterRecord>& clusters,
                 const std::vector<EstimateRecord>& estimates) {
    TFile file(path.c_str(), "RECREATE");
    if (file.IsZombie()) {
        std::fprintf(stderr, "[StageIO] 无法写 %s\n", path.c_str());
        return false;
    }

    // ---- 聚类阶段 ----
    {
        TTree tree("cluster", "clusters with their points");
        long long event = -1;
        int cluster = -1, npoints = 0, direction = -1;
        std::vector<double> *x = nullptr, *y = nullptr, *z = nullptr, *amp = nullptr,
                            *integ = nullptr, *tb = nullptr;

        tree.Branch("event", &event, "event/L");
        tree.Branch("cluster", &cluster, "cluster/I");
        tree.Branch("npoints", &npoints, "npoints/I");
        tree.Branch("direction", &direction, "direction/I");
        tree.Branch("x", &x);
        tree.Branch("y", &y);
        tree.Branch("z", &z);
        tree.Branch("amplitude", &amp);
        tree.Branch("integral", &integ);
        tree.Branch("tb", &tb);

        for (const auto& c : clusters) {
            event = c.event;
            cluster = c.cluster;
            npoints = c.npoints;
            direction = c.direction;
            x = const_cast<std::vector<double>*>(&c.x);
            y = const_cast<std::vector<double>*>(&c.y);
            z = const_cast<std::vector<double>*>(&c.z);
            amp = const_cast<std::vector<double>*>(&c.amplitude);
            integ = const_cast<std::vector<double>*>(&c.integral);
            tb = const_cast<std::vector<double>*>(&c.tb);
            tree.Fill();
        }
        tree.Write();
    }

    // ---- 估计阶段 ----
    {
        TTree tree("estimate", "estimator output");
        EstimateRecord r;
        tree.Branch("event", &r.event, "event/L");
        tree.Branch("cluster", &r.cluster, "cluster/I");
        tree.Branch("ok", &r.ok, "ok/I");
        tree.Branch("vx", &r.vx, "vx/D");
        tree.Branch("vy", &r.vy, "vy/D");
        tree.Branch("vz", &r.vz, "vz/D");
        tree.Branch("polar", &r.polar, "polar/D");
        tree.Branch("azimuthal", &r.azimuthal, "azimuthal/D");
        tree.Branch("brho", &r.brho, "brho/D");
        tree.Branch("radius", &r.radius, "radius/D");
        tree.Branch("dedx", &r.dEdx, "dedx/D");
        tree.Branch("arclength", &r.arclength, "arclength/D");
        tree.Branch("charge", &r.charge, "charge/D");
        tree.Branch("direction", &r.direction, "direction/I");

        for (const auto& e : estimates) {
            r = e;
            tree.Fill();
        }
        tree.Write();
    }

    file.Close();
    return true;
}

bool WriteTracePeaks(const std::string& path, const std::vector<TracePeakRecord>& peaks) {
    TFile file(path.c_str(), "RECREATE");
    if (file.IsZombie()) {
        std::fprintf(stderr, "[StageIO] 无法写 %s\n", path.c_str());
        return false;
    }

    TTree tree("peaks", "per-trace peaks (template fit or local maxima)");
    TracePeakRecord r;
    tree.Branch("event", &r.event, "event/L");
    tree.Branch("pad", &r.pad, "pad/I");
    tree.Branch("npeaks", &r.npeaks, "npeaks/I");
    tree.Branch("index", &r.index, "index/I");
    tree.Branch("time", &r.time, "time/D");
    tree.Branch("amplitude", &r.amplitude, "amplitude/D");
    tree.Branch("charge", &r.charge, "charge/D");
    tree.Branch("chi2", &r.chi2, "chi2/D");
    tree.Branch("npoints", &r.npoints, "npoints/I");
    tree.Branch("saturated", &r.saturated, "saturated/I");

    for (const auto& p : peaks) {
        r = p;
        tree.Fill();
    }
    tree.Write();
    file.Close();
    return true;
}

} // namespace rk4
