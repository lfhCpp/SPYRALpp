#include "Simulator.h"

#include "../core/DataPath.h"
#include "../core/MassData.h"
#include "../core/Units.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace rk4 {

namespace {
constexpr long long kTbFactor = 4096; // pad * 4096 + tb
// A separate deterministic draw per (seed,event,pad,time bin) keeps detector
// charge statistics and parallel scheduling independent of time quantization.
double TimeDither(uint64_t seed, uint64_t event, uint64_t key) {
    uint64_t value = seed ^ (event * 0x9e3779b97f4a7c15ULL) ^ key;
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    value ^= value >> 31;
    return (value >> 11) * (1.0 / 9007199254740992.0);
}
}

bool DetectorSimulator::Init(const DetectorConfig& cfg, const Target* gas) {
    ready_ = false;
    cfg_ = cfg;

    if (!MassData::IsLoaded()) {
        std::fprintf(stderr,
                     "[DetectorSimulator] 核素质量表未加载"
                     "（把 data/ 放到工作目录，或设置 SPYRALPP_DATA_DIR）\n");
        return false;
    }

    if (cfg_.pad_dir.empty()) cfg_.pad_dir = DataDir();
    if (!pads_.Load(cfg_.pad_dir)) {
        std::fprintf(stderr, "[DetectorSimulator] pad 数据加载失败（目录 %s）\n",
                     cfg_.pad_dir.c_str());
        return false;
    }
    response_.Init(cfg_);

    if (gas != nullptr) {
        if (gas->compound().empty()) {
            std::fprintf(stderr, "[DetectorSimulator] 靶材组成为空\n");
            return false;
        }
        dynamics_.ConfigureTarget(*gas);
    }
    dynamics_.ConfigureFields(cfg_.bfield_t, cfg_.efield_vm);
    dynamics_.SetStopZ(0.0, cfg_.length_m);

    rng_.Seed(cfg_.seed);
    ready_ = true;
    return true;
}

void DetectorSimulator::Deposit(double x_mm, double y_mm, double tb, double electrons,
                                double sigma_mm) {
    const int tb_i = static_cast<int>(tb);
    if (tb_i < 0 || tb_i >= cfg_.num_tb) return;

    auto add = [&](int pad, double e) {
        if (pad < 0 || pads_.IsBeamPad(pad)) return;
        charge_[static_cast<long long>(pad) * kTbFactor + tb_i] += e;
    };

    if (sigma_mm < 0.05) {
        add(pads_.FindPad(x_mm, y_mm), electrons);
        return;
    }

    // 在 ±3 sigma 上铺 10x10 网格，按二维高斯分配电子
    const int steps = 10;
    const double half = 3.0 * sigma_mm;
    const double step = 2.0 * half / (steps - 1);
    const double norm = 1.0 / (2.0 * Units::PI * sigma_mm * sigma_mm);

    for (int i = 0; i < steps; ++i) {
        const double px = x_mm - half + i * step;
        for (int j = 0; j < steps; ++j) {
            const double py = y_mm - half + j * step;
            const double d2 = (px - x_mm) * (px - x_mm) + (py - y_mm) * (py - y_mm);
            const double w = norm * std::exp(-0.5 * d2 / (sigma_mm * sigma_mm)) *
                             step * step * electrons;
            if (w <= 0.0) continue;
            add(pads_.FindPad(px, py), w);
        }
    }
}

bool DetectorSimulator::Simulate(const Event& ev, PointCloud& pc,
                                 const std::vector<int>& indices) {
    pc.event = ev.id;
    pc.points.clear();
    charge_.clear();
    total_electrons_ = 0.0;
    if (!ready_) return false;

    std::vector<int> use = indices;
    if (use.empty()) {
        if (ev.particles.size() > 2) use.push_back(2);
        if (ev.particles.size() > 3) use.push_back(3);
    }

    const double dv = cfg_.length_m / (cfg_.windows_edge - cfg_.micromegas_edge);
    std::vector<TrackPoint> track;

    for (int idx : use) {
        if (idx < 0 || idx >= static_cast<int>(ev.particles.size())) continue;
        const Particle& part = ev.particles[idx];
        if (part.Z == 0) continue; // 中子不产生信号

        const double p = part.P();
        if (p <= 0.0) continue;
        double cos_theta = part.pz / p;
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
        const double theta = std::acos(cos_theta);
        // atan2 returns radians; ParticleDynamics::Propagate accepts radians.
        const double phi = std::atan2(part.py, part.px);

        double vtx[3] = {ev.vx, ev.vy, ev.vz};
        dynamics_.SetParticle(part.A, part.Z);
        dynamics_.Propagate(part.Kinetic(), vtx, theta, phi, track);
        if (track.size() < 2) continue;

        for (size_t i = 0; i + 1 < track.size(); ++i) {
            const double dE = track[i].ke - track[i + 1].ke;
            if (dE <= 0.0) continue;

            double n_electrons = dE * 1.0e6 / cfg_.w_value_ev;
            if (cfg_.fano > 0.0) {
                n_electrons = rng_.Gauss(n_electrons, std::sqrt(cfg_.fano * n_electrons));
            }
            if (n_electrons < 1.0) continue;
            n_electrons *= cfg_.mpgd_gain;

            const double x_mm = 0.5 * (track[i].x + track[i + 1].x) * 1000.0;
            const double y_mm = 0.5 * (track[i].y + track[i + 1].y) * 1000.0;
            const double z_m = 0.5 * (track[i].z + track[i + 1].z);
            if (z_m <= 0.0 || z_m >= cfg_.length_m) continue;

            const double tb = (cfg_.length_m - z_m) / dv + cfg_.micromegas_edge;
            if (tb < 0.0 || tb >= cfg_.num_tb) continue;

            const double sigma_m =
                std::sqrt(2.0 * cfg_.diffusion * dv * tb / cfg_.efield_vm);
            Deposit(x_mm, y_mm, tb, n_electrons, sigma_m * 1000.0);
        }
    }

    double total = 0.0;
    for (const auto& kv : charge_) {
        const int pad = static_cast<int>(kv.first / kTbFactor);
        const int tb = static_cast<int>(kv.first % kTbFactor);
        const double electrons = kv.second;
        total += electrons;

        double amplitude = 0.0, integral = 0.0;
        response_.Apply(electrons, amplitude, integral);
        if (amplitude <= cfg_.adc_threshold) continue;

        CloudPoint cp;
        cp.pad = pad;
        cp.tb = tb + (cfg_.dither_time_buckets ? TimeDither(cfg_.seed, ev.id, kv.first) : 0.0);
        cp.x = pads_.X(pad);
        cp.y = pads_.Y(pad);
        cp.z = (cfg_.windows_edge - cp.tb) /
               (cfg_.windows_edge - cfg_.micromegas_edge) * cfg_.length_m * 1000.0;
        cp.amplitude = amplitude;
        cp.integral = integral;
        cp.scale = pads_.Scale(pad);
        pc.points.push_back(cp);
    }
    total_electrons_ = total;

    std::sort(pc.points.begin(), pc.points.end(),
              [](const CloudPoint& a, const CloudPoint& b) { return a.z < b.z; });

    return true;
}

} // namespace rk4
