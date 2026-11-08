//----------------------------------------------------------------------
// M1 I/O 测试：事件写入 ROOT TTree 再读回，逐字段比对
//   root -l -b -q test/testIO.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Event.h"
#include "io/TreeIO.h"
#include "sim/KinematicsGenerator.h"

#include "TSystem.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rk4;

static int g_fail = 0;

static void Check(bool ok, const char* msg) {
    if (!ok) {
        std::printf("  [FAIL] %s\n", msg);
        ++g_fail;
    }
}

static void CheckClose(double a, double b, double tol, const char* msg) {
    if (std::fabs(a - b) > tol) {
        std::printf("  [FAIL] %s: %.9g vs %.9g\n", msg, a, b);
        ++g_fail;
    }
}

void testIO() {
    std::printf("========== TTree I/O 测试 ==========\n");

    GeneratorConfig cfg;
    cfg.reaction.beam_A = 12;
    cfg.reaction.beam_Z = 4;
    cfg.reaction.target_A = 1;
    cfg.reaction.target_Z = 1;
    cfg.reaction.ejectile_A = 1;
    cfg.reaction.ejectile_Z = 1;
    cfg.beam_ke_mev = 150.0;
    cfg.ex_centroid = 1.0;
    cfg.ex_fwhm = 0.2;
    cfg.apply_target_loss = false;

    KinematicsGenerator gen;
    Check(gen.Init(cfg, nullptr, 2024), "Init");

    std::vector<Event> written;
    for (int i = 0; i < 100; ++i) {
        Event ev;
        if (!gen.Generate(ev)) {
            Check(false, "Generate 返回 false");
            gSystem->Exit(1);
        }
        written.push_back(ev);
    }

    const std::string path = "/tmp/rk4_kinematics_io_test.root";
    Check(WriteKinematics(path, written), "WriteKinematics");

    std::vector<Event> read_back;
    Check(ReadKinematics(path, read_back), "ReadKinematics");
    Check(read_back.size() == written.size(), "事件数一致");

    if (read_back.size() == written.size()) {
        int mismatched = 0;
        for (size_t i = 0; i < written.size(); ++i) {
            const Event& a = written[i];
            const Event& b = read_back[i];
            if (a.id != b.id || a.particles.size() != b.particles.size()) {
                ++mismatched;
                continue;
            }
            if (std::fabs(a.vx - b.vx) > 1e-12 || std::fabs(a.vy - b.vy) > 1e-12 ||
                std::fabs(a.vz - b.vz) > 1e-12 ||
                std::fabs(a.beam_energy - b.beam_energy) > 1e-12) {
                ++mismatched;
                continue;
            }
            for (size_t j = 0; j < a.particles.size(); ++j) {
                const Particle& pa = a.particles[j];
                const Particle& pb = b.particles[j];
                if (pa.A != pb.A || pa.Z != pb.Z || pa.label != pb.label ||
                    std::fabs(pa.px - pb.px) > 1e-12 ||
                    std::fabs(pa.py - pb.py) > 1e-12 ||
                    std::fabs(pa.pz - pb.pz) > 1e-12 ||
                    std::fabs(pa.E - pb.E) > 1e-12) {
                    ++mismatched;
                    break;
                }
            }
        }
        Check(mismatched == 0, "逐事件字段一致");
        if (mismatched != 0) {
            std::printf("    %d 条事件不一致\n", mismatched);
        }
    }

    // 抽查第一条
    if (!read_back.empty() && !read_back[0].particles.empty()) {
        const Particle& p = read_back[0].particles[2];
        CheckClose(p.Mass(), gen.Reaction().EjectileMass(), 1e-6, "读回后仍在壳");
    }

    gSystem->Unlink(path.c_str());

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
