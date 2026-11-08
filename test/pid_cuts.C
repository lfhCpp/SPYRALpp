// 估计阶段的质量切割能不能救 PID？（切割前后 Bethe 残差宽度 + 门控后的图）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "io/CloudIO.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"

#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace rk4;

struct Row { double brho, dedx, len, npts, vrho; };

static void Report(const char* tag, const std::vector<Row>& rows) {
    std::vector<double> r;
    int in15 = 0;
    for (const auto& x : rows) {
        const double e = ParticleId::ExpectedDedx(x.brho, 1, 1);
        if (!(e > 0)) continue;
        r.push_back(std::log10(x.dedx / e));
        const double ratio = x.dedx / e;
        if (std::fabs(ratio - 1.0) < 0.15) ++in15;
    }
    if (r.empty()) { std::printf("%-28s 空\n", tag); return; }
    std::sort(r.begin(), r.end());
    const double med = r[r.size() / 2];
    double s2 = 0;
    for (double v : r) s2 += (v - med) * (v - med);
    std::printf("%-28s n=%zu 中位=%+.2f RMS=%.3f dex  ±15%%:%.0f%%\n", tag, r.size(), med,
                std::sqrt(s2 / r.size()), 100.0 * in15 / r.size());
}

void pid_cuts() {
    gStyle->SetOptStat(0);
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/pidB_cloud.root", clouds)) return;

    Clusterer clusterer;
    std::vector<Row> all, good;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& cl : clusters) {
            if (cl.points.size() < 12) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            EstimateResult r;
            if (!est.Estimate(cl, cloud, 2.85, r) || !(r.brho > 0) || !(r.dEdx > 0)) continue;
            Row row{r.brho, r.dEdx, r.arclength, (double)cl.points.size(),
                    std::hypot(r.vertex_x, r.vertex_y)};
            all.push_back(row);
            // Spyral 式质量切割
            if (row.npts < 20) continue;
            if (row.len < 40.0) continue;
            if (row.vrho > 10.0) continue;       // 顶点必须在束流轴附近
            if (r.radius > 2.0 * row.len) continue; // 短弧约束不住曲率 → Bρ 不可信
            if (std::fabs(std::log10(row.dedx / ParticleId::ExpectedDedx(row.brho, 1, 1))) > 0.5)
                continue;                         // 极端离群
            good.push_back(row);
        }
    }
    std::printf("切割前 %zu 簇，切割后 %zu 簇（保留 %.0f%%）\n\n", all.size(), good.size(),
                100.0 * good.size() / std::max<size_t>(1, all.size()));
    Report("切割前（全部）", all);
    Report("切割后（质量+门控）", good);

    // 参考图风格：x = Bρ×100，y = √(dEdx) 归一化到质子峰
    auto norm = [&](const std::vector<Row>& v) {
        std::vector<double> ratio;
        for (const auto& r : v) {
            const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
            if (e > 0) ratio.push_back(std::sqrt(e) / std::sqrt(r.dedx));
        }
        std::sort(ratio.begin(), ratio.end());
        return ratio.empty() ? 1.0 : ratio[ratio.size() / 2];
    };
    const double ka = norm(all), kg = norm(good);

    TCanvas c("c", "PID cuts", 1200, 500);
    c.Divide(2, 1);
    TH2D* h1 = new TH2D("h1", "before cuts;B#rho x100;sqrt(dE/dx) norm", 120, 0, 120, 100, 0, 3);
    TH2D* h2 = new TH2D("h2", "after quality cuts;B#rho x100;sqrt(dE/dx) norm", 120, 0, 120, 100, 0, 3);
    for (const auto& r : all) h1->Fill(r.brho * 100.0, ka * std::sqrt(r.dedx));
    for (const auto& r : good) h2->Fill(r.brho * 100.0, kg * std::sqrt(r.dedx));
    c.cd(1);
    h1->Draw("colz");
    c.cd(2);
    h2->Draw("colz");
    c.SaveAs("/tmp/pid_cuts.png");
    std::printf("\n图：/tmp/pid_cuts.png\n");
}
