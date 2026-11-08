// 估计阶段（已补齐 Spyral 的样条弧长 + 顶点/方向 cut）的 PID 质量 + 参考图风格
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

void pid_est() {
    gStyle->SetOptStat(0);
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/pidB_cloud.root", clouds)) return;

    Clusterer clusterer;
    std::vector<double> brho, dedx;
    int n_clusters = 0;

    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& cl : clusters) {
            if (cl.points.size() < 10) continue;
            ++n_clusters;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(10.0);
            EstimateResult r;
            if (!est.Estimate(cl, cloud, 2.85, r)) continue;   // 已含样条+顶点+方向 cut
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            brho.push_back(r.brho);
            dedx.push_back(r.dEdx);
        }
    }
    if (brho.empty()) { std::printf("没有通过估计的簇\n"); return; }

    std::vector<double> res;
    int in15 = 0, in25 = 0;
    double k = 0;
    {
        std::vector<double> ratio;
        for (size_t i = 0; i < brho.size(); ++i) {
            const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
            if (!(e > 0)) continue;
            ratio.push_back(std::sqrt(e) / std::sqrt(dedx[i]));
            const double rr = dedx[i] / e;
            res.push_back(std::log10(rr));
            if (std::fabs(rr - 1.0) < 0.15) ++in15;
            if (std::fabs(rr - 1.0) < 0.25) ++in25;
        }
        std::sort(ratio.begin(), ratio.end());
        k = ratio.empty() ? 1.0 : ratio[ratio.size() / 2];
    }
    std::sort(res.begin(), res.end());
    const double med = res[res.size() / 2];
    double s2 = 0;
    for (double v : res) s2 += (v - med) * (v - med);
    std::printf("簇 %d → 通过 %zu（%.0f%%）\n", n_clusters, brho.size(),
                100.0 * brho.size() / std::max(1, n_clusters));
    std::printf("Bethe 残差：中位 %+.2f dex，RMS %.3f dex，±15%% %.0f%%，±25%% %.0f%%\n",
                med, std::sqrt(s2 / res.size()), 100.0 * in15 / res.size(),
                100.0 * in25 / res.size());

    TH2D* h = new TH2D("h", "estimate stage (Spyral-style cuts);B#rho x100;#sqrt{dE/dx} norm",
                       120, 0, 120, 100, 0, 3);
    for (size_t i = 0; i < brho.size(); ++i) h->Fill(brho[i] * 100.0, k * std::sqrt(dedx[i]));
    TCanvas c("c", "pid", 700, 520);
    h->Draw("colz");
    c.SaveAs("/tmp/pid_est_style.png");
    std::printf("图：/tmp/pid_est_style.png\n");
}
