// PID 图（参考图风格）：x = Bρ×100，y = √(dE/dx) 归一化（质子峰=1），叠加质子门
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "io/CloudIO.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLine.h"
#include "TText.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace rk4;

void pid_final() {
    gStyle->SetOptStat(0);
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/spy_full_cloud.root", clouds)) return;

    Clusterer clusterer;
    ClusterParams cp;
    cp.eps_mm = 20;
    cp.min_points = 3;
    cp.min_cluster_size = 50;
    clusterer.SetParams(cp);

    std::vector<double> brho, sq, sqn;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> cs;
        clusterer.Run(cloud, cs);
        for (const auto& c : cs) {
            if (c.points.size() < 30) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            brho.push_back(r.brho);
            sq.push_back(std::sqrt(r.dEdx));
        }
    }
    if (brho.empty()) { std::printf("无数据\n"); return; }

    std::vector<double> ratio;
    for (size_t i = 0; i < brho.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
        if (e > 0) ratio.push_back(std::sqrt(e) / sq[i]);
    }
    std::sort(ratio.begin(), ratio.end());
    std::vector<double> sqc = sq;
    std::sort(sqc.begin(), sqc.end());
    const double k = 1.0 / std::max(1e-9, sqc[sqc.size() / 2]); // 质子峰归一化到 1
    for (size_t i = 0; i < brho.size(); ++i) sqn.push_back(k * sq[i]);
    std::printf("事件 %zu，刻度 k=%.3g\n", brho.size(), k);

    TH2D* h = new TH2D("h", "Particle ID (estimate stage, ADC threshold 40);brho x100;#sqrt{dE/dx} (norm)",
                       120, 0, 120, 100, 0, 3);
    for (size_t i = 0; i < brho.size(); ++i) h->Fill(brho[i] * 100.0, sqn[i]);

    TCanvas c("c", "pid", 760, 560);
    h->SetMinimum(0.5);
    gPad->SetLogz();
    h->Draw("colz");

    // 质子门：Bethe 曲线 ±25%
    for (int side = 0; side < 2; ++side) {
        TGraph* g = new TGraph();
        int n = 0;
        for (int i = 0; i < 240; ++i) {
            const double b = 0.05 + 0.01 * i;
            const double e = ParticleId::ExpectedDedx(b, 1, 1);
            if (!(e > 0)) continue;
            const double y = (side == 0 ? 0.75 : 1.25) * k * std::sqrt(e);
            g->SetPoint(n++, b * 100.0, y);
        }
        g->SetLineColor(kRed + 1);
        g->SetLineWidth(2);
        g->Draw("L same");
    }
    TText t;
    t.SetNDC();
    t.SetTextSize(0.045);
    t.DrawText(0.45, 0.80, Form("Total Events: %.0f", (double)brho.size()));
    t.SetTextColor(kRed + 1);
    t.DrawText(0.45, 0.72, "proton gate: Bethe #pm 25%");
    c.SaveAs("/tmp/pid_final.png");
    std::printf("图：/tmp/pid_final.png\n");
}
