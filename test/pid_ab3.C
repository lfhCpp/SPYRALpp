// PID 平面对比（A=FFT+局部极大，B=FFT+电子学模板拟合）
// 与 8993 那张图同一套约定：
//   y = k·sqrt(dE/dx) / S，红线 = 0.75/1.00/1.25 · sqrt(Bethe_p(ρ)) / S
//   k = median( sqrt(Bethe) / sqrt(dEdx) )   （把带贴到 Bethe 曲线上）
//   S = median( sqrt(Bethe) )               （把纵轴落到 0..3）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "io/CloudIO.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TText.h"
#include "TStyle.h"
#include "TPad.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace rk4;

struct Data { std::vector<double> brho, sq; size_t ncloud = 0; };

static double Median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

static Data Collect(const char* path) {
    Data d;
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds(path, clouds)) { std::printf("读不到 %s\n", path); return d; }
    d.ncloud = clouds.size();
    Clusterer clusterer;
    ClusterParams cp;
    cp.eps_mm = 20.0;
    cp.min_points = 3;
    cp.min_cluster_size = 50;
    clusterer.SetParams(cp);
    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& c : clusters) {
            if (c.points.size() < 12) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            d.brho.push_back(r.brho);
            d.sq.push_back(std::sqrt(r.dEdx));
        }
    }
    return d;
}

void pid_ab3() {
    gStyle->SetOptStat(0);
    const Data A = Collect("/tmp/abA_cloud.root");
    const Data B = Collect("/tmp/abB_cloud.root");
    if (A.brho.empty() || B.brho.empty()) { std::printf("空\n"); return; }

    // k：把带贴到 Bethe 曲线；S：把纵轴落到 0..3
    std::vector<double> ratio, sb;
    for (size_t i = 0; i < B.brho.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(B.brho[i], 1, 1);
        if (!(e > 0)) continue;
        ratio.push_back(std::sqrt(e) / B.sq[i]);
        sb.push_back(std::sqrt(e));
    }
    const double k = Median(ratio);
    const double S = std::max(1e-9, Median(sb));
    std::printf("A: 点云 %zu → 估计 %zu | B: 点云 %zu → 估计 %zu | k=%.3f S=%.2f\n",
                A.ncloud, A.brho.size(), B.ncloud, B.brho.size(), k, S);
    std::printf("√dEdx 中位（归一后）：A=%.2f  B=%.2f\n", k * Median(A.sq) / S,
                k * Median(B.sq) / S);

    auto make = [&](const char* title) {
        auto* h = new TH2D(NULL, title, 120, 0, 120, 100, 0, 3);
        return h;
    };
    TH2D* ha = make("A: FFT baseline + local maxima;brho x100;#sqrt{dE/dx} (norm)");
    TH2D* hb = make("B: FFT baseline + GET template fit;brho x100;#sqrt{dE/dx} (norm)");
    for (size_t i = 0; i < A.brho.size(); ++i) ha->Fill(A.brho[i] * 100.0, k * A.sq[i] / S);
    for (size_t i = 0; i < B.brho.size(); ++i) hb->Fill(B.brho[i] * 100.0, k * B.sq[i] / S);

    TCanvas c("c", "PID A/B", 1300, 560);
    c.Divide(2, 1);
    for (int idx = 0; idx < 2; ++idx) {
        c.cd(idx + 1);
        gPad->SetLogz();
        TH2D* h = (idx == 0) ? ha : hb;
        h->SetMinimum(0.5);
        h->Draw("colz");
        for (int side = 0; side < 3; ++side) {
            const double fac = (side == 0 ? 0.75 : (side == 1 ? 1.0 : 1.25));
            TGraph* g = new TGraph();
            int n = 0;
            for (int i = 0; i < 240; ++i) {
                const double b = 0.05 + 0.01 * i;
                const double e = ParticleId::ExpectedDedx(b, 1, 1);
                if (!(e > 0)) continue;
                g->SetPoint(n++, b * 100.0, fac * std::sqrt(e) / S);
            }
            g->SetLineColor(kRed + 1);
            g->SetLineWidth(2);
            g->Draw("L same");
        }
        TText t;
        t.SetNDC();
        t.SetTextSize(0.045);
        t.DrawText(0.40, 0.82, Form("Total: %.0f", idx == 0 ? (double)A.brho.size()
                                                             : (double)B.brho.size()));
    }
    c.SaveAs("/tmp/pid_ab.png");
    std::printf("图：/tmp/pid_ab.png\n");
}
