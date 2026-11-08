// 两种点云方法的 PID 平面并排：A = FFT+局部极大、B = FFT+电子学模板拟合
//   x = Bρ × 100，y = √(dE/dx)（用同一刻度归一化），红线 = 质子门 Bethe ±25%
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

struct Data {
    std::vector<double> brho, sq;
    size_t ncloud = 0;
};

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
            if (c.points.size() < 12) continue; // min_total_trajectory_points
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

static void Draw(TH2D* h, double k) {
    h->Draw("colz");
    for (int side = 0; side < 2; ++side) {
        TGraph* g = new TGraph();
        int n = 0;
        for (int i = 0; i < 240; ++i) {
            const double b = 0.05 + 0.01 * i;
            const double e = ParticleId::ExpectedDedx(b, 1, 1);
            if (!(e > 0)) continue;
            g->SetPoint(n++, b * 100.0, (side == 0 ? 0.75 : 1.25) * k * std::sqrt(e));
        }
        g->SetLineColor(kRed + 1);
        g->SetLineWidth(2);
        g->Draw("L same");
    }
}

void pid_ab() {
    gStyle->SetOptStat(0);
    const Data A = Collect("/tmp/abA_cloud.root");
    const Data B = Collect("/tmp/abB_cloud.root");
    if (A.brho.empty() || B.brho.empty()) { std::printf("数据为空\n"); return; }

    // 刻度：sqrt(dEdx) 中位归一（与 8993 那张图一致的取法）
    std::vector<double> sqb = B.sq;
    std::sort(sqb.begin(), sqb.end());
    const double k = 1.0 / std::max(1e-9, sqb[sqb.size() / 2]);

    std::printf("A: 点云 %zu → 估计 %zu  |  B: 点云 %zu → 估计 %zu | 刻度 k=%.4g\n", A.ncloud,
                A.brho.size(), B.ncloud, B.brho.size(), k);
    std::printf("√dEdx 中位（归一后）：A=%.2f  B=%.2f\n", k * Median(A.sq), k * Median(B.sq));

    TH2D* ha = new TH2D("ha", "A: FFT baseline + local maxima;brho x100;#sqrt{dE/dx} (norm)", 120,
                        0, 120, 100, 0, 3);
    TH2D* hb = new TH2D("hb", "B: FFT baseline + GET template fit;brho x100;#sqrt{dE/dx} (norm)",
                        120, 0, 120, 100, 0, 3);
    for (size_t i = 0; i < A.brho.size(); ++i) ha->Fill(A.brho[i] * 100.0, k * A.sq[i]);
    for (size_t i = 0; i < B.brho.size(); ++i) hb->Fill(B.brho[i] * 100.0, k * B.sq[i]);

    TCanvas c("c", "PID A/B", 1300, 560);
    c.Divide(2, 1);
    c.cd(1);
    gPad->SetLogz();
    ha->SetMinimum(0.5);
    Draw(ha, k);
    TText t1;
    t1.SetNDC();
    t1.SetTextSize(0.045);
    t1.DrawText(0.40, 0.82, Form("Total: %.0f", (double)A.brho.size()));
    c.cd(2);
    gPad->SetLogz();
    hb->SetMinimum(0.5);
    Draw(hb, k);
    TText t2;
    t2.SetNDC();
    t2.SetTextSize(0.045);
    t2.DrawText(0.40, 0.82, Form("Total: %.0f", (double)B.brho.size()));
    c.SaveAs("/tmp/pid_ab.png");
    std::printf("图：/tmp/pid_ab.png\n");
}
