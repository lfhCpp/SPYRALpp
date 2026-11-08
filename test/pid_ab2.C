// 两种点云方法的 PID 平面并排（A=局部极大，B=电子学模板拟合）
// 关键：纵轴先按"沿路径 Bethe 平均"归一再回到 Bethe 曲线，数据才会落在红线上。
//   y = sqrt( dEdx / dedx_pred ) * sqrt(Bethe_p(Bρ)) / sqrt(median Bethe)
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
    std::vector<double> brho, ratio; // ratio = dEdx / dedx_pred
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
            if (c.points.size() < 30) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            est.SetTarget(1, 1, 2, 6.6168e-05); // 沿路径 Bethe 平均（H2 靶）
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0) || !(r.dedx_pred > 0)) continue;
            d.brho.push_back(r.brho);
            d.ratio.push_back(r.dEdx / r.dedx_pred);
        }
    }
    return d;
}

void pid_ab2() {
    gStyle->SetOptStat(0);
    const Data A = Collect("/tmp/abA_cloud.root");
    const Data B = Collect("/tmp/abB_cloud.root");
    if (A.brho.empty() || B.brho.empty()) { std::printf("数据为空\n"); return; }

    // 纵向刻度：以"质子峰"为 1（用 B 组的比值中位定）
    // y = sqrt(ratio/median(ratio)) * sqrt(Bethe)/median(sqrt Bethe)
    const double k = 1.0 / std::sqrt(std::max(1e-30, Median(B.ratio)));
    std::printf("A: 点云 %zu → 估计 %zu  |  B: 点云 %zu → 估计 %zu | 刻度 k=%.4g\n", A.ncloud,
                A.brho.size(), B.ncloud, B.brho.size(), k);

    // Bethe 归一化的公共刻度
    std::vector<double> sb;
    for (size_t i = 0; i < B.brho.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(B.brho[i], 1, 1);
        if (e > 0) sb.push_back(std::sqrt(e));
    }
    const double sbm = std::max(1e-9, Median(sb));

    auto make = [&](const char* title) {
        return new TH2D(Form("h%p", (void*)title), title, 120, 0, 120, 100, 0, 3);
    };
    TH2D* ha = make("A: FFT baseline + local maxima;brho x100;#sqrt{dE/dx} (norm)");
    TH2D* hb = make("B: FFT baseline + GET template fit;brho x100;#sqrt{dE/dx} (norm)");
    for (size_t i = 0; i < A.brho.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(A.brho[i], 1, 1);
        if (!(e > 0)) continue;
        ha->Fill(A.brho[i] * 100.0, k * std::sqrt(A.ratio[i]) * std::sqrt(e) / sbm);
    }
    for (size_t i = 0; i < B.brho.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(B.brho[i], 1, 1);
        if (!(e > 0)) continue;
        hb->Fill(B.brho[i] * 100.0, k * std::sqrt(B.ratio[i]) * std::sqrt(e) / sbm);
    }

    TCanvas c("c", "PID A/B", 1300, 560);
    c.Divide(2, 1);
    for (int idx = 0; idx < 2; ++idx) {
        c.cd(idx + 1);
        gPad->SetLogz();
        TH2D* h = (idx == 0) ? ha : hb;
        h->SetMinimum(0.5);
        h->Draw("colz");
        for (int side = 0; side < 3; ++side) {
            TGraph* g = new TGraph();
            int n = 0;
            const double fac = (side == 0 ? 0.75 : (side == 1 ? 1.0 : 1.25));
            for (int i = 0; i < 240; ++i) {
                const double b = 0.05 + 0.01 * i;
                const double e = ParticleId::ExpectedDedx(b, 1, 1);
                if (!(e > 0)) continue;
                g->SetPoint(n++, b * 100.0, fac * std::sqrt(e) / sbm);
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
