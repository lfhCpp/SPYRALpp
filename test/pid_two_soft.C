// 两软件 PID 对比：左 = Spyral(Python) estimate，右 = 我们(SpyralPP) estimate
// 数据：同一 300 事件子样本。纵轴 √(dE/dx) 各自归一化到自己的质子峰=1
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
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace rk4;

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

static void ReadCsv(const char* path, std::vector<double>& brho, std::vector<double>& dedx) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> vals;
        while (std::getline(ss, f, ',')) vals.push_back(std::atof(f.c_str()));
        if (vals.size() >= 3) { brho.push_back(vals[1]); dedx.push_back(vals[2]); }
    }
}

void pid_two_soft() {
    gStyle->SetOptStat(0);
    std::vector<double> sb, sd, ob, od;

    ReadCsv("/tmp/spyral_est.csv", sb, sd);

    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_2000.root", clouds)) { std::printf("our cloud fail\n"); return; }
    Clusterer cl;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> cs;
        cl.Run(cloud, cs);
        for (const auto& c : cs) {
            if (c.points.size() < 30) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            ob.push_back(r.brho);
            od.push_back(r.dEdx);
        }
    }
    std::printf("Spyral: %zu 条  |  SpyralPP: %zu 条\n", sb.size(), ob.size());
    if (sb.empty() || ob.empty()) return;

    std::vector<double> s_sq, o_sq;
    for (double d : sd) s_sq.push_back(std::sqrt(d));
    for (double d : od) o_sq.push_back(std::sqrt(d));
    const double ks = 1.0 / std::max(1e-9, Med(s_sq));
    const double ko = ks; // 两组共用同一刻度，直接可比
    std::printf("Spyral Bρ 中位 %.3f dEdx 中位 %.1f | 我们 Bρ 中位 %.3f dEdx 中位 %.1f\n",
                Med(sb), Med(sd), Med(ob), Med(od));

    TCanvas c("c", "PID comparison", 1300, 560);
    c.Divide(2, 1);

    auto panel = [&](int idx, const char* title, const std::vector<double>& b,
                     const std::vector<double>& sq, double k) {
        c.cd(idx + 1);
        gPad->SetLogz();
        TH2D* h = new TH2D(NULL, Form("%s;brho x100;#sqrt{dE/dx} (norm)", title), 120, 0, 120, 100,
                           0, 4);
        for (size_t i = 0; i < b.size(); ++i) h->Fill(b[i] * 100.0, k * sq[i]);
        h->SetMinimum(0.5);
        h->Draw("colz");
        TGraph g;
        int n = 0;
        for (int i = 0; i < 240; ++i) {
            const double bb = 0.05 + 0.01 * i;
            const double e = ParticleId::ExpectedDedx(bb, 1, 1);
            if (e > 0) g.SetPoint(n++, bb * 100.0, 0.25 * k * std::sqrt(e));
        }
        g.SetLineColor(kRed + 1);
        g.SetLineWidth(1);
        g.Draw("L same");
        TText t;
        t.SetNDC();
        t.SetTextSize(0.05);
        t.DrawText(0.35, 0.85, Form("N = %.0f", (double)b.size()));
    };

    panel(0, "Spyral (Python)", sb, s_sq, ks);
    panel(1, "SpyralPP (C++)", ob, o_sq, ko);
    c.SaveAs("/tmp/pid_two_soft.png");
    std::printf("图：/tmp/pid_two_soft.png\n");
}
