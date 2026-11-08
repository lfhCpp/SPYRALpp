// 最终对比图：同一批簇（Spyral 聚出来的）分别用两个软件估计
//   左 = Spyral (Python) 的 estimate（/tmp/spyral_est.csv）
//   右 = SpyralPP (C++) 的 estimate（对 Spyral 的簇点集跑我们的 Estimator）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
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
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() >= 3) { brho.push_back(v[1]); dedx.push_back(v[2]); }
    }
}

void same_clusters() {
    gStyle->SetOptStat(0);
    std::vector<double> sb, sd, ob, od;
    ReadCsv("/tmp/spyral_est.csv", sb, sd);

    // 我们的估计：跑在 Spyral 的簇点集上（同一份 CSV）
    std::ifstream in("/tmp/spy_clusters.csv");
    std::string line;
    std::getline(in, line);
    std::vector<CloudPoint> cur;
    int cur_ev = -1, cur_cl = -1;
    auto flush = [&]() {
        if (cur.size() < 30) { cur.clear(); return; }
        Cluster c;
        c.direction = Direction::Forward;
        c.points = cur;
        Estimator est;
        est.SetSmallPadRadius(152.0);
        est.SetBeamRegionRadius(25.0);
        EstimateResult r;
        if (est.Estimate(c, PointCloud(), 2.85, r)) {
            if (r.brho > 0 && r.dEdx > 0) { ob.push_back(r.brho); od.push_back(r.dEdx); }
        }
        cur.clear();
    };
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() < 6) continue;
        const int ev = (int)v[0], cl = (int)v[1];
        if (ev != cur_ev || cl != cur_cl) { flush(); cur_ev = ev; cur_cl = cl; }
        CloudPoint p;
        p.x = v[2]; p.y = v[3]; p.z = v[4]; p.integral = v[5];
        cur.push_back(p);
    }
    flush();

    std::printf("同一批簇：Spyral N=%zu（Bρ中位 %.3f） | SpyralPP N=%zu（Bρ中位 %.3f）\n",
                sb.size(), Med(sb), ob.size(), Med(ob));
    if (sb.empty() || ob.empty()) return;

    std::vector<double> s_sq, o_sq;
    for (double d : sd) s_sq.push_back(std::sqrt(d));
    for (double d : od) o_sq.push_back(std::sqrt(d));
    const double k = 1.0 / std::max(1e-9, Med(s_sq));

    TCanvas c("c", "same clusters", 1300, 560);
    c.Divide(2, 1);
    auto panel = [&](int idx, const char* title, const std::vector<double>& b,
                     const std::vector<double>& sq) {
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
        g.Draw("L same");
        TText t;
        t.SetNDC();
        t.SetTextSize(0.05);
        t.DrawText(0.32, 0.85, Form("N = %.0f, B#rho_{med} = %.3f", (double)b.size(), Med(b)));
    };
    panel(0, "Spyral (Python) on its clusters", sb, s_sq);
    panel(1, "SpyralPP (C++) on the SAME clusters", ob, o_sq);
    c.SaveAs("/tmp/pid_same_clusters.png");
    std::printf("图：/tmp/pid_same_clusters.png\n");
}
