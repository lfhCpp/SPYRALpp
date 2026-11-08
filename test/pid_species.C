// 多粒子种类的 PID 图：叠加 p/d/t/3He/α 的 Bethe 曲线（两软件并排）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/ParticleId.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLegend.h"
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

static void Load(const char* path, std::vector<double>& b, std::vector<double>& d) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() < 5) continue;
        b.push_back(v[2]);
        d.push_back(v[3]);
    }
}

static void DrawCurves(double k, int colorScheme) {
    const int A[5] = {1, 2, 3, 3, 4};
    const int Z[5] = {1, 1, 1, 2, 2};
    const int cols[5] = {kRed + 1, kGreen + 1, kCyan + 1, kMagenta + 1, kOrange + 1};
    for (int s = 0; s < 5; ++s) {
        TGraph* g = new TGraph();
        int n = 0;
        for (int i = 0; i < 240; ++i) {
            const double bb = 0.05 + 0.01 * i;
            const double e = ParticleId::ExpectedDedx(bb, A[s], Z[s]);
            if (e > 0) g->SetPoint(n++, bb * 100.0, k * std::sqrt(e));
        }
        g->SetLineColor(cols[s]);
        g->SetLineWidth(2);
        g->Draw("L same");
        (void)colorScheme;
    }
}

void pid_species() {
    gStyle->SetOptStat(0);
    std::vector<double> sb, sd, ob, od;
    Load("/tmp/spy_est5k.csv", sb, sd);
    Load("/tmp/port_reinsch5k.csv", ob, od);
    std::vector<double> ssq, osq;
    for (double d : sd) ssq.push_back(std::sqrt(d));
    for (double d : od) osq.push_back(std::sqrt(d));
    const double k = 1.0 / std::max(1e-9, Med(ssq));
    std::printf("Spyral N=%zu | SpyralPP N=%zu | k=%.4g\n", sb.size(), ob.size(), k);

    TCanvas c("c", "species", 1300, 560);
    c.Divide(2, 1);
    auto panel = [&](int idx, const char* title, const std::vector<double>& b,
                     const std::vector<double>& sq, double kk) {
        c.cd(idx + 1);
        gPad->SetLogz();
        TH2D* h = new TH2D(NULL, Form("%s;brho x100;#sqrt{dE/dx}", title), 120, 0, 120, 100, 0, 4);
        for (size_t i = 0; i < b.size(); ++i) h->Fill(b[i] * 100.0, kk * sq[i]);
        h->SetMinimum(0.5);
        h->Draw("colz");
        DrawCurves(kk, idx);
        TText t;
        t.SetNDC();
        t.SetTextSize(0.045);
        t.DrawText(0.30, 0.87, Form("N=%.0f", (double)b.size()));
        t.SetTextSize(0.035);
        t.DrawText(0.60, 0.55, "p,d,t,3He,#alpha");
    };
    panel(0, "Spyral (Python)", sb, ssq, k);
    panel(1, "SpyralPP (C++)", ob, osq, k);
    c.SaveAs("/tmp/pid_species.png");
    std::printf("图：/tmp/pid_species.png\n");

    // 统计有多少点靠近各条曲线（判别多带结构）
    const int A[5] = {1, 2, 3, 3, 4};
    const int Z[5] = {1, 1, 1, 2, 2};
    for (int s = 0; s < 5; ++s) {
        int cs = 0, co = 0;
        for (size_t i = 0; i < sb.size(); ++i) {
            const double e = ParticleId::ExpectedDedx(sb[i], A[s], Z[s]);
            if (e > 0 && std::fabs(std::log(sb.size() ? k * ssq[i] / (k * std::sqrt(e)) : 1.0)) < 0.15) ++cs;
        }
        for (size_t i = 0; i < ob.size(); ++i) {
            const double e = ParticleId::ExpectedDedx(ob[i], A[s], Z[s]);
            if (e > 0 && std::fabs(std::log(k * osq[i] / (std::sqrt(e) * 1.0 / 1.0))) < 0.15) ++co;
        }
        std::printf("  A=%d Z=%d: Spyral 靠近=%d, 我们靠近=%d\n", A[s], Z[s], cs, co);
    }
}
