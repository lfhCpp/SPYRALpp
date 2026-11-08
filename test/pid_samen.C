// 同 N 对比：只取两软件都接受的簇
//   左：Spyral 平面；中：SpyralPP 平面；右：叠在一张图上（黑=Spyral, 红=我们）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/ParticleId.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TText.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace rk4;

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

struct Row { double brho, dedx; };

static std::map<std::pair<int, int>, Row> Load(const char* path) {
    std::map<std::pair<int, int>, Row> m;
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() < 4) continue;
        m[{(int)v[0], (int)v[1]}] = Row{v[2], v[3]};
    }
    return m;
}

void pid_samen() {
    gStyle->SetOptStat(0);
    auto spy = Load("/tmp/spy_est5k.csv");
    auto ours = Load("/tmp/port_est5k.csv");
    std::vector<double> sb, sd, ob, od;
    for (const auto& kv : spy) {
        auto it = ours.find(kv.first);
        if (it == ours.end()) continue;
        sb.push_back(kv.second.brho); sd.push_back(std::sqrt(kv.second.dedx));
        ob.push_back(it->second.brho); od.push_back(std::sqrt(it->second.dedx));
    }
    std::printf("同 N 配对: %zu 簇\n", sb.size());
    if (sb.empty()) return;
    const double k = 1.0 / std::max(1e-9, Med(sd));

    TCanvas c("c", "same N", 1500, 500);
    c.Divide(3, 1);
    auto panel = [&](int idx, const char* title, const std::vector<double>& b,
                     const std::vector<double>& sq) {
        c.cd(idx + 1);
        gPad->SetLogz();
        TH2D* h = new TH2D(NULL, Form("%s;brho x100;#sqrt{dE/dx}", title), 120, 0, 120, 100, 0, 4);
        for (size_t i = 0; i < b.size(); ++i) h->Fill(b[i] * 100.0, k * sq[i]);
        h->SetMinimum(0.5);
        h->Draw("colz");
        TText t;
        t.SetNDC();
        t.SetTextSize(0.05);
        t.DrawText(0.28, 0.85, Form("N=%.0f", (double)b.size()));
    };
    panel(0, "Spyral (Python)", sb, sd);
    panel(1, "SpyralPP (C++)", ob, od);

    // 叠加
    c.cd(3);
    TH2D* frame = new TH2D("f", "overlay (black=Spyral, red=SpyralPP);brho x100;#sqrt{dE/dx}",
                           120, 0, 120, 100, 0, 4);
    TGraph gs, go;
    int n = 0;
    for (size_t i = 0; i < sb.size(); ++i) {
        gs.SetPoint(n, sb[i] * 100.0, k * sd[i]);
        go.SetPoint(n, ob[i] * 100.0, k * od[i]);
        ++n;
    }
    frame->Draw();
    gs.SetMarkerStyle(20);
    gs.SetMarkerSize(0.35);
    gs.SetMarkerColor(kBlack);
    gs.Draw("P same");
    go.SetMarkerStyle(20);
    go.SetMarkerSize(0.35);
    go.SetMarkerColor(kRed + 1);
    go.Draw("P same");
    TLegend leg(0.55, 0.75, 0.88, 0.88);
    leg.AddEntry(&gs, "Spyral", "p");
    leg.AddEntry(&go, "SpyralPP", "p");
    leg.Draw();
    c.SaveAs("/tmp/pid_samen.png");
    std::printf("图：/tmp/pid_samen.png\n");
}
