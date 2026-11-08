// "一模一样"验证：逐簇对比 移植版 vs Spyral
//   左：dEdx_port vs dEdx_spyral（log-log，红线=1:1）
//   右：log10(比值) 分布
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TText.h"
#include "TStyle.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct R3 { double a, b, c; };

static void Load(const char* path, std::vector<R3>& v) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> r;
        while (std::getline(ss, f, ',')) r.push_back(std::atof(f.c_str()));
        if (r.size() < 5) continue;
        v.push_back({r[0] * 1e6 + r[1], r[2], r[3]}); // key, brho, dedx
    }
}

void same_check5k() {
    gStyle->SetOptStat(0);
    std::vector<R3> port, spy;
    Load("/tmp/port_est5k.csv", port);
    Load("/tmp/spy_est5k.csv", spy);
    std::sort(port.begin(), port.end(), [](const R3& x, const R3& y) { return x.a < y.a; });
    std::sort(spy.begin(), spy.end(), [](const R3& x, const R3& y) { return x.a < y.a; });

    std::vector<double> ratioL, bratio;
    TGraph g;
    int n = 0;
    for (const auto& p : port) {
        auto it = std::lower_bound(spy.begin(), spy.end(), p.a,
                                   [](const R3& x, double key) { return x.a < key; });
        if (it == spy.end() || std::fabs(it->a - p.a) > 0.5) continue;
        if (p.c <= 0 || it->c <= 0) continue;
        ratioL.push_back(std::log10(p.c / it->c));
        if (it->b > 0) bratio.push_back(p.b / it->b);
        g.SetPoint(n++, std::log10(it->c), std::log10(p.c));
    }
    std::sort(ratioL.begin(), ratioL.end());
    const double med = ratioL[ratioL.size() / 2];
    std::printf("配对 %d 个 | log10(dEdx比) 中位 %.4f  p16/84 %.3f/%.3f\n", n, med,
                ratioL[ratioL.size() * 16 / 100], ratioL[ratioL.size() * 84 / 100]);
    std::sort(bratio.begin(), bratio.end());
    std::printf("Bρ 比中位 %.3f\n", bratio[bratio.size() / 2]);

    TCanvas c("c", "same", 1300, 560);
    c.Divide(2, 1);
    c.cd(1);
    g.SetMarkerStyle(20);
    g.SetMarkerSize(0.5);
    g.SetMarkerColor(kBlue + 1);
    g.Draw("AP");
    g.GetXaxis()->SetTitle("log10(dE/dx) Spyral");
    g.GetYaxis()->SetTitle("log10(dE/dx) SpyralPP (port)");
    double lo = 200, hi = 1e4;
    TLine l1(std::log10(lo), std::log10(lo), std::log10(hi), std::log10(hi));
    l1.SetLineColor(kRed + 1);
    l1.SetLineWidth(2);
    l1.Draw("same");
    c.cd(2);
    TH1D* h = new TH1D("h", "log_{10}(dEdx_{port} / dEdx_{Spyral});log_{10} ratio;counts", 60, -0.3, 0.3);
    for (double v : ratioL) h->Fill(v);
    h->SetLineWidth(2);
    h->Draw("hist");
    TLine l2(0, 0, 0, h->GetMaximum() * 1.1);
    l2.SetLineColor(kRed + 1);
    l2.SetLineStyle(2);
    l2.Draw("same");
    c.SaveAs("/tmp/same_check5k.png");
    std::printf("图：/tmp/same_check.png\n");
}
