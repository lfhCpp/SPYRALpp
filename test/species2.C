// 粒子假设判别：每簇在 5 个物种假设下取最小 |log(dE/dx / Bethe)|，按物种着色
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

static const int kA[5] = {1, 2, 3, 3, 4};
static const int kZ[5] = {1, 1, 1, 2, 2};
static const char* kName[5] = {"p", "d", "t", "3He", "alpha"};

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

static void Analyse(const char* tag, const std::vector<double>& b, const std::vector<double>& d,
                    double k, TGraph* g[5]) {
    // 先把 dEdx 换算到 Bethe 单位：k_dedx = median(dEdx / Bethe_p)
    std::vector<double> rr;
    for (size_t i = 0; i < b.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(b[i], 1, 1);
        if (e > 0) rr.push_back(d[i] / e);
    }
    std::sort(rr.begin(), rr.end());
    const double kd = rr.empty() ? 1.0 : rr[rr.size() / 2];
    int counts[5] = {0, 0, 0, 0, 0};
    int n[5] = {0, 0, 0, 0, 0};
    for (size_t i = 0; i < b.size(); ++i) {
        double best = 1e9;
        int bi = 0;
        for (int s = 0; s < 5; ++s) {
            const double e = ParticleId::ExpectedDedx(b[i], kA[s], kZ[s]);
            if (!(e > 0)) continue;
            const double chi = std::fabs(std::log((d[i] / kd) / e));
            if (chi < best) { best = chi; bi = s; }
        }
        counts[bi]++;
        g[bi]->SetPoint(n[bi]++, b[i] * 100.0, k * std::sqrt(d[i]));
    }
    std::printf("%s: ", tag);
    for (int s = 0; s < 5; ++s) std::printf("%s=%d (%.0f%%) ", kName[s], counts[s],
                                            100.0 * counts[s] / std::max<size_t>(1, b.size()));
    std::printf("\n");
}

void species2() {
    gStyle->SetOptStat(0);
    std::vector<double> sb, sd, ob, od;
    Load("/tmp/spy_est5k.csv", sb, sd);
    Load("/tmp/port_reinsch5k.csv", ob, od);
    std::vector<double> ssq;
    for (double d : sd) ssq.push_back(std::sqrt(d));
    const double k = 1.0 / std::max(1e-9, Med(ssq));

    TGraph* ga[5];
    TGraph* gb[5];
    const int cols[5] = {kRed + 1, kGreen + 2, kCyan + 1, kMagenta + 1, kOrange + 1};
    for (int s = 0; s < 5; ++s) {
        ga[s] = new TGraph();
        ga[s]->SetMarkerColor(cols[s]);
        ga[s]->SetMarkerStyle(20);
        ga[s]->SetMarkerSize(0.5);
        gb[s] = new TGraph();
        gb[s]->SetMarkerColor(cols[s]);
        gb[s]->SetMarkerStyle(20);
        gb[s]->SetMarkerSize(0.5);
    }
    Analyse("Spyral", sb, sd, k, ga);
    Analyse("SpyralPP", ob, od, k, gb);

    TCanvas c("c", "species", 1300, 560);
    c.Divide(2, 1);
    auto panel = [&](int idx, const char* title, TGraph* g[5]) {
        c.cd(idx + 1);
        TH2D* fr = new TH2D(NULL, Form("%s;brho x100;#sqrt{dE/dx}", title), 120, 0, 120, 100, 0, 4);
        fr->Draw();
        for (int s = 0; s < 5; ++s) g[s]->Draw("P same");
        TLegend leg(0.62, 0.62, 0.88, 0.88);
        for (int s = 0; s < 5; ++s) leg.AddEntry(g[s], kName[s], "p");
        leg.Draw();
    };
    panel(0, "Spyral (Python)", ga);
    panel(1, "SpyralPP (C++)", gb);
    c.SaveAs("/tmp/pid_species2.png");
    std::printf("图：/tmp/pid_species2.png\n");
}
