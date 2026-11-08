// 最终对比：左 Spyral；右 = 我们（Bρ 来自我们的估计器，dE/dx 用与 Spyral 相同的 FITPACK 弧长）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
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
        if (v.size() < 2) continue;
        b.push_back(v[0]);
        d.push_back(v[1]);
    }
}

void pid_hybrid() {
    gStyle->SetOptStat(0);
    std::vector<double> sb, sd, ob, od;
    Load("/tmp/spy_est5k.csv", sb, sd);   // 注意：spy 文件前两列是 event,cluster → 需按 5 列读
    // 重新读取 spy（event,cluster,brho,dedx,...）
    sb.clear(); sd.clear();
    {
        std::ifstream in("/tmp/spy_est5k.csv");
        std::string line;
        std::getline(in, line);
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string f;
            std::vector<double> v;
            while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
            if (v.size() < 5) continue;
            sb.push_back(v[2]); sd.push_back(v[3]);
        }
    }
    Load("/tmp/ours_scipy_hybrid.csv", ob, od);
    std::vector<double> ssq, osq;
    for (double d : sd) ssq.push_back(std::sqrt(d));
    for (double d : od) osq.push_back(std::sqrt(d));
    const double k = 1.0 / std::max(1e-9, Med(ssq));
    std::printf("Spyral N=%zu | hybrid N=%zu\n", sb.size(), ob.size());

    TCanvas c("c", "hybrid", 1300, 560);
    c.Divide(2, 1);
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
        t.DrawText(0.30, 0.85, Form("N=%.0f", (double)b.size()));
    };
    panel(0, "Spyral (Python)", sb, ssq);
    panel(1, "SpyralPP (our B#rho + FITPACK arclength)", ob, osq);
    c.SaveAs("/tmp/pid_hybrid.png");
    std::printf("图：/tmp/pid_hybrid.png\n");
}
