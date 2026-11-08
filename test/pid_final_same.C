// 最终对比：同一批簇（我们的 sklearn HDBSCAN，与 Spyral 同参数）
//   左 = Spyral (Python) 估计器   右 = SpyralPP (C++) 估计器
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

static void Read(const char* path, std::vector<double>& b, std::vector<double>& d) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() >= 4) { b.push_back(v[2]); d.push_back(v[3]); }        // ours: ev,cl,brho,dedx
        else if (v.size() == 3) { b.push_back(v[0]); d.push_back(v[1]); }   // spyral: brho,dedx,polar
    }
}

void pid_final_same() {
    gStyle->SetOptStat(0);
    std::vector<double> sb, sd, ob, od;
    Read("/tmp/spy_est_ourclusters.csv", sb, sd);
    Read("/tmp/ours_est_sklearn.csv", ob, od);
    if (sb.empty() || ob.empty()) { std::printf("空数据\n"); return; }

    std::vector<double> s_sq, o_sq;
    for (double d : sd) s_sq.push_back(std::sqrt(d));
    for (double d : od) o_sq.push_back(std::sqrt(d));
    const double k = 1.0 / std::max(1e-9, Med(s_sq));          // 统一刻度
    std::printf("同一批簇：Spyral N=%zu Bρ中位=%.3f | SpyralPP N=%zu Bρ中位=%.3f\n", sb.size(),
                Med(sb), ob.size(), Med(ob));

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
        t.DrawText(0.30, 0.85, Form("N=%.0f  B#rho_{med}=%.3f", (double)b.size(), Med(b)));
    };
    panel(0, "Spyral (Python)", sb, s_sq);
    panel(1, "SpyralPP (C++)", ob, o_sq);
    c.SaveAs("/tmp/pid_final_same.png");
    std::printf("图：/tmp/pid_final_same.png\n");
}
