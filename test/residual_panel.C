// 残差平面（沿参考线展开）：y = dEdx/(k·Bethe(Bρ))，band 应该是一条水平带
//   左 = Spyral (Python)，右 = SpyralPP (C++)；同一批簇、同一套极角切割
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/ParticleId.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLine.h"
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

static void Load(const char* path, bool ours, double th, std::vector<double>& b,
                 std::vector<double>& d) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        double bb, dd, pp;
        if (ours) {
            if (v.size() < 5) continue;
            bb = v[2]; dd = v[3]; pp = v[4];
        } else {
            if (v.size() < 3) continue;
            bb = v[0]; dd = v[1]; pp = v[2];
        }
        if (pp < th || !(bb > 0) || !(dd > 0)) continue;
        b.push_back(bb);
        d.push_back(dd);
    }
}

static void Panel(TCanvas& c, int idx, const char* title, const std::vector<double>& b,
                  const std::vector<double>& d) {
    std::vector<double> ratio;
    for (size_t i = 0; i < b.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(b[i], 1, 1);
        if (e > 0) ratio.push_back(d[i] / e);
    }
    const double k = 1.0 / std::max(1e-30, Med(ratio));
    c.cd(idx + 1);
    gPad->SetLogy();
    gPad->SetLogz();
    TH2D* h = new TH2D(NULL, Form("%s;brho x100;dEdx / (k#upointBethe)", title), 120, 0, 120, 80,
                       0.1, 10);
    std::vector<double> res;
    for (size_t i = 0; i < b.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(b[i], 1, 1);
        if (!(e > 0)) continue;
        const double rr = d[i] / (k * e);
        h->Fill(b[i] * 100.0, rr);
        res.push_back(std::log10(rr));
    }
    h->SetMinimum(0.5);
    h->Draw("colz");
    TLine* l = new TLine(0, 1, 120, 1);
    l->SetLineColor(kRed + 1);
    l->SetLineWidth(2);
    l->Draw("same");
    std::sort(res.begin(), res.end());
    const double q1 = res[res.size() / 4], q3 = res[3 * res.size() / 4];
    TText t;
    t.SetNDC();
    t.SetTextSize(0.045);
    t.DrawText(0.40, 0.84, Form("N=%.0f", (double)b.size()));
    t.DrawText(0.40, 0.76, Form("IQR = %.3f dex", q3 - q1));
}

void residual_panel() {
    gStyle->SetOptStat(0);
    std::vector<double> b1, d1, b2, d2;
    Load("/tmp/ours_est_sklearn.csv", true, 30.0, b1, d1);
    Load("/tmp/spy_est_ourclusters.csv", false, 30.0, b2, d2);
    if (b1.empty() || b2.empty()) { std::printf("空\n"); return; }
    TCanvas c("c", "residual", 1300, 560);
    c.Divide(2, 1);
    Panel(c, 1, "SpyralPP (C++)", b1, d1);
    Panel(c, 0, "Spyral (Python)", b2, d2);
    c.SaveAs("/tmp/residual_panel.png");
    std::printf("图：/tmp/residual_panel.png\n");
}
