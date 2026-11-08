// PID 平面：横轴 sqrt(dE/dx)，纵轴 Bρ —— 两种点云方法并排
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "reco/ParticleId.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace rk4;

struct Fits {
    std::vector<double> brho, dedx;
    std::vector<double> sq;
};

static Fits Read(const char* path) {
    Fits f;
    TFile file(path);
    if (file.IsZombie()) return f;
    TTree* t = (TTree*)file.Get("fit");
    if (!t) return f;
    double brho = 0, dedx = 0;
    t->SetBranchAddress("brho", &brho);
    t->SetBranchAddress("dedx", &dedx);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        if (!(brho > 0) || !(dedx > 0)) continue;
        f.brho.push_back(brho);
        f.dedx.push_back(dedx);
        f.sq.push_back(std::sqrt(dedx));
    }
    return f;
}

static double Median(std::vector<double> v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// 单常数刻度：让质子曲线穿过数据的 sqrt(dEdx) 中位
static double Calibrate(const Fits& f) {
    std::vector<double> ratio;
    for (size_t i = 0; i < f.brho.size(); ++i) {
        const double exp = ParticleId::ExpectedDedx(f.brho[i], 1, 1);
        if (exp > 0) ratio.push_back(std::sqrt(exp) / f.sq[i]);
    }
    return 1.0 / std::max(1e-9, Median(ratio)); // exp_sqrt = k * sqrt(Bethe)
}

static void DrawCurve(TH2D* frame, double k, int A, int Z, int color) {
    TGraph* g = new TGraph();
    int n = 0;
    for (int i = 0; i < 160; ++i) {
        const double brho = 0.05 + 0.02 * i; // T·m
        const double e = ParticleId::ExpectedDedx(brho, A, Z);
        if (!(e > 0)) continue;
        g->SetPoint(n++, brho, k * std::sqrt(e));
    }
    g->SetLineColor(color);
    g->SetLineWidth(2);
    g->Draw("L same");
}

static void Panel(const Fits& f, const char* tag) {
    if (f.brho.empty()) { std::printf("%s: 空\n", tag); return; }
    const double k = Calibrate(f);
    const double sqmax = 2.0 * Median(f.sq);
    std::printf("%s: 拟合=%zu  sqrt(dEdx) 中位=%.1f  刻度 k=%.3g\n", tag, f.brho.size(),
                Median(f.sq), k);

    TH2D* h = new TH2D("h", Form("%s;B#rho [T m];#sqrt{dE/dx} [a.u.]", tag), 160, 0.1, 2.4, 160,
                       0.0, sqmax);
    for (size_t i = 0; i < f.brho.size(); ++i) h->Fill(f.brho[i], f.sq[i]);
    h->Draw("colz");
    DrawCurve(h, k, 1, 1, kYellow);   // p
    DrawCurve(h, k, 2, 1, kGreen);    // d
    DrawCurve(h, k, 3, 1, kCyan);     // t
    DrawCurve(h, k, 3, 2, kMagenta);  // 3He
    DrawCurve(h, k, 4, 2, kRed);      // alpha
}

void pid_sqrt() {
    gStyle->SetOptStat(0);
    const Fits fa = Read("/tmp/pidA_fit.root");
    const Fits fb = Read("/tmp/pidB_fit.root");
    if (fa.brho.empty() || fb.brho.empty()) { std::printf("缺输入\n"); return; }

    TCanvas c("c", "PID sqrt(dEdx) vs Brho", 1300, 560);
    c.Divide(2, 1);
    c.cd(1);
    Panel(fa, "A: FFT baseline + local maxima");
    c.cd(2);
    Panel(fb, "B: FFT baseline + GET template fit");
    c.SaveAs("/tmp/pid_sqrt.png");
    std::printf("\n曲线：黄=p 绿=d 青=t 品红=3He 红=alpha（单常数刻度由质子中位定）\n");
}
