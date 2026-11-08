//----------------------------------------------------------------------
// PID 对比：FFT 基线 + 局部极大（方法 A） vs FFT 基线 + 电子学响应模板拟合（方法 B）
//
//   root -l -b -q test/pidCompare.C
//
// 画两件事：
//   1) dE/dx–Bρ PID 平面（叠加 p/d/t/3He/alpha 的 Bethe 曲线）
//   2) 质子假设下的残差 log10(dE/dx ÷ Bethe_p)：越窄说明 PID 越干净
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "reco/ParticleId.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace rk4;

struct Fits {
    std::vector<double> brho, dedx, chi2ndf;
};

static Fits ReadFits(const char* path) {
    Fits f;
    TFile file(path);
    if (file.IsZombie()) return f;
    TTree* t = (TTree*)file.Get("fit");
    if (!t) return f;
    double brho = 0.0, dedx = 0.0, chi2 = 0.0;
    t->SetBranchAddress("brho", &brho);
    t->SetBranchAddress("dedx", &dedx);
    t->SetBranchAddress("chi2ndf", &chi2);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        if (!(brho > 0.0) || !(dedx > 0.0)) continue;
        f.brho.push_back(brho);
        f.dedx.push_back(dedx);
        f.chi2ndf.push_back(chi2);
    }
    return f;
}

static void Fill(TH2D* h, const Fits& f) {
    for (size_t i = 0; i < f.brho.size(); ++i) h->Fill(f.brho[i], f.dedx[i]);
}

// 质子假设的残差：log10(dedx / Bethe_p(Bρ))，只做形状比较（绝对刻度由数据定）
static void Residual(TH1D* h, const Fits& f) {
    for (size_t i = 0; i < f.brho.size(); ++i) {
        const double exp = ParticleId::ExpectedDedx(f.brho[i], 1, 1);
        if (!(exp > 0.0)) continue;
        h->Fill(std::log10(f.dedx[i] / exp));
    }
}

static double Median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

static void Summary(const char* tag, TH1D* h, const Fits& f) {
    // 用 ±2.5σ 的窗口估宽（先粗估中心，再算 RMS）
    const double center = h->GetBinCenter(h->GetMaximumBin());
    double sum = 0.0, sum2 = 0.0;
    int n = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double x = h->GetBinCenter(b);
        if (std::fabs(x - center) > 0.25) continue;
        const double c = h->GetBinContent(b);
        sum += c * x;
        sum2 += c * x * x;
        n += static_cast<int>(c);
    }
    const double mean = (n > 0) ? sum / n : 0.0;
    const double rms = (n > 0) ? std::sqrt(std::max(0.0, sum2 / n - mean * mean)) : 0.0;
    // 落在中心 ±0.1 dex（≈±26%）内的比例
    int inside = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        if (std::fabs(h->GetBinCenter(b) - center) <= 0.10) inside += (int)h->GetBinContent(b);
    }
    std::printf("%-22s 拟合=%zu  残差中心=%.2f  RMS=%.3f dex  ±0.1dex 内=%.0f%%  χ²/ndf 中位=%.2f\n",
                tag, f.brho.size(), mean, rms, 100.0 * inside / std::max(1, n),
                Median(f.chi2ndf));
}

void pidCompare() {
    gStyle->SetOptStat(0);

    const char* path_a = "/tmp/pid_a_fft.root";   // 方法 A：局部极大
    const char* path_b = "/tmp/pid_b_tf.root";    // 方法 B：电子学响应模板拟合
    const Fits fa = ReadFits(path_a);
    const Fits fb = ReadFits(path_b);
    if (fa.brho.empty() || fb.brho.empty()) {
        std::printf("读不到拟合结果：%s / %s\n", path_a, path_b);
        return;
    }

    // ---- 1) dE/dx - Bρ 平面（范围按数据 1%–99% 分位定）----
    std::vector<double> all;
    for (const auto& v : fa.dedx) all.push_back(v);
    for (const auto& v : fb.dedx) all.push_back(v);
    std::sort(all.begin(), all.end());
    const double lo = all[all.size() / 100];
    const double hi = all[all.size() * 99 / 100];
    std::printf("dE/dx 范围（1%%–99%%）：%.3g – %.3g charge/mm\n", lo, hi);

    TH2D* ha = new TH2D("ha", "A: FFT baseline + local maxima;B#rho [T m];dE/dx [charge/mm]",
                        120, 0.1, 2.2, 100, lo, hi);
    TH2D* hb = new TH2D("hb", "B: FFT baseline + GET template fit;B#rho [T m];dE/dx [charge/mm]",
                        120, 0.1, 2.2, 100, lo, hi);
    Fill(ha, fa);
    Fill(hb, fb);
    std::printf("PID 平面落入范围的点：A=%g  B=%g\n", ha->GetEntries(), hb->GetEntries());

    TCanvas c2("c2", "PID planes", 1200, 500);
    c2.Divide(2, 1);
    c2.cd(1);
    ha->Draw("colz");
    c2.cd(2);
    hb->Draw("colz");
    c2.SaveAs("/tmp/pid_planes.png");

    // ---- 2) 质子假设残差 ----
    TH1D* ra = new TH1D("ra", "A: local maxima;log_{10}(dE/dx / Bethe_{p});counts", 120, -0.8, 0.8);
    TH1D* rb = new TH1D("rb", "B: template fit;log_{10}(dE/dx / Bethe_{p});counts", 120, -0.8, 0.8);
    Residual(ra, fa);
    Residual(rb, fb);

    TCanvas c3("c3", "PID residual", 1000, 500);
    ra->SetLineColor(kAzure + 1);
    rb->SetLineColor(kRed + 1);
    ra->SetLineWidth(2);
    rb->SetLineWidth(2);
    ra->Draw("hist");
    rb->Draw("hist same");
    TLegend leg(0.62, 0.72, 0.88, 0.88);
    leg.AddEntry(ra, "A: 局部极大", "l");
    leg.AddEntry(rb, "B: 电子学模板拟合", "l");
    leg.Draw();
    c3.SaveAs("/tmp/pid_residual.png");

    std::printf("\n===== PID 质量对比（质子假设）=====\n");
    Summary("A: 局部极大", ra, fa);
    Summary("B: 模板拟合", rb, fb);
    std::printf("（残差 = log10(dE/dx ÷ Bethe_p(Bρ))；RMS 越小、±0.1dex 内比例越高，PID 越干净）\n");
}
