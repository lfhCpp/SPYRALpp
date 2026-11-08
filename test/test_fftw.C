R__ADD_INCLUDE_PATH(src)
R__ADD_INCLUDE_PATH(extern/fftw/api)
#include <fftw3.h>
#include <cmath>
#include <iostream>

// ROOT 头文件
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TPad.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"

int main(int argc, char **argv) {
    // 初始化 ROOT 应用程序
    TApplication app("FFT Visualization", &argc, argv);
    
    int n = 64;  // 信号长度
    double fs = 100.0;  // 采样频率 (Hz)
    double df = fs / n;  // 频率分辨率

    fftw_complex *in, *out;
    fftw_plan p;

    // 创建输入和输出数组
    in = fftw_alloc_complex(n);
    out = fftw_alloc_complex(n);

    // 初始化输入数组 - 包含两个频率成分的余弦信号
    double f1 = 10.0;  // 第一个频率 (Hz)
    double f2 = 30.0;  // 第二个频率 (Hz)
    double amplitude1 = 1.0;
    double amplitude2 = 0.5;
    
    for (int i = 0; i < n; i++) {
        double t = i / fs;  // 时间 (秒)
        // 创建包含两个频率成分的信号
        in[i][0] = amplitude1 * cos(2 * M_PI * f1 * t) 
                 + amplitude2 * cos(2 * M_PI * f2 * t);
        in[i][1] = 0;  // 虚部为0
    }

    // 创建FFT计划
    p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // 执行FFT变换
    fftw_execute(p);

    // 计算幅度谱
    double *magnitude = new double[n/2 + 1];
    double *frequency = new double[n/2 + 1];
    
    for (int i = 0; i <= n/2; i++) {
        frequency[i] = i * df;  // 频率轴
        magnitude[i] = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        
        // 对实信号，幅度需要乘以2（除了直流分量和Nyquist频率）
        if (i > 0 && i < n/2) {
            magnitude[i] *= 2.0;
        }
        
        // 归一化
        magnitude[i] /= n;
    }

    // 创建 ROOT 画布
    TCanvas *c1 = new TCanvas("c1", "FFT Analysis", 1200, 800);
    c1->Divide(2, 2);  // 分割为2x2四个区域

    // 1. 绘制时域信号
    c1->cd(1);
    TGraph *timeGraph = new TGraph(n);
    for (int i = 0; i < n; i++) {
        timeGraph->SetPoint(i, i/fs, in[i][0]);
    }
    timeGraph->SetTitle("Time Domain Signal");
    timeGraph->GetXaxis()->SetTitle("Time (s)");
    timeGraph->GetYaxis()->SetTitle("Amplitude");
    timeGraph->SetLineColor(kBlue);
    timeGraph->SetLineWidth(2);
    timeGraph->Draw("AL");

    // 2. 绘制幅度谱
    c1->cd(2);
    TGraph *freqGraph = new TGraph(n/2 + 1, frequency, magnitude);
    freqGraph->SetTitle("Frequency Domain (Magnitude Spectrum)");
    freqGraph->GetXaxis()->SetTitle("Frequency (Hz)");
    freqGraph->GetYaxis()->SetTitle("Magnitude");
    freqGraph->SetLineColor(kRed);
    freqGraph->SetLineWidth(2);
    freqGraph->SetMarkerStyle(20);
    freqGraph->SetMarkerSize(0.8);
    freqGraph->SetMarkerColor(kRed);
    freqGraph->Draw("APL");
    


    // 3. 绘制功率谱
    c1->cd(3);
    TGraph *powerGraph = new TGraph(n/2 + 1);
    for (int i = 0; i <= n/2; i++) {
        powerGraph->SetPoint(i, frequency[i], magnitude[i] * magnitude[i]);
    }
    powerGraph->SetTitle("Power Spectrum");
    powerGraph->GetXaxis()->SetTitle("Frequency (Hz)");
    powerGraph->GetYaxis()->SetTitle("Power");
    powerGraph->SetLineColor(kGreen+2);
    powerGraph->SetLineWidth(2);
    powerGraph->SetFillColor(kGreen-9);
    powerGraph->Draw("AB");

    // 4. 绘制相位谱
    c1->cd(4);
    TGraph *phaseGraph = new TGraph(n/2 + 1);
    for (int i = 0; i <= n/2; i++) {
        double phase = 0;
        if (magnitude[i] > 0) {
            phase = atan2(out[i][1], out[i][0]) * 180 / M_PI;
        }
        phaseGraph->SetPoint(i, frequency[i], phase);
    }
    phaseGraph->SetTitle("Phase Spectrum");
    phaseGraph->GetXaxis()->SetTitle("Frequency (Hz)");
    phaseGraph->GetYaxis()->SetTitle("Phase (degrees)");
    phaseGraph->SetLineColor(kMagenta);
    phaseGraph->SetLineWidth(2);
    phaseGraph->SetMarkerStyle(22);
    phaseGraph->SetMarkerColor(kMagenta);
    phaseGraph->Draw("APL");

    // 添加主标题
    TLatex *title = new TLatex();
    title->SetNDC();
    title->SetTextSize(0.03);
    title->DrawLatex(0.1, 0.97, "FFT Analysis of Multi-tone Signal (f1=10Hz, f2=30Hz)");

    // 在控制台输出结果摘要
    std::cout << "========================================" << std::endl;
    std::cout << "FFT Analysis Results:" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Signal length: " << n << " points" << std::endl;
    std::cout << "Sampling frequency: " << fs << " Hz" << std::endl;
    std::cout << "Frequency resolution: " << df << " Hz" << std::endl;
    std::cout << std::endl;
    std::cout << "Detected peaks:" << std::endl;
    
    // 查找峰值
    for (int i = 1; i < n/2; i++) {
        if (magnitude[i] > magnitude[i-1] && magnitude[i] > magnitude[i+1] && magnitude[i] > 0.1) {
            std::cout << "  Peak at " << frequency[i] << " Hz: Magnitude = " 
                      << magnitude[i] << std::endl;
        }
    }
    std::cout << "========================================" << std::endl;

    // 创建并拟合频域峰值
    c1->Update();
    
    // 创建一个新画布用于拟合
    TCanvas *c2 = new TCanvas("c2", "FFT Peak Fitting", 800, 600);
    
    // 使用TH1D来绘制频域以便于拟合
    TH1D *hist = new TH1D("hist", "Frequency Spectrum for Fitting", n/2, 0, fs/2);
    for (int i = 0; i <= n/2; i++) {
        if (i < n/2) {
            hist->SetBinContent(i+1, magnitude[i]);
        }
    }
    hist->SetTitle("FFT Magnitude Spectrum with Gaussian Fit");
    hist->GetXaxis()->SetTitle("Frequency (Hz)");
    hist->GetYaxis()->SetTitle("Magnitude");
    hist->SetLineColor(kBlue);
    hist->SetLineWidth(2);
    hist->SetFillColor(kBlue-9);
    hist->Draw("HIST");
    
    // 创建高斯函数来拟合峰值
    TF1 *fitFunc = new TF1("fitFunc", "gaus(0) + gaus(3)", 0, fs/2);
    fitFunc->SetParameters(0.5, 10, 2, 0.3, 30, 2);  // 初始参数估计
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
    
    // 执行拟合
    hist->Fit("fitFunc", "Q");
    
    // 添加图例
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(hist, "FFT Data", "f");
    leg->AddEntry(fitFunc, "Gaussian Fit", "l");
    leg->Draw();

    // 输出拟合结果
    std::cout << std::endl;
    std::cout << "Gaussian Fit Results:" << std::endl;
    std::cout << "Peak 1: Center = " << fitFunc->GetParameter(1) 
              << " Hz, Amplitude = " << fitFunc->GetParameter(0) << std::endl;
    std::cout << "Peak 2: Center = " << fitFunc->GetParameter(4) 
              << " Hz, Amplitude = " << fitFunc->GetParameter(3) << std::endl;
    std::cout << "========================================" << std::endl;

    // 运行ROOT应用程序
    app.Run();

    // 清理内存
    // delete[] magnitude;
    // delete[] frequency;
    // delete timeGraph;
    // delete freqGraph;
    // delete powerGraph;
    // delete phaseGraph;
    // delete title;
    // delete hist;
    // delete fitFunc;
    // delete leg;
    // delete c1;
    // delete c2;

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return 0;
}