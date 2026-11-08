R__ADD_INCLUDE_PATH(src)
R__ADD_INCLUDE_PATH(extern/fftw/api)
#include "src/reco/WaveProcessor.h"
#include "src/reco/PeakFinder.h"
#include <iostream>
#include <random>
#include <vector>
#include <thread>
#include <chrono>
#include <fstream>
#include <sstream>
#include <math.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMarker.h>


// 从 event_waveform.txt 中读取指定事件的单个波形（默认 event 0, 1903_channel0_waveform）
static std::vector<double> readWaveformPad(const std::string& filename,
                                       int target_event = 0,
                                       const std::string& target_waveform = "1903_channel0_waveform",
                                       int expected_samples = 512) {
    std::vector<double> waveform_data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file '" << filename << "'" << std::endl;
        return waveform_data;
    }
    
    std::string line;
    int current_event = -1;
    bool in_target_event = false;
    bool found_target_waveform = false;  // 已找到目标波形名称
    bool reading_waveform_data = false;  // 正在读取波形数据
    
    while (std::getline(file, line) && (int)waveform_data.size() < expected_samples) {
        if (line.empty()) continue;
        
       if (line.find("==========") != std::string::npos && line.find("event") != std::string::npos) {
            // 解析事件号
            size_t start = line.find("event ") + 6;
            size_t end = line.find(" ", start);
            if (end == std::string::npos) end = line.length();
            
            std::string event_str = line.substr(start, end - start);
            current_event = std::stoi(event_str);
            std::getline(file, line); // 读取下一行，通常是 "waveform: 1903_channel0_waveform"
            std::getline(file, line);
            std::getline(file, line);
            std::getline(file, line);
            std::istringstream iss(line); 
            //std::cout << "Debug: Read line for event " << current_event << ": " << line << std::endl;
            std::string s;
            double v;
            iss >> s;iss >> s;
            iss >> v;iss >> v;
            iss >> v;iss >> v;
            iss >> v;
            while (iss >> v && (int)waveform_data.size() < expected_samples) {
                waveform_data.push_back(v);
            }
        break;  
        }
        
    
    }
    
    if ((int)waveform_data.size() != expected_samples) {
        std::cerr << "Warning: Expected " << expected_samples << " samples, got "
                  << waveform_data.size() << "\n";
    } else {
        std::cout << "Successfully read " << waveform_data.size() << " samples from event " 
                  << target_event << ", waveform: " << target_waveform << std::endl;
    }
    
    return waveform_data;
}


static void plotWaveforms(const std::vector<double>& input,
                          const std::vector<double>& output,
                          const std::vector<double>& baseline,
                          const std::vector<Peak>& peaks,
                          const std::string& title,
                          const std::string& out_file = "waveform_debug.png") {
    if (input.empty() || input.size() != output.size() || input.size() != baseline.size()) {
        std::cerr << "plotWaveforms: invalid input/output/baseline size" << std::endl;
        return;
    }

    // Ensure ROOT does not try to open a GUI window (useful for headless environments)
    gROOT->SetBatch(kTRUE);

    const int n = static_cast<int>(input.size());
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) x[i] = i;

    TGraph g_in(n, x.data(), input.data());
    TGraph g_out(n, x.data(), output.data());
    TGraph g_base(n, x.data(), baseline.data());

    g_in.SetLineColor(kRed);
    g_out.SetLineColor(kBlue);
    g_base.SetLineColor(kGreen+2);
    g_base.SetLineStyle(2);

    // 自动计算 y 轴范围（加一点边距）
    double ymin = input[0];
    double ymax = input[0];
    for (int i = 0; i < n; ++i) {
        ymin = std::min({ymin, input[i], output[i], baseline[i]});
        ymax = std::max({ymax, input[i], output[i], baseline[i]});
    }
    double yrange = ymax - ymin;
    if (yrange <= 0) {
        yrange = (std::abs(ymax) > 0 ? std::abs(ymax) : 1.0);
    }
    ymin -= 0.1 * yrange;
    ymax += 0.1 * yrange;

    g_in.SetTitle(title.c_str());
    g_in.GetXaxis()->SetTitle("Sample");
    g_in.GetYaxis()->SetTitle("Value");
    g_in.GetYaxis()->SetRangeUser(ymin, ymax);

    TCanvas c("c", "Waveform", 1200, 800);
    g_in.Draw("AL");
    g_out.Draw("L SAME");
    g_base.Draw("L SAME");

    TLegend leg(0.65, 0.70, 0.92, 0.88);
    leg.SetBorderSize(0);
    leg.AddEntry(&g_in, "Input", "l");
    leg.AddEntry(&g_out, "Output", "l");
    leg.AddEntry(&g_base, "Baseline", "l");
    TMarker dummy_peak(0, 0, 20);
    dummy_peak.SetMarkerColor(kBlack);
    leg.AddEntry(&dummy_peak, "Peaks", "p");
    TMarker dummy_inflection(0, 0, 21);
    dummy_inflection.SetMarkerColor(kRed);
    leg.AddEntry(&dummy_inflection, "Inflection Points", "p");
    leg.Draw();

    // Draw peaks
    for (const auto& p : peaks) {
        int idx = static_cast<int>(std::round(p.centroid));
        if (idx >= 0 && idx < n) {
            TMarker* marker = new TMarker(x[idx], output[idx], 20);
            marker->SetMarkerColor(kBlack);
            marker->SetMarkerSize(1.5);
            marker->Draw("SAME");
        }
    }

    // Draw inflection points
    for (const auto& p : peaks) {
        if (p.left_ip >= 0 && p.left_ip < n) {
            TMarker* left_marker = new TMarker(x[p.left_ip], output[p.left_ip], 21);
            left_marker->SetMarkerColor(kRed);
            left_marker->SetMarkerSize(1.0);
            left_marker->Draw("SAME");
        }
        if (p.right_ip >= 0 && p.right_ip < n) {
            TMarker* right_marker = new TMarker(x[p.right_ip], output[p.right_ip], 21);
            right_marker->SetMarkerColor(kRed);
            right_marker->SetMarkerSize(1.0);
            right_marker->Draw("SAME");
        }
    }

    c.SaveAs(out_file.c_str());
}

int main(int argc, char** argv) {

    const int trace_length = 512;
    const double window_scale =10;
    
    // 使用工厂函数创建预处理器
    auto processor = WaveProcessor::createUnique(trace_length, window_scale);
    if (!processor) {
        std::cerr << "Failed to create preprocessor" << std::endl;
        return 1;
    }
    
    std::cout << "Created preprocessor with trace_length=" << trace_length 
              << ", window_scale=" << window_scale << std::endl;


    // 1. 准备测试数据
    // 从 event_waveform.txt 读取波形并测试处理
    const std::string waveform_file = (argc > 1 ? argv[1] : "data/event_waveform.txt");
    int event_id = (argc > 2 ? std::stoi(argv[2]) : 0);
    event_id=4; // 强制测试 event 0，1903_channel0_waveform

    auto file_waveform = readWaveformPad(waveform_file, event_id, "1903_channel0_waveform", trace_length);
    std::vector<double> file_output(trace_length);
    std::vector<double> file_baseline(trace_length);

    processor->process(file_waveform.data(), file_output.data(), file_baseline.data());
    
    auto trace=file_output; // 直接使用处理后的输出作为测试数据

    // 2. 设置参数
    PeakFinderParameters params;
    params.peak_separation = 50.0;
    params.peak_prominence = 20.0;
    params.peak_max_width = 50.0;
    params.peak_threshold = 40.0;
    params.rel_height = 0.95;
    
    // 3. 创建峰值查找器
    PeakFinder finder(params);
    finder.setTrace(trace);
    
    // 4. 创建随机数生成器
    std::random_device rd;
    std::mt19937 rng(rd());
    
    // 5. 查找峰值
    std::vector<Peak> peaks = finder.findPeaks(rng);

    // 使用 ROOT 绘图做调试（新增基线曲线）
    plotWaveforms(file_waveform, file_output, file_baseline, peaks, "WaveProcessor", "waveform_debug.png");
    
    // 6. 输出结果
    std::cout << "找到 " << peaks.size() << " 个峰值:\n";
    for (size_t i = 0; i < peaks.size(); i++) {
        const Peak& p = peaks[i];
        std::cout << "  峰值 " << i << ":\n";
        std::cout << "  质心: " << p.centroid << "\n";
        std::cout << "  振幅: " << p.amplitude << "\n";
        std::cout << "  积分: " << p.integral << "\n";
        std::cout << "  左拐点: " << p.left_ip << "\n";
        std::cout << "  右拐点: " << p.right_ip << "\n";
        std::cout << "  宽度: " << (p.right_ip - p.left_ip) << "\n";
        std::cout << "  突出度: " << p.prominence << "\n";
    }
    
    return 0;
}