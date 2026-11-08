R__ADD_INCLUDE_PATH(src)
R__ADD_INCLUDE_PATH(extern/fftw/api)
// test_waveprocessor.cpp
#include "src/reco/WaveProcessor.h"
#include <iostream>
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

// 从 event_waveform.txt 中读取指定事件的单个波形（默认 event 0, 1903_channel0_waveform）
static std::vector<double> readWaveform(const std::string& filename,
                                       int target_event = 0,
                                       const std::string& target_waveform = "1903_channel0_waveform",
                                       int expected_samples = 256) {
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
        
        // 检查事件分隔符
        if (line.find("==========") != std::string::npos && line.find("event") != std::string::npos) {
            // 解析事件号
            size_t start = line.find("event ") + 6;
            size_t end = line.find(" ", start);
            if (end == std::string::npos) end = line.length();
            
            std::string event_str = line.substr(start, end - start);
            current_event = std::stoi(event_str);
            
            in_target_event = (current_event == target_event);
            found_target_waveform = false;
            reading_waveform_data = false;
            continue;
        }
        
        if (!in_target_event) continue;
        
        // 检查是否是目标波形名称行
        if (!found_target_waveform) {
            if (line.find(target_waveform) != std::string::npos && line.find(":") != std::string::npos) {
                found_target_waveform = true;
                // 现在准备读取下一行的数据
                reading_waveform_data = true;
                continue;
            }
        }
        
        // 如果已找到目标波形名称，且这是数据行
        if (reading_waveform_data) {
            std::istringstream iss(line);
            double v;
            while (iss >> v && (int)waveform_data.size() < expected_samples) {
                waveform_data.push_back(v);
            }
            
            // 读取完一行后，重置状态
            reading_waveform_data = false;
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
    leg.Draw();


    c.SaveAs(out_file.c_str());
}

void workerThread(const WaveProcessor* processor, 
                  int thread_id, int num_traces) {
    std::cout << "Thread " << thread_id << " started, processing " 
              << num_traces << " traces" << std::endl;
    
    std::vector<double> input(processor->getTraceLength());
    std::vector<double> output(processor->getTraceLength());
    
    for (int i = 0; i < num_traces; i++) {
        // 生成测试数据
        for (size_t j = 0; j < input.size(); j++) {
            double t = (double)j / input.size();
            input[j] = sin(2 * M_PI * (5 + thread_id) * t) + 0.1 * t;
        }
        
        // 处理径迹
        if (!processor->process(input.data(), output)) {
            std::cerr << "Thread " << thread_id << ": Failed to process trace " << i << std::endl;
            continue;
        }
    }
    
    std::cout << "Thread " << thread_id << " completed" << std::endl;
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
    
    
    // 从 event_waveform.txt 读取波形并测试处理
    const std::string waveform_file = (argc > 1 ? argv[1] : "data/event_waveform.txt");
    int event_id = (argc > 2 ? std::stoi(argv[2]) : 0);
    event_id=4; // 强制测试 event 0，1903_channel0_waveform

    auto file_waveform = readWaveformPad(waveform_file, event_id, "1903_channel0_waveform", trace_length);
    std::vector<double> file_output(trace_length);
    std::vector<double> file_baseline(trace_length);
    if (file_waveform.size() == (size_t)trace_length) {
        std::cout << "First 5 input samples: ";
        for (int i = 0; i < 5 && i < trace_length; i++) {
            std::cout << file_waveform[i] << (i + 1 < 5 ? ", " : "\n");
        }

        if (processor->process(file_waveform.data(), file_output.data(), file_baseline.data())) {
            std::cout << "File waveform processing successful (" << waveform_file 
                      << ", event " << event_id << ")" << std::endl;
            std::cout << "First 5 output samples: ";
            for (int i = 0; i < 5 && i < trace_length; i++) {
                std::cout << file_output[i] << (i + 1 < 5 ? ", " : "\n");
            }

            // 检查输出统计
            int neg_count = 0;
            double min_val = *std::min_element(file_output.begin(), file_output.end());
            double max_val = *std::max_element(file_output.begin(), file_output.end());
            for (double v : file_output) {
                if (v < 0) neg_count++;
            }
            std::cout << "Output stats: min=" << min_val << ", max=" << max_val 
                      << ", negative_count=" << neg_count << "/" << trace_length << std::endl;

            // 使用 ROOT 绘图做调试（新增基线曲线）
            plotWaveforms(file_waveform, file_output, file_baseline, "WaveProcessor", "waveform_debug.png");
            std::cout << "Saved ROOT plot: waveform_debug.png" << std::endl;
        } else {
            std::cerr << "File waveform processing failed" << std::endl;
        }
    } else {
        std::cerr << "Skipping file waveform processing due to invalid sample length (" 
                  << file_waveform.size() << ")" << std::endl;
    }


    /*
    // 测试单线程处理
    // std::vector<double> test_input(trace_length);
    // std::vector<double> test_output;
    
    // for (int i = 0; i < trace_length; i++) {
    //     double t = (double)i / trace_length;
    //     test_input[i] =sin(2 * M_PI * 5 * t) + 0.2 * t;
    // }
    
    // if (processor->process(test_input.data(), test_output)) {
    //     std::cout << "Single trace processing successful" << std::endl;
    //     std::cout << "Output size: " << test_output.size() << std::endl;
    // } else {
    //     std::cerr << "Single trace processing failed" << std::endl;
    // }



    // 测试多线程处理
    const int num_threads = 4;
    const int traces_per_thread = 100;
    
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++) {
        threads.emplace_back(workerThread, processor.get(), i, traces_per_thread);
    }
    
    // 等待所有线程完成
    for (auto& t : threads) {
        t.join();
    }
    
    // 测试批量处理
    const int batch_size = 10;
    std::vector<const double*> inputs(batch_size);
    std::vector<double*> outputs(batch_size);
    std::vector<std::vector<double>> batch_data(batch_size, std::vector<double>(trace_length));
    std::vector<std::vector<double>> batch_results(batch_size, std::vector<double>(trace_length));
    
    for (int i = 0; i < batch_size; i++) {
        for (int j = 0; j < trace_length; j++) {
            double t = (double)j / trace_length;
            batch_data[i][j] =sin(2 * M_PI * (3 + i) * t) + 0.15 * t;
        }
        inputs[i] = batch_data[i].data();
        outputs[i] = batch_results[i].data();
    }
    
    int processed = processor->processBatch(inputs.data(), outputs.data(), batch_size);
    std::cout << "Batch processing: " << processed << "/" << batch_size 
              << " traces processed successfully" << std::endl;
    */
    return 0;
}