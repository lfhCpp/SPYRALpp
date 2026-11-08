#include "WaveProcessor.h"

#include <mutex>

namespace {
// FFTW 只有 fftw_execute 是线程安全的：建计划/销毁计划都会改全局 planner 状态。
// 这里把"规划"串行化，允许多个线程各自持有计划并并行 execute。
std::mutex g_fftw_plan_mutex;
} // namespace
#include <math.h>
#include <cstring>
#include <algorithm>
#include <iostream>

// 工厂函数实现
WaveProcessor* WaveProcessor::create(int trace_length, 
                                                                 double baseline_window_scale) {
    if (trace_length < 2 || baseline_window_scale <= 0) {
        return nullptr;
    }
    
    WaveProcessor* processor = new WaveProcessor(trace_length, 
                                                                            baseline_window_scale);
    if (!processor || !processor->initialize()) {
        delete processor;
        return nullptr;
    }
    
    return processor;
}

std::unique_ptr<WaveProcessor> 
WaveProcessor::createUnique(int trace_length, double baseline_window_scale) {
    return std::unique_ptr<WaveProcessor>(create(trace_length, baseline_window_scale));
}

// 私有构造函数
WaveProcessor::WaveProcessor(int trace_length, 
                                                        double baseline_window_scale)
    : trace_length_(trace_length)
    , baseline_window_scale_(baseline_window_scale)
    , half_length_(trace_length / 2 + 1)
    , filter_time_(nullptr)
    , filter_freq_(nullptr)
    , filter_plan_(nullptr)
    , forward_plan_(nullptr)
    , inverse_plan_(nullptr) {
}

// 初始化函数
bool WaveProcessor::initialize() {
    // 使用fftw_malloc分配对齐内存
    filter_time_ = (double*)fftw_malloc(trace_length_ * sizeof(double));
    filter_freq_ = (fftw_complex*)fftw_malloc(half_length_ * sizeof(fftw_complex));
    
    if (!filter_time_ || !filter_freq_) {
        return false;
    }
    
    // 创建sinc滤波器 sinc(t/T)<——> T·rect(T·f)
    createSincFilter();  
    
    // Python preprocess_traces uses ifftshift(sinc(window/scale)) directly
    // as the frequency multiplier. Do not FFT or sum-normalize this array.
    for (int i = 0; i < half_length_; ++i) {
        filter_freq_[i][0] = filter_time_[i];
        filter_freq_[i][1] = 0.0;
    }

    // 创建临时的缓冲区用于创建可重用的计划
    double* temp_time = (double*)fftw_malloc(trace_length_ * sizeof(double));
    fftw_complex* temp_freq = (fftw_complex*)fftw_malloc(half_length_ * sizeof(fftw_complex));
    
    if (!temp_time || !temp_freq) {
        if (temp_time) fftw_free(temp_time);
        if (temp_freq) fftw_free(temp_freq);
        return false;
    }
    
    // 创建可重用的FFTW计划
    {
        std::lock_guard<std::mutex> lk(g_fftw_plan_mutex);
        forward_plan_ =
            fftw_plan_dft_r2c_1d(trace_length_, temp_time, temp_freq, FFTW_ESTIMATE);
        inverse_plan_ =
            fftw_plan_dft_c2r_1d(trace_length_, temp_freq, temp_time, FFTW_ESTIMATE);
    }
    
    // 释放临时缓冲区
    fftw_free(temp_time);
    fftw_free(temp_freq);
    
    return forward_plan_ && inverse_plan_;
}

WaveProcessor::~WaveProcessor() {

    {
        std::lock_guard<std::mutex> lk(g_fftw_plan_mutex);
        if (forward_plan_) fftw_destroy_plan(forward_plan_);
        if (inverse_plan_) fftw_destroy_plan(inverse_plan_);
        if (filter_plan_) fftw_destroy_plan(filter_plan_);
    }
    
    if (filter_time_) fftw_free(filter_time_);
    if (filter_freq_) fftw_free(filter_freq_);
}

void WaveProcessor::createSincFilter() {
  
    int N = trace_length_;
    for (int i = 0; i < N; i++) {
        // 关键：直接计算ifftshift后的滤波器：FFT算法期望信号的零频在数组开头
        // 原始索引映射：i -> (i + N/2) % N - N/2
        double window_idx = static_cast<double>((i + N/2) % N) - N/2;
        double x = window_idx / baseline_window_scale_;
        
        if (std::fabs(x) < 1e-10) {
            filter_time_[i] = 1.0;
        } else {
            filter_time_[i] = sin(M_PI * x) / (M_PI * x);
        }
    }
     

}

bool WaveProcessor::computeStats(const double* data, int n, 
                                              double& mean, double& sigma) const {
    if (!data || n <= 0) {
        return false;
    }
    
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += data[i];
    }
    mean = sum / n;
    
    double sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = data[i] - mean;
        sum_sq += diff * diff;
    }
    
    sigma = std::sqrt(sum_sq / n);
    return true;
}

bool WaveProcessor::processInternal(const double* input, double* output,
                                                 double* time_buf, fftw_complex* freq_buf) const {
    if (!input || !output || !time_buf || !freq_buf) {
        return false;
    }
    
    // 1. 复制输入到输出（将同时作为工作数组）
    std::memcpy(output, input, trace_length_ * sizeof(double));
    
    // 2. 边缘平滑
    output[0] = output[1];
    output[trace_length_ - 1] = output[trace_length_ - 2];
    
    // 3. 创建基线数组并去除峰值
    std::memcpy(time_buf, output, trace_length_ * sizeof(double));
    
    double mean, sigma;
    if (!computeStats(time_buf, trace_length_, mean, sigma)) {
        return false;
    }
    
    // 找出非异常点并计算均值
    double valid_sum = 0.0;
    int valid_count = 0;
    for (int i = 0; i < trace_length_; i++) {
        if (time_buf[i] - mean <= 1.5 * sigma) {
            valid_sum += time_buf[i];
            valid_count++;
        }
    }
    
    if (valid_count == 0) {
        return false;  // 所有点都是异常点
    }
    
    double valid_mean = valid_sum / valid_count;
    
    // 替换异常点
    for (int i = 0; i < trace_length_; i++) {
        if (time_buf[i] - mean > 1.5 * sigma) {
            time_buf[i] = integer_input_ ? std::trunc(valid_mean) : valid_mean;
        }
    }
    
  
    // 4. 频域滤波
    // a) 执行正向FFT
    fftw_execute_dft_r2c(forward_plan_, time_buf, freq_buf);
    
    // b) 频域复数乘法
    for (int i = 0; i < half_length_; i++) {
        double a = freq_buf[i][0];
        double b = freq_buf[i][1];
        double c = filter_freq_[i][0];
        double d = filter_freq_[i][1];
        
        double real_part = a * c - b * d;
        double imag_part = a * d + b * c;
        
        freq_buf[i][0] = real_part;
        freq_buf[i][1] = imag_part;
    }

    // c) 执行逆向FFT
    fftw_execute_dft_c2r(inverse_plan_, freq_buf, time_buf);

    // d) 归一化（FFTW 逆变换不自动除以 N）
    for (int i = 0; i < trace_length_; i++) {
        time_buf[i] /= trace_length_;  
        //std::cout << "time_buf[" << i << "] = " << time_buf[i] << std::endl;  // 调试输出
    }
    
    // e) 从原始信号中减去滤波后的基线
    for (int i = 0; i < trace_length_; i++) {
        output[i] -= time_buf[i];
    }
    
    return true;
}

bool WaveProcessor::process(const double* input, double* output)  const {
    return process(input, output, nullptr);
}

bool WaveProcessor::process(const double* input, double* output, double* baseline) const {
    if (!input || !output) {
        return false;
    }
    
    // 为每个线程调用分配独立的缓冲区
    double* time_buf = (double*)fftw_malloc(trace_length_ * sizeof(double));
    fftw_complex* freq_buf = (fftw_complex*)fftw_malloc(half_length_ * sizeof(fftw_complex));
    
    if (!time_buf || !freq_buf) {
        if (time_buf) fftw_free(time_buf);
        if (freq_buf) fftw_free(freq_buf);
        return false;
    }
    
    bool result = processInternal(input, output, time_buf, freq_buf);

    if (result && baseline) {
        std::memcpy(baseline, time_buf, trace_length_ * sizeof(double));
    }
    
    fftw_free(time_buf);
    fftw_free(freq_buf);
    
    return result;
}

bool WaveProcessor::process(const double* input, std::vector<double>& output) const {
    if (!input) {
        return false;
    }
    
    output.resize(trace_length_);
    return process(input, output.data());
}

int WaveProcessor::processBatch(const double** inputs, double** outputs, 
                                             int num_traces) const {
    if (!inputs || !outputs || num_traces <= 0) {
        return 0;
    }
    
    int success_count = 0;
    
    // 为每个波形独立处理
    for (int i = 0; i < num_traces; i++) {
        if (process(inputs[i], outputs[i])) {
            success_count++;
        }
    }
    
    return success_count;
}
