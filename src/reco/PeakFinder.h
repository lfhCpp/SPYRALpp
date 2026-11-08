#ifndef PEAK_FINDER_H
#define PEAK_FINDER_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <functional>

// ROOT 头文件
#include "TH1D.h"
#include "TSpectrum.h"
#include "TMath.h"

class Peak {
public:
    double centroid;           // 峰值质心位置
    double amplitude;         // 峰值振幅
    double integral;          // 积分电荷
    int left_ip;              // 左拐点
    int right_ip;             // 右拐点
    double prominence;        // 峰值突出度
    double width;             // 峰值宽度
    size_t dither_index = 0;   // Draw ordinal, including amplitude-rejected peaks.
    
    Peak() : centroid(0), amplitude(0), integral(0), 
             left_ip(0), right_ip(0), prominence(0), width(0) {}
};

class PeakFinderParameters {
public:
    double peak_separation = 20;       // 最小峰值间距
    double peak_prominence = 50.0;  // 最小突出度
    double peak_max_width = 100.0;  // 最大宽度
    double peak_threshold = 40.0;   // 振幅阈值
    double rel_height = 0.95;       // 宽度计算相对高度
    
    double peak_min_width = 1.0;

    // 构造函数
    PeakFinderParameters() = default;
};

class PeakFinder {
private:
    std::vector<double> trace_;      // 迹线数据
    std::vector<Peak> peaks_;        // 检测到的峰值
    PeakFinderParameters params_;    // 参数
    size_t dither_draws_ = 0;
    
    // 辅助函数
    std::vector<int> findLocalMaxima(const std::vector<double>& data);
    double calculateProminence(const std::vector<double>& data, int peak_idx, 
                             int& left_base, int& right_base);
    void calculateWidth(const std::vector<double>& data, int peak_idx, 
                       double prominence, int left_base, int right_base,
                       double& left_ip, double& right_ip);
    
public:
    // 构造函数
    PeakFinder() = default;
    explicit PeakFinder(const PeakFinderParameters& params) : params_(params) {}
    
    // 设置迹线
    void setTrace(const std::vector<double>& trace) { trace_ = trace; }
    
    // 设置参数
    void setParameters(const PeakFinderParameters& params) { params_ = params; }
    
    // 主函数：检测峰值
    std::vector<Peak> findPeaks(std::mt19937& rng);
    std::vector<Peak> findPeaksWithDither(const std::function<double()>& draw);
    size_t DitherDraws() const { return dither_draws_; }
    
    // 获取峰值数量
    int getNumberOfPeaks() const { return peaks_.size(); }
    
    // 获取峰值列表
    const std::vector<Peak>& getPeaks() const { return peaks_; }
    
    // 清空峰值
    void clearPeaks() { peaks_.clear(); }
    
    // 使用ROOT TSpectrum的备选方法
    // std::vector<Peak> findPeaksWithRoot(const std::vector<double>& trace, 
    //                                     int n_peaks = 20, 
    //                                     double sigma = 2.0,
    //                                     double threshold = 0.1);
};

#endif
