#ifndef WAVE_PREPROCESSOR_H
#define WAVE_PREPROCESSOR_H

#include <fftw3.h>
#include <vector>
#include <memory>
#include <cstdint>

class WaveProcessor {
public:
   
    ~WaveProcessor();


    /**
     * 创建预处理器的工厂函数（无异常版本）
     * @param trace_length 波形长度
     * @param baseline_window_scale 滤波器缩放参数
     * @return 成功返回预处理器指针，失败返回nullptr
     */
    static WaveProcessor* create(int trace_length, double baseline_window_scale);
    
    /**
     * 智能指针版本的工厂函数
     */
    static std::unique_ptr<WaveProcessor> createUnique(int trace_length, 
                                                                    double baseline_window_scale);
    
    /**
     * 处理单个波形 - 完全线程安全
     * @param input 输入波形数据指针
     * @param output 输出缓冲区，必须预先分配至少trace_length大小
     * @return 成功返回true，失败返回false
     */
    bool process(const double* input, double* output) const;
    bool process(const double* input, double* output, double* baseline) const;

    /**
     * 处理单个波形，返回vector
     * @param input 输入波形
     * @param output 输出向量
     * @return 成功返回true，失败返回false
     */
    bool process(const double* input, std::vector<double>& output) const;
    
    /**
     * 批量处理多条波形 - 完全线程安全
     * @param inputs 输入波形数组
     * @param outputs 输出波形数组
     * @param num_traces 波形数量
     * @return 成功处理的波形数量
     */
    int processBatch(const double** inputs, double** outputs, int num_traces) const;
    
    int getTraceLength() const { return trace_length_; }
    double getWindowScale() const { return baseline_window_scale_; }
    // GET stores int16 samples: replacement means are truncated before the FFT.
    void setIntegerInput(bool value) { integer_input_ = value; }
    
    // 禁止拷贝
    WaveProcessor(const WaveProcessor&) = delete;
    WaveProcessor& operator=(const WaveProcessor&) = delete;
    
private:
    // 私有构造函数，只能通过工厂函数创建
    WaveProcessor(int trace_length, double baseline_window_scale);
    
    // 初始化函数
    bool initialize();
    
    // 内部处理函数
    bool processInternal(const double* input, double* output, 
                        double* time_buf, fftw_complex* freq_buf) const;
    
    // 计算均值和标准差
    bool computeStats(const double* data, int n, double& mean, double& sigma) const;
    
    // 创建sinc滤波器
    void createSincFilter();
    
    int trace_length_;
    double baseline_window_scale_;
    int half_length_;
    bool integer_input_ = false;
    
    // 只读的滤波器数据（线程安全）
    double* filter_time_;
    fftw_complex* filter_freq_;
    fftw_plan filter_plan_;
    
    // FFTW计划（线程安全，可重用）
    fftw_plan forward_plan_;
    fftw_plan inverse_plan_;
};

#endif 
