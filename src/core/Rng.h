#ifndef RK4_CORE_RNG_H
#define RK4_CORE_RNG_H

//----------------------------------------------------------------------
// 可复现的随机数发生器（PCG32）
//
// 为什么不用 TRandom3：
//   1. 物理内核 core/ 不依赖 ROOT，便于单独编译与单元测试；
//   2. PCG32 整数序列可跨平台复现；高斯变换还依赖平台数学库。
// 原始 GET 的 NumPy PCG64 参考随机流见 reco/NumpyPcg64.h。
//
// 用法：
//   Rng rng(12345);
//   double u = rng.Uniform();            // [0,1)
//   double g = rng.Gauss(0.0, 1.0);      // 正态
//   int    i = rng.Integer(10);       // [0,10)
//----------------------------------------------------------------------

#include <cstdint>
#include <cmath>

namespace core {

class Rng {
public:
    explicit Rng(uint64_t seed = 1, uint64_t stream = 1) { Seed(seed, stream); }

    void Seed(uint64_t seed, uint64_t stream = 1) {
        state_ = 0u;
        inc_ = (stream << 1u) | 1u;
        NextU32();
        state_ += seed;
        NextU32();
        has_gauss_ = false;
    }

    // 32 位无符号均匀分布
    uint32_t NextU32() {
        const uint64_t old = state_;
        state_ = old * 6364136223846793005ULL + inc_;
        const uint32_t xorshifted = static_cast<uint32_t>(((old >> 18u) ^ old) >> 27u);
        const uint32_t rot = static_cast<uint32_t>(old >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31u));
    }

    // [0,1) 均匀分布
    double Uniform() { return (NextU32() >> 8) * 0x1.0p-24; }

    // [a,b) 均匀分布
    double Uniform(double a, double b) { return a + (b - a) * Uniform(); }

    // [0,n) 整数
    int Integer(int n) { return static_cast<int>(Uniform() * n); }

    // 标准正态（Marsaglia 极坐标法，成对生成，缓存第二个）
    double Gauss(double mu = 0.0, double sigma = 1.0) {
        if (has_gauss_) {
            has_gauss_ = false;
            return mu + sigma * gauss_cache_;
        }
        double u, v, s;
        do {
            u = 2.0 * Uniform() - 1.0;
            v = 2.0 * Uniform() - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        const double f = std::sqrt(-2.0 * std::log(s) / s);
        gauss_cache_ = v * f;
        has_gauss_ = true;
        return mu + sigma * u * f;
    }

    // 指数分布（平均值为 mean）
    double Exp(double mean) { return -mean * std::log(1.0 - Uniform()); }

private:
    uint64_t state_ = 0;
    uint64_t inc_ = 1;
    bool has_gauss_ = false;
    double gauss_cache_ = 0.0;
};

} // namespace core

#endif // RK4_CORE_RNG_H
