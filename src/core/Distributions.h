#ifndef RK4_CORE_DISTRIBUTIONS_H
#define RK4_CORE_DISTRIBUTIONS_H

//----------------------------------------------------------------------
// 抽样分布：激发能、极角
//
// 与 attpc_engine 对应：
//   ExcitationGaussian / ExcitationUniform / ExcitationBreitWigner
//   PolarUniform（对 cos(theta) 均匀，等价于对立体角均匀）
//
// 所有分布都只依赖 core/Rng.h，便于单元测试。
//----------------------------------------------------------------------

#include "Rng.h"

#include <cmath>

namespace rk4 {

//----------------------------------------------------------------------
// 激发能分布（MeV）
//----------------------------------------------------------------------
class ExcitationDistribution {
public:
    virtual ~ExcitationDistribution() = default;
    virtual double Sample(core::Rng& rng) const = 0;
};

class ExcitationGaussian : public ExcitationDistribution {
public:
    // centroid 单位 MeV，fwhm 为半高全宽（MeV）
    ExcitationGaussian(double centroid, double fwhm)
        : centroid_(centroid), fwhm_(fwhm), sigma_(fwhm / 2.354820045) {}

    double Sample(core::Rng& rng) const override {
        return rng.Gauss(centroid_, sigma_);
    }

    double GetCentroid() const { return centroid_; }
    double GetFwhm() const { return fwhm_; }

private:
    double centroid_;
    double fwhm_;
    double sigma_;
};

class ExcitationUniform : public ExcitationDistribution {
public:
    ExcitationUniform(double min_value, double max_value)
        : min_(min_value), max_(max_value) {}

    double Sample(core::Rng& rng) const override { return rng.Uniform(min_, max_); }

private:
    double min_;
    double max_;
};

// 非相对论 Breit-Wigner：P(E) ∝ (Γ/2)^2 / ((E-E0)^2 + (Γ/2)^2)
// 用逆变换法抽样：E = E0 + (Γ/2) tan(pi (u - 0.5))
class ExcitationBreitWigner : public ExcitationDistribution {
public:
    ExcitationBreitWigner(double centroid, double width)
        : centroid_(centroid), width_(width) {}

    double Sample(core::Rng& rng) const override {
        const double u = rng.Uniform();
        return centroid_ + 0.5 * width_ * std::tan(3.14159265358979323846 * (u - 0.5));
    }

private:
    double centroid_;
    double width_;
};

//----------------------------------------------------------------------
// 极角分布（质心系，弧度）
//----------------------------------------------------------------------
class PolarDistribution {
public:
    virtual ~PolarDistribution() = default;
    virtual double SamplePolar(core::Rng& rng) const = 0;
};

// 对 cos(theta) 均匀：立体角均匀
class PolarUniform : public PolarDistribution {
public:
    // 输入为角度（弧度），内部转换成 cos 的上下限
    PolarUniform(double angle_min, double angle_max)
        : cos_min_(std::cos(angle_max)), cos_max_(std::cos(angle_min)) {}

    double SamplePolar(core::Rng& rng) const override {
        return std::acos(rng.Uniform(cos_min_, cos_max_));
    }

private:
    double cos_min_;
    double cos_max_;
};

} // namespace rk4

#endif // RK4_CORE_DISTRIBUTIONS_H
