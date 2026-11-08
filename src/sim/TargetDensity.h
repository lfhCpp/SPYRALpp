#ifndef TARGET_DENSITY_HPP
#define TARGET_DENSITY_HPP

//----------------------------------------------------------------------
// 靶材：气体靶、固体靶、混合气体靶
//
// 依赖：core/Element.h（元素项）、core/Units.h（常数与换算）
// 说明：Constants 命名空间保留，内部全部指向 Units，避免两套常数。
//----------------------------------------------------------------------

#include "../core/Element.h"
#include "../core/Units.h"

#include <vector>
#include <string>
#include <memory>
#include <map>
#include <utility>
#include <stdexcept>


// 物理常数（兼容旧代码的 Constants::xxx 写法，实际定义在 core/Units.h）
namespace Constants {
    constexpr double gasConstant = Units::gasConstant;    // cm³·Torr/(mol·K)
    constexpr double roomTemp    = Units::roomTemp;       // K (20°C)
    constexpr double avogadro    = Units::avogadro;       // 阿伏伽德罗常数
    constexpr double ug2g        = Units::ug2g;           // μg 转 g
    constexpr double boltzmann   = Units::boltzmann;      // J/K

    constexpr double me          = Units::me_u;           // 电子静止质量 (u)
    constexpr double u_MeV       = Units::u_MeV;          // 1u -> MeV/c^2
    constexpr double C           = Units::C;              // 光速 (m/s)
    constexpr double E_CHARGE    = Units::E_CHARGE;       // 基本电荷 (C)
    constexpr double MEV_2_KG    = Units::MEV_2_KG;       // MeV/c² to kg
    constexpr double MEV_2_JOULE = Units::MEV_2_JOULE;    // MeV to J
}

// 异常类
class TargetError : public std::runtime_error {
public:
    explicit TargetError(const std::string& msg) : std::runtime_error(msg) {}
};

// 靶材基类
class Target {
    public:
        virtual ~Target() = default;

        // 获取密度 (g/cm³)
        virtual double density() const = 0;

        // 获取分子数密度 (molecules/cm³)
        virtual double numberDensity() const = 0;

        // 获取面密度 (g/cm²)
        virtual double arealDensity() const = 0;

        // 获取摩尔质量 (g/mol)
        virtual double molarMass() const = 0;

        // 获取化合物
        virtual std::vector<Element> compound() const = 0;

        // 字符串表示
        virtual std::string str() const = 0;

    protected:
        // 计算摩尔质量g/mol=u
        static double calcMolarMass(const std::vector<Element>& compound);
};

// 气体靶
class GasTarget : public Target {

    private:
        std::vector<Element> comp;

        double pressure;    // Torr
        double massDensity; // g/cm³
        double molarMassVal; // g/mol

        void calcDensity();

    public:
        GasTarget(const std::vector<Element>& comp, double pressure);

        double density() const override { return massDensity; } //常成员函数 重写基类的虚函数
        double numberDensity() const override;
        double arealDensity() const override { return 0.0; }
        double molarMass() const override { return molarMassVal; }
        std::vector<Element> compound() const override { return comp; }

        std::string str() const override;

        double pressureTorr() const { return pressure; }
};

// // 固体靶
class SolidTarget : public Target {

    private:

        std::vector<Element> comp;
        double thickness;   // μg/cm²
        double arealDens;   // g/cm²
        double molarMassVal; // g/mol

        void calcArealDensity();

    public:

        SolidTarget(const std::vector<Element>& comp, double thickness);

        double density() const override { return 0.0; }
        double numberDensity() const override { return 0.0; }
        double arealDensity() const override { return arealDens; }
        double molarMass() const override { return molarMassVal; }
        std::vector<Element> compound() const override { return comp; }

        std::string str() const override;

        double thicknessUgCm2() const { return thickness; }
        double thicknessCm(double massDensity) const;

};

// 气体混合物靶
class GasMixTarget : public Target {

private:
    struct MixComp {
        std::vector<Element> comp;
        double fraction;  // 体积分数
    };

    std::vector<MixComp> components;
    double pressure;       // Torr
    double massDensity;    // g/cm³
    double molarMassVal;   // g/mol
    std::vector<Element> equivComp;  // Mole-fraction weighted atomic counts.

    void calcMixProps();
    double compMolarMass(const std::vector<Element>& comp) const;

public:

    GasMixTarget(const std::vector<std::vector<Element>>& components,
                const std::vector<double>& fractions,
                double pressure);

    double density() const override { return massDensity; }
    double numberDensity() const override;
    double arealDensity() const override { return 0.0; }
    double molarMass() const override { return molarMassVal; }
    std::vector<Element> compound() const override { return equivComp; }

    std::string str() const override;

    const std::vector<MixComp>& comps() const { return components; }
    double pressureTorr() const { return pressure; }

    double partialPressure(int idx) const;
    double compNumberDensity(int idx) const;
};

// 工具函数
namespace Utils {
    // Catima treats counts below 1 as mass fractions. Scale all counts together
    // so that fractional mole composition remains stoichiometric input.
    double stoichiometricScale(const std::vector<Element>& comp);

    // 验证化合物
    bool validCompound(const std::vector<Element>& comp);

    // 计算气体密度 (g/cm³)
    double gasDensity(double molarMass, double pressure, double temp = Constants::roomTemp);

    // 计算分子数密度 (molecules/cm³)
    double numberDensity(double pressure, double temp = Constants::roomTemp);

    // 计算面密度 (g/cm²)
    double arealDensity(double thickness);
}

#endif // TARGET_DENSITY_HPP
