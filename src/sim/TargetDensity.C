#include "TargetDensity.h"
#include "../core/MassData.h"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>

struct Element;

//--------------------->Target 基类<-----------------------//

// Target 基类实现
double Target::calcMolarMass(const std::vector<Element>& compound) {
    return MassData::compoundMass(compound);
}

// GasTarget 实现
GasTarget::GasTarget(const std::vector<Element>& Comp, double pressure)
    : comp(Comp), pressure(pressure) {
    
    if (!Utils::validCompound(Comp)) {
        throw TargetError("无效化合物");
    }
    
    if (!std::isfinite(pressure) || pressure <= 0) {
        throw TargetError("压力需为正数");
    }
    
    molarMassVal = calcMolarMass(Comp);
    
    calcDensity();
}

void GasTarget::calcDensity() {
    // 使用精确原子质量计算密度
    massDensity = Utils::gasDensity(molarMassVal, pressure, Constants::roomTemp);
}

double GasTarget::numberDensity() const {
    if (molarMassVal <= 0) return 0.0;
    return (massDensity * Constants::avogadro) / molarMassVal;
}

std::string GasTarget::str() const {
    std::stringstream ss;
    ss << "Gas[";
    for (const auto& elem : comp) {
        ss << MassData::GetNuclearName(elem.a, elem.z)  << elem.s;
        if (&elem != &comp.back()) ss << " ";
    }
    ss << "] P=" << pressure << " Torr";
    return ss.str();
}

// // SolidTarget 实现
SolidTarget::SolidTarget(const std::vector<Element>& Comp, double thickness)
    : comp(Comp), thickness(thickness) {
    
    if (!Utils::validCompound(Comp)) {
        throw TargetError("无效化合物");
    }
    
    if (!std::isfinite(thickness) || thickness <= 0) {
        throw TargetError("厚度需为正数");
    }
    
    molarMassVal = MassData::compoundAvgMass(Comp);
    calcArealDensity();
}

void SolidTarget::calcArealDensity() {
    arealDens = Utils::arealDensity(thickness);
}

double SolidTarget::thicknessCm(double massDensity) const {
    
    if (massDensity <= 0) {
        throw TargetError("质量密度需为正数");
    }
    return arealDens / massDensity;

}

std::string SolidTarget::str() const {
   
    std::stringstream ss;
    ss << "Solid[";
    for (const auto& elem : comp) {
        ss << MassData::GetNuclearName(elem.a, elem.z)  << elem.s;
        if (&elem != &comp.back()) ss << " ";
    }
    ss << "] T=" << thickness << " μg/cm²";
    return ss.str();

}

// GasMixTarget 实现
GasMixTarget::GasMixTarget(const std::vector<std::vector<Element>>& comps,
                         const std::vector<double>& fractions,
                         double pressure)
    : pressure(pressure) {
    
    if (comps.size() != fractions.size()) {
        throw TargetError("组分数量不匹配");
    }
    
    if (comps.empty()) {
        throw TargetError("至少需要一个组分");
    }
    
    double total = 0.0;
    for (size_t i = 0; i < comps.size(); ++i) {
        if (!Utils::validCompound(comps[i])) {
            throw TargetError("组分" + std::to_string(i) + "无效");
        }
        
        if (!std::isfinite(fractions[i]) || fractions[i] < 0.0 || fractions[i] > 1.0) {
            throw TargetError("分数需在0-1之间");
        }
        
        components.push_back({comps[i], fractions[i]});
        total += fractions[i];
    }
    
    if (std::abs(total - 1.0) > 1e-6) {
        throw TargetError("分数之和需为1.0");
    }
    
    if (!std::isfinite(pressure) || pressure <= 0) {
        throw TargetError("压力需为正数");
    }
    
    for (auto& component : components) component.fraction /= total;
    calcMixProps();
}

void GasMixTarget::calcMixProps() {
   
    // 计算平均摩尔质量
    molarMassVal = 0.0;
    for (const auto& comp : components) {
        double compMass = compMolarMass(comp.comp);
        molarMassVal += compMass * comp.fraction;
    }
    
    // 计算密度
    massDensity = Utils::gasDensity(molarMassVal, pressure, Constants::roomTemp);
    
    // 构建等效化合物
    std::map<std::pair<int,int>,double> counts;
    for (const auto& comp : components) {
        for (const auto& elem : comp.comp) {
            counts[{elem.a,elem.z}] += elem.s * comp.fraction;
        }
    }
    equivComp.clear();
    for (const auto& entry : counts) if (entry.second > 0)
        equivComp.push_back({entry.first.first,entry.first.second,entry.second});
    
}

double GasMixTarget::compMolarMass(const std::vector<Element>& comp) const {
    return Target::calcMolarMass(comp);
}

double GasMixTarget::numberDensity() const {
    if (molarMassVal <= 0) return 0.0;
    return (massDensity * Constants::avogadro) / molarMassVal;
}

double GasMixTarget::partialPressure(int idx) const {
    if (idx < 0 || idx >= static_cast<int>(components.size())) {
        throw TargetError("无效索引");
    }
    return pressure * components[idx].fraction;
}

double GasMixTarget::compNumberDensity(int idx) const {
    double pp = partialPressure(idx);
    return Utils::numberDensity(pp, Constants::roomTemp);
}

std::string GasMixTarget::str() const {
    std::stringstream ss;
    ss << "GasMix[";
    for (size_t i = 0; i < components.size(); ++i) {
        ss << std::fixed << std::setprecision(0) 
           << components[i].fraction * 100.0 << "% ";
        for (const auto& elem : components[i].comp) {
            ss << MassData::GetNuclearName(elem.a, elem.z)  << elem.s;
        }
        if (i < components.size() - 1) ss << " + ";
    }
    ss << "] P=" << pressure << " Torr";
    return ss.str();
}

// Utils 实现
namespace Utils {
    double stoichiometricScale(const std::vector<Element>& comp) {
        if (!validCompound(comp)) throw TargetError("Invalid atomic composition");
        double smallest = 1.0;
        for (const auto& elem : comp) smallest = std::min(smallest, elem.s);
        const double scale = 1.0 / smallest;
        for (const auto& elem : comp)
            if (!std::isfinite(elem.s * scale)) throw TargetError("Atomic count ratio overflows");
        return scale;
    }
    bool validCompound(const std::vector<Element>& comp) {
        if (comp.empty()) return false;
        for (const auto& elem : comp) {
            if (!MassData::validIsotope(elem.z, elem.a)) return false;
            if (!std::isfinite(elem.s) || elem.s <= 0) return false;
        }
        return true;
    }
    
    double gasDensity(double molarMass, double pressure, double temp) {
        // PV = nRT, n = m/M, m = rho * V, so rho = (P * M) / (R * T)
        return (pressure * molarMass) / (Constants::gasConstant * temp);
    }
    
    double numberDensity(double pressure, double temp) {
        // n = P / (kT)
        return (pressure * 133.322) / (Constants::boltzmann * temp) * 1e-6;  // Torr to Pa, m³ to cm³
    }
    
    double arealDensity(double thickness) {
        return thickness * Constants::ug2g;  // μg/cm² to g/cm²
    }
}
