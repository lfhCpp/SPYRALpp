#include "src/sim/TargetDensity.h"
#include "src/core/MassData.h"
#include <iostream>
#include <iomanip>

struct Element;

void testGasTarget() {

    //原子质量（以u为单位）在数值上等于摩尔质量（以g/mol为单位）
    std::cout << "\n=== 测试气体靶 ===" << std::endl;
    
    
    // H₂
    std::vector<Element> h2 = {{1, 1, 2}};
    GasTarget h2Target(h2, 760.0);
    
    std::cout << h2Target.str() << std::endl;
    std::cout << "质量: " << h2Target.molarMass() << " u" << std::endl; 
    std::cout << "密度: " << h2Target.density() << " g/cm³" << std::endl;
    std::cout << "分子数密度: " << h2Target.numberDensity() << " /cm³" << std::endl;
    
    // D₂
    std::vector<Element> d2 = {{2, 1, 2}};
    GasTarget d2Target(d2, 600.0);
    
    std::cout << "\n" << d2Target.str() << std::endl;
    std::cout << "质量: " << d2Target.molarMass() << " u" << std::endl;
    std::cout << "密度: " << d2Target.density() << " g/cm³" << std::endl;
    
    // CO₂
    std::vector<Element> co2 = {{12, 6, 1}, {16, 8, 2}};
    GasTarget co2Target(co2, 760.0);
    
    std::cout << "\n" << co2Target.str() << std::endl;
    std::cout << "质量: " << co2Target.molarMass() << " u" << std::endl;
    std::cout << "密度: " << co2Target.density() << " g/cm³" << std::endl;
 }

void testSolidTarget() {
    std::cout << "\n=== 测试固体靶 ===" << std::endl;
    
    // 金箔
    std::vector<Element> au = {{197, 79, 1}};
    SolidTarget auTarget(au, 100.0);
    
    std::cout << auTarget.str() << std::endl;
    std::cout << "面密度: " << auTarget.arealDensity() << " g/cm²" << std::endl;
    std::cout << "厚度(密度19.3g/cm³时): " 
              << auTarget.thicknessCm(19.3) * 1e4 << " μm" << std::endl;
    
    // 聚乙烯 (C₂H₄)ₙ
    std::vector<Element> pe = {{12, 6, 2}, {1, 1, 4}};
    SolidTarget peTarget(pe, 50.0);
    
    std::cout << "\n" << peTarget.str() << std::endl;
    std::cout << "质量: " << peTarget.molarMass() << " u" << std::endl;
    std::cout << "面密度: " << peTarget.arealDensity() << " g/cm²" << std::endl;
}

void testGasMixTarget() {
    std::cout << "\n=== 测试混合气体靶 ===" << std::endl;
    
    // 90% Ar  + 10% CH4
    std::vector<std::vector<Element>> comps = {
        {{40, 18, 1}},  // Ar
        {{12, 6, 1}, {1, 1, 4}}   // CH4
    };
    std::vector<double> fracs = {0.9, 0.1};
    
    GasMixTarget mixTarget(comps, fracs, 600.0);
    
    std::cout << mixTarget.str() << std::endl;
    std::cout << "平均质量: " << mixTarget.molarMass() << " u" << std::endl;
    std::cout << "密度: " << mixTarget.density() << " g/cm³" << std::endl;
    std::cout << "分子数密度: " << mixTarget.numberDensity() << " /cm³" << std::endl;
    
    // 组分信息
    for (size_t i = 0; i < mixTarget.comps().size(); ++i) {
        std::cout << "组分" << i << "分压: " << mixTarget.partialPressure(i) 
                  << " Torr" << std::endl;
    }
}

void testMassData() {
    std::cout << "\n=== 测试原子质量数据库 ===" << std::endl;
    
    std::cout << "H-1质量: " << MassData::GetNuclearMass(1, 1) << " u" << std::endl;
    std::cout << "C-12质量: " << MassData::GetNuclearMass(12, 6) << " u" << std::endl;
    std::cout << "O-16质量: " << MassData::GetNuclearMass(16, 8) << " u" << std::endl;
    std::cout << "Au-197质量: " << MassData::GetNuclearMass(197, 79) << " u" << std::endl;
    
    
    // 测试化合物质量
    std::vector<Element> h2o = {{1, 1, 2}, {16, 8, 1}};
    std::cout << "H₂O摩尔质量: " << MassData::compoundMass(h2o) << " u" << std::endl;
}

int main() {
    std::cout << std::fixed << std::setprecision(8);
    
    try {
        // 测试原子质量数据库
        testMassData();
        
        // 测试各种靶材
        testGasTarget();
        testSolidTarget();
        testGasMixTarget();
        
        // 测试工具函数
        std::cout << "\n=== 测试工具函数 ===" << std::endl;
        double h2Density = Utils::gasDensity(2.016, 760.0);
        std::cout << "H₂理论密度: " << h2Density << " g/cm³" << std::endl;
        
        double nDensity = Utils::numberDensity(760.0);
        std::cout << "标准状态分子数密度: " << nDensity << " /cm³" << std::endl;
        
    } catch (const TargetError& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "异常: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n=== 所有测试通过 ===" << std::endl;
    return 0;
}

