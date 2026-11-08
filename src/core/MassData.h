#ifndef MASSDATA_H
#define MASSDATA_H

//----------------------------------------------------------------------
// 核素质量表（数据来自 AMDC 2020，文件 data/amdc_2020）
//
// 设计要点：
//   1. 文件只读一次，结果缓存在静态 map 里；
//      旧实现在 GetNuclearMass() 里每次重读 2550 行，是主要性能瓶颈之一。
//   2. 表里的质量是"原子质量"（含 Z 个电子）。做反应运动学时
//      需要的是原子核质量，用 GetNuclearMassNuclear() 减掉 Z*me。
//   3. 路径自动搜索，宏在 RK4/ 或 RK4/test/ 下运行都能找到数据文件。
//----------------------------------------------------------------------

#include "Element.h"
#include "Units.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

class MassData {
public:
    // 原子质量（含电子），单位 u
    static double GetNuclearMass(int A, int Z);

    // 原子核质量（原子质量 - Z*me），单位 u
    static double GetNuclearMassNuclear(int A, int Z);

    // 核素名称，例如 "12C"；找不到返回 "Unknown"
    static std::string GetNuclearName(int A, int Z);

    static bool ValidIsotope(int Z, int A);

    // 质量表是否成功加载 / 实际使用的文件路径（调试用）
    static bool IsLoaded();
    static std::string GetTablePath();

    // 化合物的总质量 / 平均质量，单位 u
    static double CompoundMass(const std::vector<Element>& comp);
    static double CompoundAvgMass(const std::vector<Element>& comp);

    // ---- 兼容旧接口（原 MassData 的命名）----
    static void initMassTable() { (void)Load(); }
    static bool validIsotope(int z, int a) { return ValidIsotope(z, a); }
    static double compoundMass(const std::vector<Element>& comp) { return CompoundMass(comp); }
    static double compoundAvgMass(const std::vector<Element>& comp) { return CompoundAvgMass(comp); }

private:
    static bool Load();
    static bool FindTable(std::string& path);

    static std::map<std::pair<int, int>, double> massTable;
    static std::map<std::pair<int, int>, std::string> massNameTable;
    static bool loaded;
    static std::string tablePath;
};

#endif // MASSDATA_H
