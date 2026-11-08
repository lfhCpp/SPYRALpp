#include "MassData.h"

#include "DataPath.h"

#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>

//----------------------------------------------------------------------
// 静态成员定义
//----------------------------------------------------------------------
std::map<std::pair<int, int>, double> MassData::massTable;
std::map<std::pair<int, int>, std::string> MassData::massNameTable;
bool MassData::loaded = false;
std::string MassData::tablePath = "";

//----------------------------------------------------------------------
// 按优先级搜索质量表文件
//----------------------------------------------------------------------
bool MassData::FindTable(std::string& path) {
    // 先问 DataPath（环境变量 / 当前目录 / 安装后的 share 目录），
    // 再退回旧的相对路径写法，保证在源码树里跑宏也能工作
    const std::string candidates[] = {
        rk4::DataFile("amdc_2020"),
        "data/amdc_2020",
        "amdc_2020",
        "../data/amdc_2020",
        "../../data/amdc_2020",
    };

    for (const std::string& candidate : candidates) {
        std::ifstream inf(candidate);
        if (inf.is_open()) {
            path = candidate;
            return true;
        }
    }
    return false;
}

//----------------------------------------------------------------------
// 读入质量表（只执行一次）
//----------------------------------------------------------------------
bool MassData::Load() {
    // 多线程第一次同时用到质量表时，不能并发往静态 map 里写（会踩堆）。
    // 用函数内静态锁把整个加载过程串行化。
    static std::mutex load_mutex;
    std::lock_guard<std::mutex> lock(load_mutex);
    if (loaded) return true;

    std::string path;
    if (!FindTable(path)) {
        std::cerr << "[MassData] 找不到质量表 amdc_2020（已尝试 data/、当前目录、上级目录）"
                  << std::endl;
        return false;
    }
    tablePath = path;

    std::ifstream inf(path.c_str());
    if (!inf.is_open()) {
        std::cerr << "[MassData] 无法打开质量表: " << path << std::endl;
        return false;
    }

    // 第一行是表头 "Z A EL ATOMIC MASS(u)"
    std::string header;
    std::getline(inf, header);

    // 逐行读入，直到文件结束（原实现只读 2250 行，会漏掉重核素）
    int z = 0, a = 0;
    std::string name;
    double mass = 0.0;
    int count = 0;
    while (inf >> z >> a >> name >> mass) {
        massTable[std::make_pair(a, z)] = mass;  // 注意 Map 的键是 (A, Z)
        massNameTable[std::make_pair(a, z)] = name;
        ++count;
    }

    if (count == 0) {
        std::cerr << "[MassData] 质量表解析失败: " << path << std::endl;
        return false;
    }

    loaded = true;
    std::cout << "[MassData] 已加载 " << count << " 个核素 (" << path << ")" << std::endl;
    return true;
}

bool MassData::IsLoaded() { return Load(); }

std::string MassData::GetTablePath() {
    (void)Load();
    return tablePath;
}

//----------------------------------------------------------------------
// 查询
//----------------------------------------------------------------------
double MassData::GetNuclearMass(int A, int Z) {
    if (!Load()) return -1.0;
    const auto it = massTable.find(std::make_pair(A, Z));
    if (it == massTable.end()) {
        std::cerr << "[MassData] 质量表中没有 (Z=" << Z << ", A=" << A << ")" << std::endl;
        return -1.0;
    }
    return it->second;
}

double MassData::GetNuclearMassNuclear(int A, int Z) {
    const double atomic = GetNuclearMass(A, Z);
    if (atomic < 0.0) return atomic;
    return atomic - Z * Units::me_u;  // 扣除 Z 个电子
}

std::string MassData::GetNuclearName(int A, int Z) {
    if (!Load()) return "Unknown";
    const auto it = massNameTable.find(std::make_pair(A, Z));
    if (it == massNameTable.end()) return "Unknown";
    return it->second;
}

bool MassData::ValidIsotope(int Z, int A) {
    if (Z <= 0 || A <= 0) return false;
    if (A < Z) return false;
    return true;
}

//----------------------------------------------------------------------
// 化合物
//----------------------------------------------------------------------
double MassData::CompoundMass(const std::vector<Element>& comp) {
    double total = 0.0;
    for (const auto& elem : comp) {
        const double m = GetNuclearMass(elem.a, elem.z);
        if (m < 0.0) return -1.0;
        total += m * elem.s;
    }
    return total;
}

double MassData::CompoundAvgMass(const std::vector<Element>& comp) {
    double totalMass = 0.0;
    double totalAtoms = 0.0;
    for (const auto& elem : comp) {
        const double m = GetNuclearMass(elem.a, elem.z);
        if (m < 0.0) return -1.0;
        totalMass += m * elem.s;
        totalAtoms += elem.s;
    }
    return totalAtoms > 0 ? totalMass / totalAtoms : 0.0;
}
