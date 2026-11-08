//----------------------------------------------------------------------
// 数据目录解析测试：DataDir() 必须能定位到 amdc_2020
//   root -l -b -q test/testDataPath.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/DataPath.h"
#include "core/MassData.h"

#include <cstdio>
#include <fstream>
#include <string>

static int g_fail = 0;

static void Check(bool ok, const char* msg) {
    if (!ok) {
        std::printf("  [FAIL] %s\n", msg);
        ++g_fail;
    }
}

void testDataPath() {
    std::printf("========== 数据目录解析测试 ==========\n");

    const std::string dir = rk4::DataDir();
    const std::string table = rk4::DataFile("amdc_2020");
    std::printf("    DataDir()  = %s\n", dir.c_str());
    std::printf("    DataFile() = %s\n", table.c_str());

    Check(!dir.empty(), "DataDir 非空");
    std::ifstream inf(table.c_str());
    Check(inf.is_open(), "DataFile(amdc_2020) 可打开");
    Check(table == dir + "/amdc_2020", "DataFile 拼接正确");

    // 质量表能实际加载（MassData 是全局类，Load() 私有，用公开接口触发）
    Check(MassData::GetNuclearMassNuclear(12, 6) > 11.0, "12C 核质量合理");
    Check(MassData::IsLoaded(), "质量表已加载");
    Check(MassData::GetTablePath() == table, "MassData 用的是同一个文件");

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
    }
}
