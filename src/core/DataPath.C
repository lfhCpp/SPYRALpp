#include "DataPath.h"

#include <cstdlib>
#include <fstream>
#include <vector>

#if defined(__linux__)
#include <limits.h>
#include <unistd.h>
#endif

namespace rk4 {

namespace {

bool FileExists(const std::string& path) {
    std::ifstream inf(path.c_str());
    return inf.is_open();
}

// 与可执行文件位置相关的候选目录：
//   安装后： <prefix>/bin/xxx        -> <prefix>/share/spyralpp/data
//   构建树： <repo>/build/rk4_xxx    -> <repo>/data
std::vector<std::string> ExeRelativeDataDirs() {
    std::vector<std::string> out;
#if defined(__linux__)
    char buf[PATH_MAX];
    const ssize_t len = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (len <= 0) return out;
    buf[len] = '\0';
    const std::string exe(buf);
    const std::size_t slash = exe.find_last_of('/');
    if (slash == std::string::npos) return out;
    const std::string dir = exe.substr(0, slash);
    out.push_back(dir + "/../share/spyralpp/data");  // 安装后
    out.push_back(dir + "/../data");                  // 构建树：<repo>/build/xxx
    out.push_back(dir + "/../../data");               // 构建树：<repo>/build/bin/xxx
#endif
    return out;
}

} // namespace

const std::string& DataDir() {
    static const std::string cached = []() {
        std::vector<std::string> candidates;
        if (const char* env = std::getenv("SPYRALPP_DATA_DIR")) {
            if (env[0] != '\0') candidates.push_back(env);
        }
        candidates.push_back("data");
        candidates.push_back(".");
        candidates.push_back("../data");
        candidates.push_back("../../data");
        for (const std::string& exe_rel : ExeRelativeDataDirs()) candidates.push_back(exe_rel);

        for (const std::string& dir : candidates) {
            if (FileExists(dir + "/amdc_2020")) return dir;
        }
        return std::string("data");
    }();
    return cached;
}

std::string DataFile(const std::string& name) { return DataDir() + "/" + name; }

} // namespace rk4
