#include "Config.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <stdexcept>

namespace rk4 {
namespace {

std::string Lower(std::string text) {
    for (char& c : text) {
        if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
    }
    return text;
}

std::string Trim(const std::string& text) {
    const size_t begin = text.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return "";
    const size_t end = text.find_last_not_of(" \t\r\n");
    return text.substr(begin, end - begin + 1);
}

std::string StripComment(const std::string& text){
    char quote=0;bool escaped=false;
    for(size_t i=0;i<text.size();++i){
        const char c=text[i];
        if(escaped){escaped=false;continue;}
        if(quote&&c=='\\'){escaped=true;continue;}
        if(c=='\''||c=='"'){if(!quote)quote=c;else if(quote==c)quote=0;continue;}
        if(!quote&&(c=='#'||c==';')&&(i==0||text[i-1]==' '||text[i-1]=='\t'))return text.substr(0,i);
    }
    return text;
}
std::string Unquote(const std::string& value){
    if(value.size()>=2&&(value.front()=='\''||value.front()=='"')&&value.back()==value.front())
        return value.substr(1,value.size()-2);
    return value;
}

} // namespace

bool Config::Load(const std::string& path) {
    std::ifstream file(path.c_str());
    if (!file.is_open()) return false;

    values_.clear();
    std::string section;      // INI 的 [section]
    std::vector<std::pair<int, std::string>> yaml_scope; // (缩进, 前缀)
    std::string line;

    while (std::getline(file, line)) {
        // 记下缩进（YAML 用）后再 trim
        int indent = 0;
        for (char c : line) {
            if (c == ' ') {
                ++indent;
            } else if (c == '\t') {
                indent += 2;
            } else {
                break;
            }
        }
        line = Trim(StripComment(line));
        if (line.empty()) continue;
        if (line[0] == '#' || line[0] == ';') continue;

        // ---- INI 的段头 ----
        if (line.front() == '[' && line.back() == ']') {
            section = Lower(Trim(line.substr(1, line.size() - 2)));
            yaml_scope.clear();
            continue;
        }

        const size_t eq = line.find('=');
        const size_t colon = line.find(':');
        const bool is_ini = (eq != std::string::npos && (colon == std::string::npos || eq < colon));

        if (is_ini) {
            // ---- key = value ----
            std::string key = Lower(Trim(line.substr(0, eq)));
            const std::string value = Unquote(Trim(line.substr(eq + 1)));
            if (key.empty()) continue;
            if (!section.empty()) key = section + "." + key;
            values_[key] = value;
            continue;
        }

        if (colon == std::string::npos) continue;

        const std::string key = Trim(line.substr(0, colon));
        std::string value = Trim(line.substr(colon + 1));
        if (key.empty()) continue;

        if (value.empty()) {
            // ---- YAML 段头（按缩进维护嵌套前缀）----
            while (!yaml_scope.empty() && yaml_scope.back().first >= indent) yaml_scope.pop_back();
            std::string prefix=yaml_scope.empty()?"":yaml_scope.back().second;
            if (!prefix.empty()) prefix += ".";
            yaml_scope.emplace_back(indent, prefix + Lower(key));
            continue;
        }

        // ---- YAML 键值 ----
        while (!yaml_scope.empty() && yaml_scope.back().first >= indent) yaml_scope.pop_back();
        std::string prefix=yaml_scope.empty()?"":yaml_scope.back().second;
        if (!prefix.empty()) prefix += ".";
        values_[prefix + Lower(key)] = Unquote(value);
    }

    // ---- Spyral 键名 -> 本项目键名（尺寸/频率做换算）----
    struct Alias {
        const char* from;   // 小写比较
        const char* to;
        double scale;       // 值乘这个系数
    };
    static const Alias kAliases[] = {
        // Detector（Spyral：mm / V/m / T）
        {"detector.magnetic_field", "detector.bfield_t", 1.0},
        {"detector.electric_field", "detector.efield_vm", 1.0},
        {"detector.detector_length", "detector.length_m", 1e-3}, // mm -> m
        {"detector.micromegas_time_bucket", "detector.micromegas_edge", 1.0},
        {"detector.window_time_bucket", "detector.windows_edge", 1.0},
        {"detector.get_frequency", "detector.clock_mhz", 1.0},
        // Get（MHz -> ns；peak 参数同名映射）
        {"get.sampling_rate", "detector.clock_mhz", 1.0},
        {"get.peaking_time", "peaks.peaking_ns", 1.0},
        {"get.baseline_window_scale", "peaks.baseline_scale", 1.0},
        {"get.peak_separation", "peaks.separation", 1.0},
        {"get.peak_prominence", "peaks.prominence", 1.0},
        {"get.peak_max_width", "peaks.max_width", 1.0},
        {"get.peak_threshold", "peaks.threshold", 1.0},
        // Cluster
        {"cluster.radius", "cluster.eps_mm", 1.0},
        {"cluster.eps", "cluster.eps_mm", 1.0},
        {"cluster.min_points", "cluster.min_points", 1.0},
        {"cluster.min_cloud_size", "cluster.min_cloud_size", 1.0},
        {"cluster.min_cluster_size", "cluster.min_cluster_size", 1.0},
        {"cluster.direction_threshold", "cluster.direction_fraction", 1.0},
        {"estimate.smoothing_factor", "estimate.reinsch_lambda", 1.0},
        // Solver（网格参数）
        {"solver.interp_ke_min", "solver.ke_min", 1.0},
        {"solver.interp_ke_max", "solver.ke_max", 1.0},
        {"solver.interp_ke_bins", "solver.ke_bins", 1.0},
        {"solver.interp_polar_min", "solver.polar_min", 1.0},
        {"solver.interp_polar_max", "solver.polar_max", 1.0},
        {"solver.interp_polar_bins", "solver.polar_bins", 1.0},
        {"solver.n_time_steps", "solver.n_time", 1.0},
    };

    for (const Alias& a : kAliases) {
        const std::string from = a.from;
        auto it = values_.find(from);
        if (it == values_.end()) continue;
        std::string value = it->second;
        if (a.scale != 1.0) {
            char* end = nullptr;
            const double v = std::strtod(value.c_str(), &end);
            if (end != value.c_str()) {
                char buf[64];
                std::snprintf(buf, sizeof(buf), "%.10g", v * a.scale);
                value = buf;
            }
        }
        values_[a.to] = value; // 别名不删原文，方便排查
    }
    // 采样率（MHz）-> 采样周期（ns），单独算一次
    {
        auto it = values_.find("detector.clock_mhz");
        if (it != values_.end()) {
            char* end = nullptr;
            const double mhz = std::strtod(it->second.c_str(), &end);
            if (end != it->second.c_str() && mhz > 0.0) {
                char buf[64];
                std::snprintf(buf, sizeof(buf), "%.10g", 1000.0 / mhz);
                values_["peaks.sampling_ns"] = buf;
            }
        }
    }

    return true;
}

bool Config::GetBool(const std::string& key,bool def) const {
    const auto it=values_.find(key);if(it==values_.end())return def;
    const std::string value=Lower(Trim(it->second));
    if(value=="true"||value=="yes"||value=="on"||value=="1")return true;
    if(value=="false"||value=="no"||value=="off"||value=="0")return false;
    throw std::invalid_argument("Invalid boolean configuration: "+key+"="+value);
}

bool Config::Has(const std::string& key) const {
    return values_.find(key) != values_.end();
}

std::string Config::GetString(const std::string& key, const std::string& def) const {
    const auto it = values_.find(key);
    return it == values_.end() ? def : it->second;
}

double Config::GetDouble(const std::string& key, double def) const {
    const auto it = values_.find(key);
    if (it == values_.end()) return def;
    char* end = nullptr;
    const double value = std::strtod(it->second.c_str(), &end);
    return end == it->second.c_str() ? def : value;
}

int Config::GetInt(const std::string& key, int def) const {
    const auto it = values_.find(key);
    if (it == values_.end()) return def;
    char* end = nullptr;
    const long value = std::strtol(it->second.c_str(), &end, 10);
    return end == it->second.c_str() ? def : static_cast<int>(value);
}

} // namespace rk4
