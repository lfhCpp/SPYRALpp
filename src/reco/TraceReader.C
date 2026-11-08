#include "TraceReader.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace rk4 {

namespace {

bool IsEventHeader(const std::string& line, int& number) {
    // "========== event 0 =========="
    if (line.find("==========") == std::string::npos) return false;
    const size_t pos = line.find("event");
    if (pos == std::string::npos) return false;
    const char* start = line.c_str() + pos + 5;
    char* end = nullptr;
    const long value = std::strtol(start, &end, 10);
    if (end == start) return false;
    number = static_cast<int>(value);
    return true;
}

// "pads 12: 3 0 1 17 5632 <512 个采样>"
bool ParsePadLine(const std::string& line, Trace& trace) {
    const size_t colon = line.find(':');
    if (colon == std::string::npos) return false;

    // pad 号
    const size_t pads_pos = line.find("pads");
    if (pads_pos == std::string::npos) return false;
    const char* p = line.c_str() + pads_pos + 4;
    char* end = nullptr;
    const long pad = std::strtol(p, &end, 10);
    if (end == p) return false;
    trace.pad = static_cast<int>(pad);

    std::istringstream ss(line.substr(colon + 1));
    double value = 0.0;

    // 前 5 个数是元数据（长度、触发信息、基线等），跳过
    for (int i = 0; i < 5; ++i) {
        if (!(ss >> value)) return false;
    }

    trace.samples.clear();
    while (ss >> value) {
        trace.samples.push_back(static_cast<float>(value));
    }

    return !trace.samples.empty();
}

} // namespace

bool TextTraceReader::Open(const std::string& path) {
    std::ifstream file(path.c_str());
    if (!file.is_open()) {
        std::fprintf(stderr, "[TextTraceReader] 无法打开 %s\n", path.c_str());
        return false;
    }

    lines_.clear();
    event_starts_.clear();
    cursor_ = 0;

    std::string line;
    while (std::getline(file, line)) {
        int number = -1;
        if (IsEventHeader(line, number)) {
            event_starts_.push_back(lines_.size());
        }
        lines_.push_back(line);
    }

    if (event_starts_.empty()) {
        std::fprintf(stderr, "[TextTraceReader] 没找到事件头：%s\n", path.c_str());
        return false;
    }
    return true;
}

bool TextTraceReader::ReadNext(TraceEvent& ev) {
    if (cursor_ >= event_starts_.size()) return false;

    const size_t begin = event_starts_[cursor_];
    const size_t end = (cursor_ + 1 < event_starts_.size()) ? event_starts_[cursor_ + 1]
                                                            : lines_.size();
    ++cursor_;

    ev.traces.clear();
    int number = -1;
    if (!IsEventHeader(lines_[begin], number)) return false;
    ev.number = number;

    for (size_t i = begin + 1; i < end; ++i) {
        const std::string& line = lines_[i];
        if (line.compare(0, 5, "pads ") != 0) continue;
        Trace trace;
        if (ParsePadLine(line, trace)) {
            ev.traces.push_back(trace);
        }
    }

    return true;
}

} // namespace rk4
