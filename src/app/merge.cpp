//----------------------------------------------------------------------
// rk4_merge：合并多个拟合结果文件（多进程分片后的收尾步骤）
//
// 用法：rk4_merge <out.root> <in1.root> [in2.root ...]
//----------------------------------------------------------------------

#include "io/ResultIO.h"

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

using namespace rk4;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::printf("用法: rk4_merge <out.root> <in1.root> [in2.root ...]\n");
        return 1;
    }

    std::vector<FitRecord> merged;
    for (int i = 2; i < argc; ++i) {
        std::vector<FitRecord> part;
        if (!ReadFitResults(argv[i], part)) {
            std::fprintf(stderr, "无法读取 %s（跳过）\n", argv[i]);
            continue;
        }
        std::printf("  %s: %zu 条\n", argv[i], part.size());
        merged.insert(merged.end(), part.begin(), part.end());
    }

    // 按 (event, cluster) 排序并去重（分片之间不应重叠，但保险起见）
    std::sort(merged.begin(), merged.end(), [](const FitRecord& a, const FitRecord& b) {
        if (a.event != b.event) return a.event < b.event;
        return a.cluster < b.cluster;
    });
    merged.erase(std::unique(merged.begin(), merged.end(),
                             [](const FitRecord& a, const FitRecord& b) {
                                 return a.event == b.event && a.cluster == b.cluster;
                             }),
                 merged.end());

    if (!WriteFitResults(argv[1], merged)) {
        std::fprintf(stderr, "写出失败\n");
        return 1;
    }

    std::printf("合并完成：%zu 条 -> %s\n", merged.size(), argv[1]);
    return 0;
}
