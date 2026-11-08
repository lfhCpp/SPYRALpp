R__ADD_INCLUDE_PATH(test)
R__ADD_INCLUDE_PATH(extern/eigen)
#include "reinsch_lib.h"

// 小例子单元测试：n=12 的合成数据，dump P 与 f 供 Python 对拍
void reinsch_unit() {
    const int n = 12;
    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = 0.5 * i + 0.1 * i * i * 0.05;           // 非均匀
        y[i] = std::sin(0.4 * i) * 10.0 + 0.5 * i;      // 带结构
    }
    std::vector<double> f;
    if (!ReinschSmooth(x, y, 10.0, f)) { std::printf("solve failed\n"); return; }
    FILE* fp = std::fopen("/tmp/unit_cpp.csv", "w");
    std::fprintf(fp, "x,y,f\n");
    for (int i = 0; i < n; ++i) std::fprintf(fp, "%.9f,%.9f,%.9f\n", x[i], y[i], f[i]);
    std::fclose(fp);
    std::printf("unit test: n=%d f[0]=%.4f f[n-1]=%.4f\n", n, f[0], f[n - 1]);
}
