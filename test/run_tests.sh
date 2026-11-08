#!/bin/bash
#----------------------------------------------------------------------
# RK4 测试运行器
#
# 约定：
#   1. 先在 RK4/ 下构建共享库：  cmake -S . -B build && cmake --build build -j
#   2. 测试本身是 ROOT 宏，不需要单独编译，直接：
#        bash test/run_tests.sh
#
# 每个测试宏返回 0 表示通过；脚本汇总结果并以非零码退出（失败时）。
#----------------------------------------------------------------------

set -u

ROOTSYS="${ROOTSYS:-/home/lfh/work/root6.30}"
export PATH="$ROOTSYS/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"

# 切到 RK4 根目录（宏里的相对路径都以此为准）
cd "$(dirname "$0")/.." || exit 1

if [ ! -f build/librk4.so ]; then
    echo "缺少 build/librk4.so，请先构建："
    echo "  cmake -S . -B build && cmake --build build -j"
    exit 2
fi

TESTS=(
    "test/testKinematics.C"
    "test/testIO.C"
    "test/testDataPath.C"
    "test/testHdf5.C"
    "test/testH5Trace.C"
    "test/testH5Cloud.C"
    "test/testDetector.C"
    "test/testWaveform.C"
    "test/testCluster.C"
    "test/testSolver.C"
    "test/testKalman.C"
    "test/testMeshCache.C"
    "test/testParallel.C"
    "test/testShard.C"
    "test/testSimParallel.C"
    "test/testHdbscan.C"
    "test/testSpline.C"
    "test/testParticleId.C"
    "test/testPipeline.C"
    "test/testWaveformPipeline.C"
    "test/testRK4.C"
    "test/testPulseShape.C"
)

pass=0
fail=0
failed_names=()

for t in "${TESTS[@]}"; do
    name=$(basename "$t")
    printf '%-32s' "$name"
    log=$(mktemp)
    if timeout 600 root -l -b -q "$t" > "$log" 2>&1; then
        # Cling 的错误有时不影响退出码，额外扫一遍关键字
        if grep -qiE "error:|fatal error|Failed to materialize|unresolved" "$log"; then
            echo "FAIL"
            sed -n '1,10p' "$log" | sed 's/^/    /'
            fail=$((fail + 1))
            failed_names+=("$name")
        else
            echo "PASS"
            pass=$((pass + 1))
        fi
    else
        echo "FAIL (exit)"
        tail -10 "$log" | sed 's/^/    /'
        fail=$((fail + 1))
        failed_names+=("$name")
    fi
    rm -f "$log"
done

echo "----------------------------------------"
echo "通过 $pass，失败 $fail"
if [ "$fail" -ne 0 ]; then
    printf '失败用例：%s\n' "${failed_names[*]}"
    exit 1
fi
exit 0
