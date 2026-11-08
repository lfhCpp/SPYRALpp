#!/bin/bash
#----------------------------------------------------------------------
# 批量跑一组 run：读 runs 清单，逐条调用 scripts/spyralpp.sh all
#
#   scripts/run_all.sh [runs_file] [--force]
#
# runs_file 每行：<run_id> <config> [n_events]
#   - 以 # 开头的行和空行忽略
#   - run_id 相对 $SPYRALPP_OUT 建立子目录
#
# 断点续跑：已有 clouds.root / fit.root 的 run 自动跳过；--force 全部重跑。
# 结尾给出汇总；任一 run 失败则退出码非 0（方便挂在后台/定时任务里）。
#
# 环境变量同 scripts/spyralpp.sh：
#   SPYRALPP_BIN / SPYRALPP_OUT / SPYRALPP_LOG / SPYRALPP_THREADS / SPYRALPP_SHARDS
#----------------------------------------------------------------------
set -u

SELF_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SELF_DIR/.." && pwd)

runs_file=$ROOT_DIR/configs/runs.txt
force=0
for arg in "$@"; do
    case "$arg" in
        --force) force=1 ;;
        -*) echo "未知参数：$arg" >&2; exit 2 ;;
        *) runs_file=$arg ;;
    esac
done

[ -f "$runs_file" ] || { echo "找不到 runs 清单：$runs_file" >&2; exit 2; }

export SPYRALPP_FORCE=$force
total=0
failed=0
failed_ids=()

echo "== SpyralPP 批量运行：$(date '+%F %T') =="
echo "   runs 清单：$runs_file"
echo "   输出目录 ：${SPYRALPP_OUT:-$ROOT_DIR/output}"
echo "   分片/线程：${SPYRALPP_SHARDS:-1} 分片，${SPYRALPP_THREADS:-auto} 线程/进程"
echo "   force    ：$force"
echo

while IFS= read -r line || [ -n "$line" ]; do
    # 去掉注释与行尾空白
    line=${line%%#*}
    line=$(echo "$line" | sed 's/[[:space:]]*$//')
    [ -z "$line" ] && continue

    # shellcheck disable=SC2086
    set -- $line
    id=$1
    cfg=${2:-$ROOT_DIR/configs/example.cfg}
    n=${3:-100}
    total=$((total + 1))

    echo "[$total] $id  (config=$(basename "$cfg"), events=$n)"
    if ! "$SELF_DIR/spyralpp.sh" all "$id" "$cfg" "$n"; then
        failed=$((failed + 1))
        failed_ids+=("$id")
        echo "  ✗ $id 失败，继续下一个（日志在 ${SPYRALPP_LOG:-$ROOT_DIR/output/logs}）"
    fi
    echo
done <"$runs_file"

echo "== 汇总 =="
echo "   共 $total 个 run，成功 $((total - failed))，失败 $failed"
if [ "$failed" -ne 0 ]; then
    printf '   失败：%s\n' "${failed_ids[*]}"
    exit 1
fi
echo "   状态：scripts/spyralpp.sh status"
exit 0
