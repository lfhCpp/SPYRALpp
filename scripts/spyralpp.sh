#!/bin/bash
# SpyralPP 单次运行驱动
#
#   scripts/spyralpp.sh <命令> <run_id> [config] [n_events]
#
# 命令：
#   simulate   只跑模拟：clouds.root + clouds_kin.root（+ HDF5）
#   analyze    只跑重建：分片并行 → 合并成 fit.root
#   merge      只把已有分片结果合并成 fit.root
#   all        simulate + analyze + merge
#   status     打印各 run 的产物状态
#
# 环境变量（都有默认值）：
#   SPYRALPP_BIN      可执行文件目录             默认 <repo>/build
#   SPYRALPP_OUT      输出根目录                 默认 <repo>/output
#   SPYRALPP_LOG      日志目录                   默认 $SPYRALPP_OUT/logs
#   SPYRALPP_THREADS  进程内线程数（0=自动）      默认 0
#   SPYRALPP_SHARDS   重建分片数（多进程）        默认 1
#   SPYRALPP_FORCE    1 = 忽略已有产物全部重跑    默认 0
#   SPYRALPP_CLEAN_SHARDS  1 = 合并后删掉分片文件 默认 0
#
# 断点续跑：clouds.root / fit.root 存在就跳过。日志：$SPYRALPP_LOG/<run>.<stage>.log
set -u

SELF_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SELF_DIR/.." && pwd)

BIN=${SPYRALPP_BIN:-$ROOT_DIR/build}
OUT=${SPYRALPP_OUT:-$ROOT_DIR/output}
LOGDIR=${SPYRALPP_LOG:-$OUT/logs}
THREADS=${SPYRALPP_THREADS:-0}
SHARDS=${SPYRALPP_SHARDS:-1}
FORCE=${SPYRALPP_FORCE:-0}
CLEAN_SHARDS=${SPYRALPP_CLEAN_SHARDS:-0}

DEFAULT_CFG=$ROOT_DIR/configs/example.cfg

die() { echo "错误：$*" >&2; exit 2; }

usage() {
    grep '^#' "${BASH_SOURCE[0]}" | sed -n '2,27p' | sed 's/^# \{0,1\}//'
    exit 2
}

now() { date '+%F %T'; }

run_logged() {
    local log="$1"
    shift
    mkdir -p "$(dirname "$log")"
    echo "[$(now)] $*" >"$log"
    "$@" >>"$log" 2>&1
    local rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "  ✗ 失败（rc=$rc）：$*" >&2
        echo "    日志：$log" >&2
        tail -n 5 "$log" | sed 's/^/    /' >&2
    fi
    return $rc
}

check_env() {
    if [ -z "${ROOTSYS:-}" ] && ! command -v root-config >/dev/null 2>&1; then
        die "找不到 ROOT。先 source <root>/bin/thisroot.sh，或 export ROOTSYS=<root>"
    fi
    [ -x "$BIN/rk4_simulate" ] || die "找不到 $BIN/rk4_simulate（先构建：cmake -S . -B build && cmake --build build -j）"
    [ -x "$BIN/rk4_analyze" ] || die "找不到 $BIN/rk4_analyze"
}

run_dir() { echo "$OUT/$1"; }

sanitize_shards() {
    case "$SHARDS" in
        ''|*[!0-9]*) echo 1 ;;
        *) [ "$SHARDS" -ge 1 ] && echo "$SHARDS" || echo 1 ;;
    esac
}

cmd_simulate() {
    local id="$1" cfg="$2" nevents="$3"
    local dir; dir=$(run_dir "$id")
    mkdir -p "$dir"
    local clouds="$dir/clouds.root"
    if [ -s "$clouds" ] && [ "$FORCE" != "1" ]; then
        echo "  · 模拟已存在，跳过：$clouds"
        return 0
    fi
    echo "  ▸ 模拟 $id（$nevents 事件）"
    run_logged "$LOGDIR/$id.simulate.log" "$BIN/rk4_simulate" "$cfg" "$clouds" "$nevents"
}

cmd_analyze() {
    local id="$1" cfg="$2"
    local dir; dir=$(run_dir "$id")
    local clouds="$dir/clouds.root"
    [ -s "$clouds" ] || die "缺少 $clouds（先跑 simulate）"

    local fit="$dir/fit.root"
    if [ -s "$fit" ] && [ "$FORCE" != "1" ]; then
        echo "  · 重建结果已存在，跳过：$fit"
        return 0
    fi

    local n; n=$(sanitize_shards)
    echo "  ▸ 重建 $id（$n 分片，进程内线程 $THREADS）"

    local i pids=()
    for ((i = 0; i < n; ++i)); do
        local shard="$dir/fit_shard_$i.root"
        rm -f "$shard"
        if [ "$i" -eq 0 ]; then
            # 第一个分片串行跑：若配置了 solver.mesh_cache，先让它把网格缓存建好，
            # 其余分片再并行加载（配 mesh_shared = 1 即零拷贝），避免多进程同时写缓存。
            run_logged "$LOGDIR/$id.analyze.$i.log" \
                "$BIN/rk4_analyze" "$clouds" "$shard" "$cfg" "$i" "$n" || return 1
        else
            ("$BIN/rk4_analyze" "$clouds" "$shard" "$cfg" "$i" "$n" \
                >"$LOGDIR/$id.analyze.$i.log" 2>&1) &
            pids+=($!)
        fi
    done

    local failed=0
    if [ ${#pids[@]} -gt 0 ]; then
        for pid in "${pids[@]}"; do
            wait "$pid" || failed=1
        done
    fi
    if [ "$failed" -ne 0 ]; then
        echo "  ✗ 有分片失败，看 $LOGDIR/$id.analyze.*.log" >&2
        return 1
    fi
    cmd_merge "$id"
}

cmd_merge() {
    local id="$1"
    local dir; dir=$(run_dir "$id")
    local fit="$dir/fit.root"
    if [ -s "$fit" ] && [ "$FORCE" != "1" ]; then
        echo "  · 已合并，跳过：$fit"
        return 0
    fi
    local parts=("$dir"/fit_shard_*.root)
    if [ ! -e "${parts[0]}" ]; then
        echo "  ✗ 没有分片结果可合并：$dir/fit_shard_*.root" >&2
        return 1
    fi
    echo "  ▸ 合并 ${#parts[@]} 个分片 → $fit"
    run_logged "$LOGDIR/$id.merge.log" "$BIN/rk4_merge" "$fit" "${parts[@]}" || return 1
    if [ "$CLEAN_SHARDS" = "1" ]; then
        rm -f "${parts[@]}"
    fi
    return 0
}

cmd_status() {
    local id="${1:-}"
    printf '%-16s %-8s %-8s %-6s\n' "run" "clouds" "fit" "logs"
    if [ -z "$id" ]; then
        [ -d "$OUT" ] || { echo "(还没有输出目录 $OUT)"; return 0; }
        local d rid c f
        for d in "$OUT"/*/; do
            [ -d "$d" ] || continue
            rid=$(basename "$d")
            [ "$rid" = "logs" ] && continue   # 日志目录不是 run
            c="-"; f="-"
            [ -s "$d/clouds.root" ] && c="ok"
            [ -s "$d/fit.root" ] && f="ok"
            printf '%-16s %-8s %-8s %-6s\n' "$rid" "$c" "$f" \
                "$(ls "$LOGDIR/$rid".*.log 2>/dev/null | wc -l)"
        done
    else
        local c="-" f="-"
        [ -s "$(run_dir "$id")/clouds.root" ] && c="ok"
        [ -s "$(run_dir "$id")/fit.root" ] && f="ok"
        printf '%-16s %-8s %-8s %-6s\n' "$id" "$c" "$f" \
            "$(ls "$LOGDIR/$id".*.log 2>/dev/null | wc -l)"
    fi
    return 0
}

main() {
    [ $# -ge 1 ] || usage
    local cmd="$1"
    shift
    if [ "$cmd" = "status" ]; then
        cmd_status "${1:-}"
        return $?
    fi
    [ $# -ge 1 ] || usage
    local id="$1"
    local cfg="${2:-$DEFAULT_CFG}"
    local nevents="${3:-100}"
    [ -f "$cfg" ] || die "配置文件不存在：$cfg"
    check_env
    mkdir -p "$OUT" "$LOGDIR"
    case "$cmd" in
        simulate) cmd_simulate "$id" "$cfg" "$nevents" ;;
        analyze)  cmd_analyze "$id" "$cfg" ;;
        merge)    cmd_merge "$id" ;;
        all)
            cmd_simulate "$id" "$cfg" "$nevents" || exit 1
            cmd_analyze "$id" "$cfg" || exit 1
            ;;
        *) usage ;;
    esac
}

main "$@"
