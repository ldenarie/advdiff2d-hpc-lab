#!/usr/bin/env bash
# Runs the benchmark matrix and stores raw results with a description of the
# machine and the build, after measuring the memory bandwidth of the machine
# with bench/stream.cpp. See docs/benchmarking.md for the protocol and the
# contract the driver must follow.
#
# Usage: bench/run_bench.sh [options]
#   -b DIR    build directory containing advdiff2d     (default: build-release)
#   -o DIR    root directory for results               (default: results)
#   -l LABEL  machine label used in the output name    (default: short hostname)
#   -s LIST   grid sizes n (n x n grid)                (default: "128 512 2048 4096")
#   -k LIST   kernels passed as kernel=...             (default: "incremental")
#   -t LIST   values of OMP_NUM_THREADS                (default: "1")
#   -r N      timed repetitions per configuration      (default: 7)
#   -w N      untimed warm-up repetitions              (default: 1)
#   -u N      target cell updates per repetition       (default: 200000000)
#   -f        accept a non-Release or sanitized build (results flagged as such)

set -euo pipefail

build_dir=build-release
out_root=results
label=$(hostname -s 2>/dev/null || echo unknown)
sizes="128 512 2048 4096"
kernels="naive weighted incremental"
threads_list="1"
reps=7
warmup=1
target_updates=200000000
force=0

while getopts "b:o:l:s:k:t:r:w:u:fh" opt; do
    case "$opt" in
        b) build_dir=$OPTARG ;;
        o) out_root=$OPTARG ;;
        l) label=$OPTARG ;;
        s) sizes=$OPTARG ;;
        k) kernels=$OPTARG ;;
        t) threads_list=$OPTARG ;;
        r) reps=$OPTARG ;;
        w) warmup=$OPTARG ;;
        u) target_updates=$OPTARG ;;
        f) force=1 ;;
        h) sed -n '2,17p' "$0"; exit 0 ;;
        *) sed -n '2,17p' "$0" >&2; exit 2 ;;
    esac
done

driver="$build_dir/advdiff2d"
cache="$build_dir/CMakeCache.txt"
[[ -x "$driver" ]] || { echo "error: $driver not found or not executable" >&2; exit 1; }
[[ -f "$cache" ]] || { echo "error: $cache not found" >&2; exit 1; }

cache_value() { sed -n "s/^$1:[A-Z]*=//p" "$cache" | head -n 1; }

build_type=$(cache_value CMAKE_BUILD_TYPE)
sanitize=$(cache_value ADVDIFF_SANITIZE)
if [[ "$build_type" != "Release" || "$sanitize" == "ON" ]]; then
    if [[ $force -eq 0 ]]; then
        echo "error: build type is '$build_type', sanitizers '${sanitize:-OFF}'." >&2
        echo "       Benchmark a Release build without sanitizers, or pass -f." >&2
        exit 1
    fi
    echo "warning: benchmarking a '$build_type' build (sanitizers '${sanitize:-OFF}')" >&2
fi

commit=$(git rev-parse --short=12 HEAD 2>/dev/null || echo unknown)
if [[ -n "$(git status --porcelain --untracked-files=no 2>/dev/null)" ]]; then
    commit="${commit}-dirty"
fi

stamp=$(date -u +%Y%m%dT%H%M%SZ)
out_dir="$out_root/${stamp}_${label}"
mkdir -p "$out_dir"
raw="$out_dir/raw.csv"
env_file="$out_dir/env.txt"

# --- Description of the machine and the build -------------------------------
{
    echo "date_utc: $stamp"
    echo "label: $label"
    echo "commit: $commit"
    echo "uname: $(uname -srvmo)"
    echo "nproc: $(nproc)"
    echo "build_dir: $build_dir"
    for key in CMAKE_BUILD_TYPE CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_RELEASE \
               ADVDIFF_SANITIZE ADVDIFF_WARNINGS_AS_ERRORS; do
        echo "$key: $(cache_value "$key")"
    done
    compiler=$(cache_value CMAKE_CXX_COMPILER)
    echo "compiler_version: $("$compiler" --version 2>/dev/null | head -n 1)"
    echo "benchmark: sizes='$sizes' kernels='$kernels' threads='$threads_list'" \
         "reps=$reps warmup=$warmup target_updates=$target_updates"
    echo "env: $(env | grep -E '^(OMP|GOMP|KMP)_' | sort | tr '\n' ' ')"
    echo
    echo "--- lscpu"
    lscpu 2>/dev/null || echo "lscpu unavailable"
    echo
    echo "--- free -m"
    free -m 2>/dev/null || echo "free unavailable"
} > "$env_file"

# --- Memory bandwidth ceiling (mini-STREAM) ---------------------------------
# Measured in the same job as the kernels, once per thread count, so that the
# ceiling describes the machine the kernels actually ran on.
stream_bin="$build_dir/stream"
if [[ -x "$stream_bin" ]]; then
    expected_stream_header="kernel,n,threads,best_s,median_s,stream_gbs,traffic_gbs,melem_s"
    echo "label,commit,$expected_stream_header" > "$out_dir/stream.csv"
    for threads in $threads_list; do
        echo "stream threads=$threads" >&2
        output=$(OMP_NUM_THREADS=$threads "$stream_bin")
        if [[ "$(head -n 1 <<< "$output")" != "$expected_stream_header" ]]; then
            echo "error: unexpected CSV header from $stream_bin" >&2
            exit 1
        fi
        tail -n +2 <<< "$output" | sed "s/^/$label,$commit,/" >> "$out_dir/stream.csv"
    done
else
    echo "warning: $stream_bin not found, no bandwidth ceiling measured" >&2
fi

# --- Benchmark matrix --------------------------------------------------------
expected_header="kernel,nx,ny,steps,rep,time_s,mlups,mass_rel_err"
echo "label,commit,threads,$expected_header" > "$raw"

for kernel in $kernels; do
    for threads in $threads_list; do
        for n in $sizes; do
            # Enough steps for each repetition to last long enough to be timed reliably.
            steps=$(( (target_updates + n * n - 1) / (n * n) ))
            echo "kernel=$kernel threads=$threads n=$n steps=$steps" >&2
            output=$(OMP_NUM_THREADS=$threads "$driver" \
                         kernel="$kernel" nx="$n" ny="$n" steps="$steps" \
                         reps="$reps" warmup="$warmup" output=csv)
            header=$(head -n 1 <<< "$output")
            if [[ "$header" != "$expected_header" ]]; then
                echo "error: unexpected CSV header from the driver:" >&2
                echo "  got:      $header" >&2
                echo "  expected: $expected_header" >&2
                exit 1
            fi
            rows=$(tail -n +2 <<< "$output")
            if [[ $(grep -c . <<< "$rows") -ne $reps ]]; then
                echo "error: expected $reps rows, got:" >&2
                echo "$rows" >&2
                exit 1
            fi
            sed "s/^/$label,$commit,$threads,/" <<< "$rows" >> "$raw"
        done
    done
done

echo "results written to $out_dir" >&2
echo "$out_dir"
