#!/usr/bin/env python3
"""Summarizes and plots benchmark results produced by bench/run_bench.sh.

Usage:
  bench/plot.py size    RAW.csv [RAW.csv ...] [--bytes-per-update B]
  bench/plot.py threads RAW.csv [RAW.csv ...]

`size`    : median throughput (MLUP/s) against the grid size, one curve per
            kernel and thread count.
`threads` : speedup and parallel efficiency against the thread count, one
            curve per grid size, relative to 1 thread of the same kernel.

Each repetition is one sample. Points are medians, bars span min to max.
If a stream.csv file (written by run_bench.sh) sits next to the first CSV
file, the main-memory ceiling it measured is drawn and each result is given
as a percentage of it. The figure (PNG) and a Markdown table are written next
to the first CSV file. See docs/benchmarking.md.
"""

import argparse
import csv
import pathlib
import statistics
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def load(paths):
    rows = []
    for path in paths:
        with open(path, newline="") as f:
            rows.extend(csv.DictReader(f))
    if not rows:
        sys.exit("error: no data")
    labels = {(r["label"], r["commit"]) for r in rows}
    if len(labels) > 1:
        print(f"warning: mixing several machines or commits: {sorted(labels)}",
              file=sys.stderr)
    return rows


def load_ceilings(out_dir):
    """Returns {threads: ceiling in MLUP/s} from stream.csv (scale kernel), if any."""
    path = out_dir / "stream.csv"
    if not path.exists():
        return {}
    with open(path, newline="") as f:
        return {int(r["threads"]): float(r["melem_s"])
                for r in csv.DictReader(f) if r["kernel"] == "scale"}


def aggregate(rows):
    """Groups repetitions by configuration; returns {key: stats}."""
    groups = defaultdict(list)
    for r in rows:
        key = (r["kernel"], int(r["threads"]), int(r["nx"]), int(r["ny"]))
        groups[key].append((float(r["time_s"]), float(r["mlups"]), int(r["steps"])))
    stats = {}
    for key, samples in groups.items():
        times = [t for t, _, _ in samples]
        mlups = [m for _, m, _ in samples]
        stats[key] = {
            "n": len(samples),
            "steps": samples[0][2],
            "time_median": statistics.median(times),
            "mlups_median": statistics.median(mlups),
            "mlups_min": min(mlups),
            "mlups_max": max(mlups),
            # Spread of the repetitions, relative to the median.
            "spread": (max(mlups) - min(mlups)) / statistics.median(mlups),
        }
    return stats


def title_of(rows):
    first = rows[0]
    return f"{first['label']} @ {first['commit']}"


def plot_size(rows, stats, out_dir, bytes_per_update, ceilings):
    curves = defaultdict(list)
    for (kernel, threads, nx, ny), s in sorted(stats.items()):
        curves[(kernel, threads)].append((nx * ny, s))

    fig, ax = plt.subplots(figsize=(10, 4.5))
    for (kernel, threads), points in sorted(curves.items()):
        cells = [c for c, _ in points]
        med = [s["mlups_median"] for _, s in points]
        low = [s["mlups_median"] - s["mlups_min"] for _, s in points]
        high = [s["mlups_max"] - s["mlups_median"] for _, s in points]
        ax.errorbar(cells, med, yerr=[low, high], marker="o", capsize=3,
                    label=f"{kernel}, {threads} thread(s)")
    threads_shown = [t for t in sorted({threads for _, threads in curves}) if t in ceilings]
    for threads, style in zip(threads_shown, ["--", ":", "-.", (0, (5, 1))]):
        ax.axhline(ceilings[threads], linestyle=style, color="gray",
                   label=f"main-memory ceiling, {threads} thread(s) "
                         f"(STREAM scale, {ceilings[threads]:.0f} MLUP/s)")
    ax.set_xscale("log")
    ax.set_xlabel("grid cells (nx × ny)")
    ax.set_ylabel("throughput (MLUP/s)")
    ax.set_ylim(bottom=0)
    ax.grid(True, which="both", alpha=0.3)
    # Outside the axes: with several kernels, a legend inside hides the data.
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8)
    ax.set_title(f"Throughput vs grid size — {title_of(rows)}\n"
                 "median of repetitions, bars: min–max", fontsize=10)
    fig.tight_layout()
    png = out_dir / "throughput_vs_size.png"
    fig.savefig(png, dpi=150)

    lines = [
        f"Results for {title_of(rows)}. Effective bandwidth assumes "
        f"{bytes_per_update:g} bytes per cell update.",
    ]
    if ceilings:
        lines.append(
            "Ceiling: element rate of the STREAM scale kernel in the same job "
            "(main memory). It bounds grids larger than the last-level cache; "
            "cache-resident grids can exceed it.")
    lines += [
        "",
        "| kernel | threads | grid | steps | reps | median MLUP/s | min–max | spread "
        "| eff. GB/s | % of ceiling |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for (kernel, threads, nx, ny), s in sorted(stats.items()):
        gbs = s["mlups_median"] * bytes_per_update / 1e3
        ratio = (f"{100 * s['mlups_median'] / ceilings[threads]:.0f} %"
                 if threads in ceilings else "–")
        lines.append(
            f"| {kernel} | {threads} | {nx}×{ny} | {s['steps']} | {s['n']} "
            f"| {s['mlups_median']:.1f} | {s['mlups_min']:.1f}–{s['mlups_max']:.1f} "
            f"| {100 * s['spread']:.1f} % | {gbs:.2f} | {ratio} |")
    return png, lines


def plot_threads(rows, stats, out_dir):
    by_size = defaultdict(dict)
    for (kernel, threads, nx, ny), s in stats.items():
        by_size[(kernel, nx, ny)][threads] = s

    fig, (ax_s, ax_e) = plt.subplots(1, 2, figsize=(10, 4.2))
    max_threads = 1
    lines = [
        f"Results for {title_of(rows)}. Speedup S(p) = T(1) / T(p) with median "
        "times; efficiency E(p) = S(p) / p.",
        "",
        "| kernel | grid | threads | median time (s) | speedup | efficiency |",
        "|---|---|---|---|---|---|",
    ]
    for (kernel, nx, ny), per_threads in sorted(by_size.items()):
        if 1 not in per_threads:
            print(f"warning: no 1-thread run for {kernel} {nx}x{ny}, skipped",
                  file=sys.stderr)
            continue
        t1 = per_threads[1]["time_median"]
        ps = sorted(per_threads)
        speedup = [t1 / per_threads[p]["time_median"] for p in ps]
        efficiency = [s / p for s, p in zip(speedup, ps)]
        max_threads = max(max_threads, ps[-1])
        label = f"{kernel}, {nx}×{ny}"
        ax_s.plot(ps, speedup, marker="o", label=label)
        ax_e.plot(ps, efficiency, marker="o", label=label)
        for p, s, e in zip(ps, speedup, efficiency):
            lines.append(f"| {kernel} | {nx}×{ny} | {p} "
                         f"| {per_threads[p]['time_median']:.4f} | {s:.2f} | {e:.2f} |")

    ideal = list(range(1, max_threads + 1))
    ax_s.plot(ideal, ideal, linestyle="--", color="gray", label="ideal")
    ax_e.axhline(1.0, linestyle="--", color="gray", label="ideal")
    for ax, ylabel in ((ax_s, "speedup"), (ax_e, "parallel efficiency")):
        ax.set_xlabel("threads")
        ax.set_ylabel(ylabel)
        ax.set_xticks(ideal)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(f"Strong scaling — {title_of(rows)} (median times)", fontsize=10)
    fig.tight_layout()
    png = out_dir / "scaling_threads.png"
    fig.savefig(png, dpi=150)
    return png, lines


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mode", choices=["size", "threads"])
    parser.add_argument("csv", nargs="+", type=pathlib.Path)
    parser.add_argument("--bytes-per-update", type=float, default=24.0,
                        help="memory traffic per cell update, for the effective "
                             "bandwidth (default: 24 = read c, write c_new, "
                             "write-allocate)")
    args = parser.parse_args()

    rows = load(args.csv)
    stats = aggregate(rows)
    out_dir = args.csv[0].parent

    if args.mode == "size":
        png, lines = plot_size(rows, stats, out_dir, args.bytes_per_update,
                               load_ceilings(out_dir))
        md = out_dir / "summary_size.md"
    else:
        png, lines = plot_threads(rows, stats, out_dir)
        md = out_dir / "summary_threads.md"

    md.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote {png} and {md}", file=sys.stderr)


if __name__ == "__main__":
    main()
