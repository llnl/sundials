"""Plot the completed 2D diffusion GPU weak scaling comparison."""

import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
GPUS = (1, 2, 4, 8)
METHODS = (
    ("C++ / CUDA", "cpp", "#1768AC", "o"),
    ("Python / Torch inplace", "inplace", "#D55E00", "s"),
    ("Python / Triton", "triton", "#008E83", "D"),
)


def read_timings():
    times = defaultdict(list)
    with (ROOT / "timings.csv").open(newline="") as source:
        for row in csv.DictReader(source):
            gpus = int(row["gpus"])
            mode = row["rhs_mode"]
            times[(gpus, mode)].append(float(row["python_seconds"]))
            # Each C++ run is paired with both Python modes in the CSV.
            if mode == "inplace":
                times[(gpus, "cpp")].append(float(row["cpp_seconds"]))

    for gpus in GPUS:
        for _, method, _, _ in METHODS:
            if len(times[(gpus, method)]) != 5:
                raise ValueError(f"Expected five {method} runs at {gpus} GPUs")
    return times


def main():
    times = read_timings()
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.titleweight": "semibold",
            "axes.labelsize": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.facecolor": "white",
        }
    )

    fig, runtime = plt.subplots(figsize=(10.5, 6.7))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.105, right=0.89, bottom=0.19, top=0.72)

    fig.text(
        0.105, 0.93, "2D diffusion  |  GPU weak scaling",
        fontsize=21, weight="bold", color="#17212B", ha="left",
    )
    fig.text(
        0.105, 0.868,
        "≈103,000 unknowns per GPU   •   ARKStep + PCG   •   five runs per configuration",
        fontsize=11, color="#596779", ha="left",
    )

    handles = []
    for label, method, color, marker in METHODS:
        medians = [statistics.median(times[(g, method)]) for g in GPUS]
        lows = [min(times[(g, method)]) for g in GPUS]
        highs = [max(times[(g, method)]) for g in GPUS]
        for gpus in GPUS:
            runtime.scatter(
                [gpus] * 5, times[(gpus, method)], marker=marker, s=18,
                color=color, alpha=0.28, linewidths=0, zorder=2,
            )
        line = runtime.errorbar(
            GPUS, medians,
            yerr=(
                [medians[i] - lows[i] for i in range(len(GPUS))],
                [highs[i] - medians[i] for i in range(len(GPUS))],
            ),
            fmt=marker + "-", color=color, linewidth=2.4, markersize=7.5,
            markeredgecolor="white", markeredgewidth=1.2, capsize=4,
            elinewidth=1.4, zorder=4,
        )
        handles.append(line.lines[0])
        runtime.annotate(
            f"{medians[-1]:.1f} s", (GPUS[-1], medians[-1]),
            xytext=(15, 0), textcoords="offset points", va="center",
            color=color, fontsize=10.5, weight="semibold",
        )

    fig.legend(
        handles, [item[0] for item in METHODS], ncol=3, loc="upper left",
        bbox_to_anchor=(0.099, 0.82), frameon=False, fontsize=10.5,
        columnspacing=2.7, handlelength=2.4,
    )

    runtime.set_title("Solve time across GPU counts", loc="left", pad=15, color="#17212B")
    runtime.set_ylabel("Seconds  ↓ lower is better")
    runtime.set_ylim(55, 180)
    runtime.set_yticks((60, 80, 100, 120, 140, 160, 180))
    runtime.set_xscale("log", base=2)
    runtime.set_xlim(0.85, 10.5)
    runtime.set_xticks(GPUS, [str(g) for g in GPUS])
    runtime.minorticks_off()
    runtime.set_xlabel("GPUs")
    runtime.grid(axis="y", color="#E5EAF0", linewidth=0.9)
    runtime.set_axisbelow(True)
    runtime.spines["left"].set_color("#CED6DF")
    runtime.spines["bottom"].set_color("#CED6DF")
    runtime.tick_params(colors="#405063", length=0, pad=7)

    fig.text(
        0.105, 0.075,
        "Dots: individual runs    •    Bars: full range    •    Lines: medians",
        fontsize=10, color="#596779", ha="left",
    )
    fig.text(
        0.105, 0.042,
        "Source: diffusion_2d_compare_482867/timings.csv  |  321 × (320 × GPUs + 1) grid",
        fontsize=9, color="#7A8796", ha="left",
    )

    for extension, options in (
        ("png", {"dpi": 240}),
        ("pdf", {}),
        ("svg", {}),
    ):
        output = ROOT / f"weak_scaling.{extension}"
        fig.savefig(output, **options)
        print(output)
    plt.close(fig)


if __name__ == "__main__":
    main()
