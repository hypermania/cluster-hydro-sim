#!/usr/bin/env python3
"""Compare two Heggie histories at identical times, without fitting shifts."""
import argparse
import json
from pathlib import Path
import numpy as np
from plot_heggie_comparison import load_run, plt


def compare(coarse, fine, output):
    end = min(coarse[-1, 1], fine[-1, 1])
    if end <= 1:
        raise ValueError("histories must extend past one initial relaxation time")
    times = np.geomspace(1, end, 3000)
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    metrics = {}
    for ax, col, label in zip(axes.flat, (2, 3, 5, 7),
                             (r"$\rho_{s,c}$", r"$\rho_{b,c}$", r"$r_c$", r"$M_c/M$")):
        ys = []
        for data, style, name in ((coarse, "-", "coarse"), (fine, "--", "fine")):
            y = np.interp(times, data[:, 1], data[:, col])
            ys.append(y)
            ax.loglog(times, y, style, label=name)
        metrics[f"column_{col}_max_relative_difference"] = float(np.max(np.abs(ys[0]/ys[1]-1)))
        ax.set(xlabel=r"$t/t_{rh}(0)$", ylabel=label)
        ax.legend();ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)
    output.with_suffix(".json").write_text(json.dumps(metrics, indent=2)+"\n")
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coarse", type=Path)
    parser.add_argument("fine", type=Path)
    parser.add_argument("output_pdf", type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(load_run(args.coarse)[0], load_run(args.fine)[0], args.output_pdf), indent=2))


if __name__ == "__main__":
    main()
