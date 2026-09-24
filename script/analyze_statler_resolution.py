#!/usr/bin/env python3
"""Compare raw central dispersion histories and quantify post-collapse rebounds."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_output import read_parameters


def reference_force_bias(zones, first_face=1e-6, last_face=1e4):
    """Reconstruct the frozen initial Plummer correction, not evolved gravity.

    This is the exact finite-volume quadrature used by initialize(...,true)
    for a unit central-density Plummer profile. Ratios do not depend on its
    overall density normalization. It diagnoses a candidate truncation term;
    it does not establish the cause of an observed oscillation.
    """
    edges = np.r_[0., np.geomspace(first_face, last_face, zones)]
    radius = (edges[1:]+edges[:-1])/2
    volume = (edges[1:]**3-edges[:-1]**3)/3
    rho = (1+radius**2)**(-2.5)
    pressure = rho/(18*np.sqrt(1+radius**2))
    mass = np.cumsum(volume*rho)
    eta = (radius**3-edges[:-1]**3)/(3*volume)
    gravity = (eta*mass+(1-eta)*np.r_[0., mass[:-1]])/radius**2
    right_pressure = np.r_[pressure[1:], pressure[-1]]
    correction = (edges[1:]**2*(right_pressure-pressure)/volume+rho*gravity)/rho
    locations = np.array([.01, .1, 1., 10.])
    return {"radius": locations.tolist(), "correction_over_initial_gravity":
            np.interp(np.log(locations), np.log(radius), correction/gravity).tolist()}


def rebound_metrics(time, values, start=50., end=10000.):
    time, values = np.asarray(time), np.asarray(values)
    if (time.ndim != 1 or values.shape != time.shape or len(time) < 2 or
            not np.all(np.isfinite(time)) or not np.all(np.isfinite(values)) or
            np.any(np.diff(time) <= 0) or np.any(values <= 0)):
        raise ValueError("need ordered finite times and positive dispersions")
    end = min(end, time[-1])
    if not time[0] <= start < end:
        raise ValueError("history does not cover the analysis interval")
    # A common logarithmic-time sampling avoids weighting by adaptive step count.
    grid = np.geomspace(start, end, 10000)
    log_values = np.log(np.interp(grid, time, values))
    peaks, properties = find_peaks(log_values, prominence=np.log(1.03))
    return {
        "start_trh": start, "end_trh": end,
        "positive_log_variation": float(np.maximum(np.diff(log_values), 0).sum()),
        "rebound_times_trh": grid[peaks].tolist(),
        "rebound_prominence_fraction": np.expm1(properties["prominences"]).tolist(),
        "definition": "10000 common log-time samples; peaks with log prominence >= log(1.03). No smoothing. This is a numerical diagnostic, not a physical-mode classification.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", type=Path, nargs="+")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(8, 8))
    summary = {}
    colors = plt.colormaps["rainbow"](np.linspace(.05, .9, len(args.directories)))
    for directory, color in zip(args.directories, colors):
        time = np.fromfile(directory/"history_time_trh.dat", dtype="<f8")
        values = np.fromfile(directory/"history_vms2_0.dat", dtype="<f8")
        grid = read_parameters(directory/"grid")
        zones = int(grid["zones"])
        summary[str(zones)] = rebound_metrics(time, values)
        summary[str(zones)]["initial_force_bias"] = reference_force_bias(
            zones, grid["first_face"], grid["last_face"])
        for j, ax in enumerate(axes):
            mask = (time >= (1 if j == 0 else 50)) & (time <= (10000 if j == 0 else 1000))
            ax.loglog(time[mask], values[mask], color=color, label=f"{zones} cells")
            ax.set(ylabel=r"$v_{m,s,c}^2/(GM/r_0)$", xlabel=r"$t/t_{rh}(0)$")
            ax.grid(alpha=.2); ax.legend()
    axes[0].set_xlim(1, 10000)
    axes[1].set_xlim(50, 1000)
    axes[0].set_title("Same physical times; raw accepted-step histories")
    axes[1].set_title("Post-collapse detail (no smoothing or fitted shifts)")
    fig.tight_layout()
    fig.savefig(args.output/"central_dispersion_resolution.pdf")
    plt.close(fig)
    (args.output/"central_dispersion_resolution.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
