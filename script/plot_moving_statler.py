#!/usr/bin/env python3
"""PDF overlays for moving-fluid Statler Figs. 11--15; no ejection-only proxy."""
import argparse
import json
from pathlib import Path

import fitz
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import plot_statler_comparison_atlas as atlas


def summarize(directory: Path, run: atlas.RunData) -> dict:
    budget = np.genfromtxt(directory / "budget.csv", delimiter=",", names=True)
    budget = np.atleast_1d(budget)
    if not budget.size or not np.all(np.isfinite(budget["mass_plus_outflow"])):
        raise ValueError("invalid moving-fluid budget")
    h = run.history
    peak = int(np.argmax(h["rho0"]))
    return {
        "final_time_trh": float(h["time"][-1]),
        "peak_time_trh": float(h["time"][peak]),
        "peak_density_statler_units": float(h["rho0"][peak]),
        "final_binary_number_in_domain": float(h["nb"][-1]),
        "cumulative_formed_binaries": float(budget["formed"][-1]),
        "mass_plus_outflow_max_relative_error": float(np.max(np.abs(
            budget["mass_plus_outflow"] / budget["mass_plus_outflow"][0] - 1))),
        "max_step_energy_ledger_error": float(np.max(budget["energy_ledger_error"])),
        "mass_fraction_outside": float(1 - h["mass"][-1] / h["mass"][0]),
        "capture_energy_code_units": float(budget["capture_energy"][-1]),
        "note": "Energy ledger checks gas energy plus discrete gravity work, not exact global gravitational-energy conservation.",
    }


def compare_runs(coarse: atlas.RunData, fine: atlas.RunData, output: Path) -> dict:
    """Compare on common physical times, without fitting the collapse epoch."""
    end = min(coarse.history["time"][-1], fine.history["time"][-1])
    if end <= 1:
        raise ValueError("convergence comparison needs histories extending past one trh")
    times = np.geomspace(1, end, 3000)
    figure, axes = atlas.plt.subplots(3, 2, figsize=(9, 10))
    metrics = {}
    for axis, (key, label) in zip(axes.flat, [
        ("rho0", r"$\rho_0/(M/r_0^3)$"), ("nb", r"$N_b$"),
        ("vms2", r"$v_{m,s}^2/(GM/r_0)$"), ("ratio", r"$\rho_{b,0}/\rho_{s,0}$"),
        ("rc", r"$r_c/r_0$"), ("rh", r"$r_h/r_0$"),
    ], strict=True):
        curves = []
        for run, style, name in [(coarse, "-", "coarse"), (fine, "--", "fine")]:
            values = np.interp(times, run.history["time"], run.history[key])
            curves.append(values)
            axis.loglog(times, values, style, label=name)
        for name, mask in [("early", times <= 20), ("late", times > 20)]:
            if np.any(mask):
                metrics[f"{key}_{name}_max_relative_difference"] = float(
                    np.max(np.abs(curves[0][mask]/curves[1][mask]-1)))
        axis.set(xlabel=r"$t/t_{rh}$", ylabel=label)
        axis.grid(alpha=.25);axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output / "resolution_comparison.pdf")
    atlas.plt.close(figure)
    (output / "resolution_comparison.json").write_text(json.dumps(metrics, indent=2)+"\n")
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_pdf", type=Path)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path, help="second resolution, compared at the same times")
    args = parser.parse_args()
    run = atlas.load_run(args.directory)
    output = args.output or args.directory / "plots"
    output.mkdir(parents=True, exist_ok=True)
    atlas.configure_style()
    note = ("Matched Plummer scales and POWER_LAW tidal capture; no encounter ejection. "
            "Local fluid conduction/heating differ from orbit-averaged Fokker--Planck kernels. "
            "Open-boundary mass loss is measured separately. No fitted time or density rescaling.")
    with fitz.open(args.reference_pdf) as document, PdfPages(output / "statler_11_15_overlays.pdf") as pdf:
        for number in range(11, 16):
            # Keep the four approximate reference-square epochs, not the
            # arbitrarily later 10000-trh endpoint. Never invent/interpolate
            # a profile at a time that was not saved by the simulation.
            data = run
            if number in (14, 15):
                indices = [0]
                for target in (18, 50, 350, 2000):
                    j = int(np.argmin(np.abs(run.snapshot_time-target)))
                    if abs(run.snapshot_time[j]/target-1) < 0.05:
                        indices.append(j)
                indices = np.unique(indices)
                data = atlas.RunData(run.history, run.snapshot_time[indices],
                                     {k: v[indices] for k, v in run.snapshots.items()})
            figure = atlas.make_page(document, number, data, data, overlay=True,
                                     model_label="MovingThreeFluidSim", note_override=note)
            figure.savefig(output / f"statler_fig{number}_overlay.pdf")
            pdf.savefig(figure)
            atlas.plt.close(figure)
    summary = summarize(args.directory, run)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if args.compare:
        print(json.dumps(compare_runs(run, atlas.load_run(args.compare), output), indent=2))


if __name__ == "__main__":
    main()
