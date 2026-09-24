#!/usr/bin/env python3
"""Overlay moving-fluid histories on Heggie & Aarseth (1992), Figs. 1--4.

PDF panels retain the source aspect ratio. The source is supplied by the user,
not redistributed. No fitted time, length, density or collapse-epoch shifts.
"""
import argparse
import json
from pathlib import Path

import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

FSTARS = (.001, .002, .005, .01, .02)
# Page (zero based), scan-read frame in PDF points, axis limits and scales.
PANELS = {
    1: (7, (328.5, 57.9, 524.8, 200.2), (1, 1000), (.001, .3), "log"),
    2: (7, (328.5, 297.4, 525.5, 440.4), (1, 1000), (.01, 10), "log"),
    3: (8, (332.1, 524.9, 533.3, 671.0), (1, 2.8), (.001, 100), "linear"),
    4: (9, (327.6, 60.7, 535.1, 213.4), (1, 150), (.1, 100), "log"),
}


def load_run(directory):
    directory = Path(directory)
    raw = np.fromfile(directory / "history.dat", dtype="<f8")
    if raw.size % 13 or raw.size < 26:
        raise ValueError(f"invalid history size: {directory}")
    data = raw.reshape(-1, 13)
    if not np.all(np.isfinite(data)) or np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError(f"nonfinite or unordered history: {directory}")
    status = (directory / "status.txt").read_text().splitlines()[0]
    return data, status


def diagnostics(data, status):
    return {
        "status": status, "accepted_steps": len(data)-1,
        "final_time_trh": float(data[-1, 1]),
        "max_mass_relative_error": float(np.max(np.abs(data[:, 8]/data[0, 8]-1))),
        "max_step_energy_ledger_error": float(np.max(np.abs(data[:, 11]))),
        "max_single_central_density": float(np.max(data[:, 2])),
        "minimum_core_mass_fraction": float(np.min(data[:, 7])),
    }


def segregation_curves(data):
    # Differentiate in code time; t_rc uses exactly the same time unit.
    xi = [np.gradient(np.log(data[:, col]), data[:, 0], edge_order=2)*data[:, relax]
          for col, relax in ((2, 9), (3, 10))]
    return data[:, 4]/data[0, 4], data[:, 3]/data[:, 2], *xi


def make_page(document, number, runs):
    page, frame, xlim, ylim, xscale = PANELS[number]
    rect = fitz.Rect(frame)
    pix = document[page].get_pixmap(matrix=fitz.Matrix(3, 3), clip=rect, alpha=False)
    background = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)
    fig, ax = plt.subplots(figsize=(8.5, 7))
    ax.set(xscale=xscale, yscale="log", xlim=xlim, ylim=ylim)
    ax.set_box_aspect(rect.height/rect.width)
    ax.imshow(background, extent=(0, 1, 0, 1), transform=ax.transAxes,
              aspect="auto", alpha=.35, zorder=0)
    # Reapply box aspect after imshow: never stretch the scanned reference.
    ax.set_box_aspect(rect.height/rect.width)
    colors = plt.colormaps["rainbow"](np.linspace(.05, .9, 5))
    if number in (1, 2):
        for fstar, color in zip(FSTARS, colors):
            data, status = runs[f"single_{fstar:g}"]
            label = rf"$f^*={fstar:g}$"
            if status != "TIME_LIMIT": label += f" ({status.lower()})"
            if number == 1:
                ax.plot(data[:, 1], data[:, 7], color=color, label=label)
            else:
                ax.plot(data[:, 1], data[:, 5], color=color, label=label)
                ax.plot(data[:, 1], data[:, 6], "--", color=color)
        ax.set(xlabel=r"$t/t_{rh}(0)$", ylabel=r"$M_c/M$" if number == 1 else r"$r_c$ (solid), $r_h$ (dashed)")
    elif number == 3:
        data, status = runs["segregation"]
        x, ratio, xis, xib = segregation_curves(data)
        for y, label, color in ((ratio, r"$\rho_{b,c}/\rho_{s,c}$", "purple"),
                                (xis, r"$\xi_s$", "red"), (xib, r"$\xi_b$", "blue")):
            ax.plot(x, np.where(y > 0, y, np.nan), color=color, label=label)
        ax.set(xlabel=r"$\Phi_c/\Phi_c(0)$", ylabel=r"$\rho_{b,c}/\rho_{s,c},\ \xi_s,\ \xi_b$")
    else:
        for name, style in (("bs", "--"), ("bsbb", "-")):
            data, status = runs[name]
            for col, component, color in ((2, "singles", "blue"), (3, "binaries", "red")):
                ax.plot(data[:, 1], data[:, col], style, color=color,
                        label=f"{component}, {'BS only' if name == 'bs' else 'BS+BB'}"+
                        (f" ({status.lower()})" if status != "TIME_LIMIT" else ""))
        ax.set(xlabel=r"$t/t_{rh}(0)$", ylabel=r"$\rho_c$ ($G=M=1,\ E_0=-1/4$)")
    ax.legend(fontsize=8, loc="best", framealpha=.85)
    ax.set_title(f"Heggie & Aarseth (1992), Fig. {number}: gray reference / colored HydroSim")
    fig.text(.08, .045, "Primordial-binary gas models; no formation or ejection. C=0.104; ln Λ=10 for Fig. 4.\n"
             "Wall at 10⁴ Plummer radii; standard N-body length/density units assumed for gas curves.\n"
             "Curves end at the recorded stopping condition; no extrapolation or collapse-time alignment.", fontsize=8)
    fig.subplots_adjust(bottom=.19, top=.91, left=.13, right=.97)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_pdf", type=Path)
    parser.add_argument("run_root", type=Path, help="contains single_0.001 ... single_0.02, segregation, bs, bsbb")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.run_root / "plots"
    output.mkdir(parents=True, exist_ok=True)
    names = [f"single_{f:g}" for f in FSTARS]+["segregation", "bs", "bsbb"]
    runs = {name: load_run(args.run_root/name) for name in names}
    summary = {name: diagnostics(*run) for name, run in runs.items()}
    with fitz.open(args.reference_pdf) as document, PdfPages(output/"heggie_1_4_overlays.pdf") as pdf:
        for number in range(1, 5):
            fig = make_page(document, number, runs)
            fig.savefig(output/f"heggie_fig{number}_overlay.pdf")
            pdf.savefig(fig)
            plt.close(fig)
    (output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
