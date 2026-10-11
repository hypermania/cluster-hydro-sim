#!/usr/bin/env python3
"""Regenerate Figs.1--4, using fitted initial profiles only for Fig.4."""
import argparse
import json
from pathlib import Path
import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from plot_heggie_hydro import FSTARS, load_stock_run as load_run, make_page
from heggie_reference import REFERENCE, sample


def density_table(reference, baseline, fitted):
    rows = [rf"{t:g} & {s:g}, {h:g} & {old[0]:.3f}, {old[1]:.3f} & {new[0]:.3f}, {new[1]:.3f} \\"
            for (t, s, h), old, new in zip(reference, baseline, fitted)]
    # Include the complete environment: LaTeX's protected \input inside a
    # tabular can leave a token before a following booktabs \bottomrule.
    header = (r"\begin{tabular}{rrrr}" "\n" r"\toprule" "\n"
              r"$\tau$ & Paper $(\rho_{c1},\rho_{c2})$ & Original & Fitted initial profile\\" "\n"
              r"\midrule" "\n")
    return header+"\n".join(rows)+"\n"+r"\bottomrule"+"\n"+r"\end{tabular}"+"\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_pdf", type=Path)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("calibrated", type=Path)
    parser.add_argument("--coarse", type=Path, help="same fitted profile at lower resolution")
    args = parser.parse_args()
    fit = json.loads((args.calibrated/"fit.json").read_text())
    a, b = fit["scales"]
    out = args.calibrated/"plots";out.mkdir(exist_ok=True)
    names = [f"single_{f:g}" for f in FSTARS]+["segregation", "bs", "bsbb"]
    original = {name: load_run(args.baseline/name) for name in names}
    runs = {**original, **{name: load_run(args.calibrated/name) for name in ("bs", "bsbb")}}
    report = dict(fit=fit, reference=REFERENCE.tolist(), runs={})
    for name in ("bs", "bsbb"):
        data, status = runs[name]
        report["runs"][name] = dict(status=status, final_time=float(data[-1, 1]),
            initial_density=data[0, 2:4].tolist(),
            density_at_one=sample(data, [1])[0].tolist(),
            initial_rh=float(data[0, 6]),
            initial_energy=float(data[0, 11]/(data[0, 8]**2*(3*np.pi/16))),
            max_mass_drift=float(np.max(np.abs(data[:, 8]/data[0, 8]-1))),
            peak_single_density=float(np.max(data[:, 2])),
            peak_single_time=float(data[np.argmax(data[:, 2]), 1]))
    calibrated = runs["bsbb"][0];baseline = original["bsbb"][0]
    values = sample(calibrated, REFERENCE[:, 0])
    base_values = sample(baseline, REFERENCE[:, 0])
    report["bsbb_comparison"] = dict(calibrated=values.tolist(), baseline=base_values.tolist(),
        calibrated_over_paper=(values/REFERENCE[:, 1:]).tolist(),
        baseline_over_paper=(base_values/REFERENCE[:, 1:]).tolist())
    if args.coarse:
        times = np.geomspace(1, 150, 1000)
        coarse = load_run(args.coarse)[0]
        report["refinement"] = dict(max_relative_difference=np.max(np.abs(
            sample(coarse, times)/sample(calibrated, times)-1), axis=0).tolist(),
            initial_density_coarse=coarse[0, 2:4].tolist(),
            sampling="1000 common logarithmic times 1..150; same fitted parameters, no refit")
    (out/"calibration_audit.json").write_text(json.dumps(report, indent=2))
    footer = (f"Stock hydrostatic Fig.4: initial Plummer scales a_s={a:.4f}, a_b/a_s={b:.4f}; f=0.03 fixed.\n"
              "Fit BS+BB densities only at t/trh(0)=1; identical fitted profile used for BS only.\n"
              "Fixed mass/length units; initial energy is NOT reset to −1/4. Gray: paper; color: simulation.")
    with fitz.open(args.reference_pdf) as doc, PdfPages(out/"heggie_1_4_overlays.pdf") as pdf:
        for number in range(1, 5):
            fig = make_page(doc, number, runs, "Hydrostatic (fitted profile)" if number==4 else "Hydrostatic (original profile)")
            if number==4:
                fig.axes[0].set_ylabel(r"$\rho_c$ (fixed baseline mass and length units)")
                fig.texts[0].set_text(footer)
            fig.savefig(out/f"heggie_fig{number}_overlay.pdf");pdf.savefig(fig);plt.close(fig)
        fig = make_page(doc, 4, {"bs": original["bsbb"], "bsbb": runs["bsbb"]}, "BS+BB: original vs fitted")
        ax = fig.axes[0]
        for line, name in zip(ax.lines, ("Original singles", "Original binaries", "Fitted singles", "Fitted binaries")):
            line.set_label(name)
        ax.legend(fontsize=8, loc="lower left")
        ax.set_ylabel(r"$\rho_c$ (fixed baseline mass and length units)")
        fig.texts[0].set_text("Colored curves: BS+BB only; dashed = original profile, solid = fitted profile.\n"
                             "Gray background: original paper panel (its dashed curves are BS only).\n"
                             "Only t/trh(0)=1 was fitted. Later times are predictions, not fit constraints.")
        fig.savefig(out/"heggie_fig4_before_after.pdf");plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    for col, (ax, component) in enumerate(zip(axes, ("Singles", "Binaries"))):
        ax.plot(REFERENCE[:, 0], base_values[:, col]/REFERENCE[:, col+1], "o--", label="Original profile")
        ax.plot(REFERENCE[:, 0], values[:, col]/REFERENCE[:, col+1], "s-", label="Fitted profile")
        ax.axhspan(.9, 1.1, color="gray", alpha=.15, label="Approximate scan uncertainty")
        ax.axhline(1, color="gray", lw=.8)
        ax.set(xscale="log", ylim=(0, max(1.3, 1.05*np.max(values/REFERENCE[:, 1:]))),
               ylabel="Simulation / paper", title=component)
        ax.legend(fontsize=8)
    axes[1].set_xlabel(r"$t/t_{rh}(0)$")
    fig.tight_layout();fig.savefig(out/"heggie_fig4_density_ratios.pdf");plt.close(fig)
    (out/"density_table.tex").write_text(density_table(REFERENCE, base_values, values))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
