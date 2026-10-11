#!/usr/bin/env python3
"""Fit two initial Plummer scales, never the plotted coordinates or heating.

The two target densities are approximate readings of Fig.4's solid curves
at t/trh(0)=1, with roughly 5--10% scan-reading uncertainty.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.optimize import least_squares
from plot_heggie_hydro import load_stock_run as load_run

TARGET = np.array([3.6, .6])


def at_time(data, time):
    if not data[0, 1] <= time <= data[-1, 1]:
        raise ValueError("requested epoch outside actual history")
    return np.array([np.interp(time, data[:, 1], data[:, col]) for col in (2, 3)])


def run(executable, directory, model, scales, zones, end):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError(directory)
    directory.parent.mkdir(parents=True, exist_ok=True)
    command = [str(executable), "hydro-heggie", str(directory), str(model), ".005",
               str(zones), str(end), ".001", "1000000", *map(str, scales)]
    with (directory.parent / (directory.name + ".command.json")).open("x") as record:
        record.write(json.dumps(command, indent=2))
    with (directory.parent / (directory.name + ".log")).open("x") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT,
                       env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"})
    np.testing.assert_array_equal(np.fromfile(directory/"profile"/"param.dat", dtype="<f8"), scales)
    if (directory/"profile"/"paramNames.txt").read_text().splitlines() != ["single_scale", "binary_scale_ratio"]:
        raise ValueError("profile parameter serialization changed")
    return load_run(directory)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--executable", type=Path, default=Path("main-strict"))
    parser.add_argument("--zones", type=int, default=2000)
    parser.add_argument("--full", action="store_true", help="evolve both heating models with saved fit")
    args = parser.parse_args()
    executable = args.executable.resolve()
    args.root.mkdir(parents=True, exist_ok=True)
    fit_path = args.root / "fit.json"
    if args.full:
        fit = json.loads(fit_path.read_text())
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(run, executable, args.root/name, model,
                                   fit["scales"], args.zones, 150)
                       for name, model in (("bsbb", 3), ("bs", 2))]
            for future in futures:
                future.result()
        return
    if fit_path.exists():
        raise ValueError("choose a fresh calibration directory")
    evaluations = []

    def residual(log_scales):
        scales = np.exp(log_scales)
        directory = args.root / f"fit_{len(evaluations):02d}"
        data, status = run(executable, directory, 3, scales, args.zones, 1.01)
        density = at_time(data, 1.)
        record = dict(scales=scales.tolist(), densities=density.tolist(), status=status,
                      directory=str(directory))
        evaluations.append(record)
        (args.root / "evaluations.json").write_text(json.dumps(evaluations, indent=2))
        print(record, flush=True)
        return np.log(density/TARGET)

    result = least_squares(residual, np.log([.72, 1.25]), diff_step=.003,
                           bounds=(np.log([.3, .3]), np.log([2., 3.])),
                           xtol=2e-4, ftol=2e-4, gtol=2e-4, max_nfev=20)
    error = residual(result.x)
    if np.max(np.abs(error)) > .005:
        raise RuntimeError("calibration did not match both targets within 0.5%")
    fit_path.write_text(json.dumps(dict(scales=np.exp(result.x).tolist(), target=TARGET.tolist(),
                                       relative_error=np.expm1(error).tolist(), zones=args.zones,
                                       fit_time=1., model=3, target_uncertainty="5--10% scan reading",
                                       optimizer_success=bool(result.success)), indent=2))


if __name__ == "__main__":
    main()
