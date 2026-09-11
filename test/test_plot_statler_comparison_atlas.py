#!/usr/bin/env python3
"""Regression checks for the Statler comparison-atlas data helpers."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "script" / "plot_statler_comparison_atlas.py"
SPEC = importlib.util.spec_from_file_location("plot_statler_comparison_atlas", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_synthetic_run(directory: Path) -> None:
    history = np.array([0.0, 1.0, 2.0], dtype="<f8")
    for filename in MODULE.HISTORY_FILES.values():
        values = history.copy()
        if "ratio" in filename or "rate" in filename or "rho0" in filename:
            values += 1.0
        values.tofile(directory / filename)
    np.array([4], dtype="<i8").tofile(directory / "snapshot_zone_count.dat")
    np.array([0.0, 2.0], dtype="<f8").tofile(directory / "snapshot_time_trh.dat")
    profile = np.tile(np.geomspace(1.0e-2, 1.0e1, 4), 2).astype("<f8")
    for filename in MODULE.SNAPSHOT_FILES.values():
        profile.tofile(directory / filename)


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        write_synthetic_run(directory)
        run = MODULE.load_run(directory)
        assert run.history["time"].shape == (3,)
        assert run.snapshots["radius"].shape == (2, 4)

    radius = np.geomspace(1.0e-3, 1.0e2, 100)
    density = (1.0 + radius**2) ** -2.5
    projected_radius, sigma = MODULE.projected_density(radius, density)
    analytic = (4.0 / 3.0) * (1.0 + projected_radius**2) ** -2.0
    central = projected_radius < 10.0
    maximum_relative_error = np.max(np.abs(sigma[central] / analytic[central] - 1.0))
    assert maximum_relative_error < 0.025, maximum_relative_error
    print("Statler atlas plotting helpers passed")
    return 0


class StatlerPlotTest(unittest.TestCase):
    def test_overshooting_time_is_not_clipped_to_paper_viewport(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_synthetic_run(directory)
            np.array([0.0, 1.0, 10001.0], dtype="<f8").tofile(
                directory / MODULE.HISTORY_FILES["time"])
            np.array([0.0, 10001.0], dtype="<f8").tofile(
                directory / "snapshot_time_trh.dat")
            run = MODULE.load_run(directory)
            figure, axis = MODULE.plt.subplots()
            try:
                MODULE.plot_binary_number(axis, run)
                self.assertEqual(axis.lines[0].get_xdata()[-1], 10001.0)
                self.assertEqual(axis.get_xlim()[1], 10000.0)
            finally:
                MODULE.plt.close(figure)

    def test_plotting_helpers(self):
        self.assertEqual(main(), 0)


if __name__ == "__main__":
    unittest.main()
