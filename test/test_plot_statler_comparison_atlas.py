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
    def test_all_panels_have_published_scales_and_limits(self):
        self.assertEqual(set(MODULE.PANEL_AXES), set(range(1, 24)))
        for number, specs in MODULE.PANEL_AXES.items():
            with self.subTest(figure=number):
                figure, axes = MODULE.plt.subplots(len(specs), 1, squeeze=False)
                try:
                    panels = list(axes[:, 0])
                    MODULE.apply_panel_axes(panels, number)
                    for axis, spec in zip(panels, specs, strict=True):
                        self.assertEqual(axis.get_xscale(), spec.xscale)
                        self.assertEqual(axis.get_yscale(), spec.yscale)
                        np.testing.assert_allclose(axis.get_xlim(), spec.xlim)
                        np.testing.assert_allclose(axis.get_ylim(), spec.ylim)
                finally:
                    MODULE.plt.close(figure)
        self.assertEqual(MODULE.PANEL_AXES[2][0].ylim, (1, 1e4))
        self.assertEqual(MODULE.PANEL_AXES[12][0].yscale, "log")
        self.assertEqual(MODULE.PANEL_AXES[14][0].xlim, (1e-7, 1e5))
        self.assertEqual(MODULE.PANEL_AXES[13][0].xlim, (20, 1e-2))

    def test_binary_count_masks_zero_on_log_axis(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_synthetic_run(directory)
            run = MODULE.load_run(directory)
            run.history["nb"][1] = 0
            figure, axis = MODULE.plt.subplots()
            try:
                MODULE.plot_binary_number(axis, run)
                self.assertEqual(axis.get_yscale(), "log")
                np.testing.assert_array_equal(axis.lines[0].get_xdata(), [2])
            finally:
                MODULE.plt.close(figure)

    def test_rates_use_paper_primary_unit(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_synthetic_run(directory)
            run = MODULE.load_run(directory)
            figure, axes = MODULE.plt.subplots(3, 1)
            try:
                MODULE.plot_interaction_rates(list(axes), run)
                np.testing.assert_allclose(axes[0].lines[0].get_ydata(),
                                           run.history["capture_rate"][1:2] * 225e6 / 3e5)
            finally:
                MODULE.plt.close(figure)

    def test_atlas_and_registered_overlays(self):
        # No copyrighted PDF or simulation output is required by this test.
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_synthetic_run(directory)
            run = MODULE.load_run(directory)
            reference = directory / "reference.pdf"
            with MODULE.fitz.open() as document:
                for _ in range(31):
                    document.new_page(width=612, height=792)
                # After clockwise rotation, the original bottom (blue) half
                # of Fig. 13a must be on the left, and the red half on the right.
                page_number, box, _ = MODULE.OVERLAY_FRAMES[13][0]
                left, top, right, bottom = box
                middle = (top + bottom) / 2
                page = document[page_number - 1]
                page.draw_rect((left, top, right, middle), color=None, fill=(1, 0, 0))
                page.draw_rect((left, middle, right, bottom), color=None, fill=(0, 0, 1))
                document.save(reference)
            with MODULE.fitz.open(reference) as document:
                for number in MODULE.OVERLAY_FRAMES:
                    figure = MODULE.make_page(document, number, run, run, overlay=True)
                    try:
                        figure.canvas.draw()
                        self.assertEqual(len(figure.axes), len(MODULE.PANEL_AXES[number]))
                        if number == 13:
                            pixels = figure.axes[0].images[0].get_array()
                            h, w = pixels.shape[:2]
                            np.testing.assert_array_equal(pixels[h // 2, w // 4], [0, 0, 255])
                            np.testing.assert_array_equal(pixels[h // 2, 3 * w // 4], [255, 0, 0])
                        for axis, spec in zip(figure.axes, MODULE.PANEL_AXES[number], strict=True):
                            self.assertEqual(len(axis.images), 1)
                            height, width = axis.images[0].get_array().shape[:2]
                            box = axis.get_window_extent()
                            self.assertAlmostEqual(box.height / box.width, height / width,
                                                   places=12)
                            np.testing.assert_allclose(axis.get_xlim(), spec.xlim)
                            np.testing.assert_allclose(axis.get_ylim(), spec.ylim)
                            np.testing.assert_allclose(
                                axis.images[0].get_transform().transform([[0, 0], [1, 1]]),
                                axis.transAxes.transform([[0, 0], [1, 1]]))
                    finally:
                        MODULE.plt.close(figure)
            output = directory / "atlas.pdf"
            MODULE.build_atlas(reference, run, run, output)
            for path, pages in [(output, 23), (directory / "atlas_overlays.pdf", 5)]:
                with MODULE.fitz.open(path) as document:
                    self.assertEqual(len(document), pages)

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
