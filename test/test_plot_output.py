#!/usr/bin/env python3
"""Regression tests for script/plot_output.py."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import struct
import sys
import tempfile
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "script" / "plot_output.py"
SPEC = importlib.util.spec_from_file_location("plot_output", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
plot_output = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = plot_output
SPEC.loader.exec_module(plot_output)


class PlotOutputTest(unittest.TestCase):
    def make_output(self, directory: Path) -> None:
        n_grid = 4
        (directory / "paramNames.txt").write_text("N\nc2\nthres\n")
        (directory / "paramTypes.txt").write_text("Integer64\nReal64\nReal64\n")
        (directory / "paramOffsets.txt").write_text("0\n8\n32\n")
        parameter_data = bytearray(40)
        struct.pack_into("<q", parameter_data, 0, n_grid)
        struct.pack_into("<3d", parameter_data, 8, 0.1, 0.2, 0.3)
        struct.pack_into("<d", parameter_data, 32, 1.0e-3)
        (directory / "param.dat").write_bytes(parameter_data)

        profile_time = np.array([0.0, 0.25], dtype="<f8")
        radius = np.array(
            [[0.1, 0.3, 1.0, 3.0], [0.09, 0.28, 1.0, 3.0]], dtype="<f8"
        )
        profile_time.tofile(directory / "t.dat")
        radius.tofile(directory / "R.dat")
        for component_index, component in enumerate(plot_output.COMPONENTS):
            scale = 10.0 ** (-component_index)
            rho = scale / (1.0 + radius**2)
            energy = 0.5 * scale + 0.1 / (1.0 + radius)
            rho.astype("<f8").tofile(directory / f"Rho_{component}.dat")
            energy.astype("<f8").tofile(directory / f"U_{component}.dat")

        radii_time = np.array([0.0, 0.1, 0.25], dtype="<f8")
        fractions = np.array([0.1, 0.5], dtype="<f8")
        radii_time.tofile(directory / "radii_t.dat")
        fractions.tofile(directory / "radii_fraction.dat")
        for component_index, component in enumerate(plot_output.COMPONENTS):
            values = np.array(
                [[0.2, 1.0], [0.21, 0.98], [0.22, 0.95]], dtype="<f8"
            ) * (1.0 + 0.1 * component_index)
            values.tofile(directory / f"radii_{component}.dat")

        radii_time.tofile(directory / "central_t.dat")
        for component_index, component in enumerate(plot_output.COMPONENTS):
            rho = np.array([1.0, 1.1, 1.3], dtype="<f8") * (component_index + 1)
            energy = np.array([0.4, 0.42, 0.45], dtype="<f8")
            rho.tofile(directory / f"central_Rho_{component}.dat")
            energy.tofile(directory / f"central_U_{component}.dat")

    def test_actual_times_beyond_max_time_are_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.make_output(directory)
            for name, suffix in (("paramNames.txt", "maxTime\n"),
                                 ("paramTypes.txt", "Real64\n"),
                                 ("paramOffsets.txt", "40\n")):
                path = directory / name
                path.write_text(path.read_text() + suffix)
            path = directory / "param.dat"
            path.write_bytes(path.read_bytes() + struct.pack("<d", 0.2))
            run = plot_output.load_output(directory)
            self.assertEqual(run.parameters["maxTime"], 0.2)
            for data in (run.profiles, run.radii, run.central):
                self.assertEqual(data.time[-1], 0.25)
            for figure in (plot_output._radii_figure(run.radii, "overshoot"),
                           plot_output._central_figure(run.central, "density", "overshoot")):
                try:
                    self.assertEqual(figure.axes[0].lines[0].get_xdata()[-1], 0.25)
                    self.assertGreater(figure.axes[0].get_xlim()[1], 0.25)
                finally:
                    plot_output.plt.close(figure)

    def test_load_and_generate_every_plot_family(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            self.make_output(directory)
            run = plot_output.load_output(directory)
            self.assertEqual(run.parameters["N"], 4)
            np.testing.assert_allclose(run.parameters["c2"], [0.1, 0.2, 0.3])
            self.assertEqual(run.profiles.radius.shape, (2, 4))
            self.assertEqual(run.radii.radius["s"].shape, (3, 2))
            self.assertEqual(run.central.time.size, 3)

            profile_figure = plot_output._profile_figure(
                run.profiles, "density", "encoding check", maximum_snapshots=8
            )
            lines = profile_figure.axes[0].lines
            self.assertEqual(
                [line.get_linestyle() for line in lines],
                ["None", "--", "-", "None", "--", "-"],
            )
            self.assertEqual(
                [line.get_marker() for line in lines],
                ["*", "None", "None", "*", "None", "None"],
            )
            self.assertEqual(lines[0].get_markevery(), 1)
            for line in lines[1:3]:
                np.testing.assert_allclose(line.get_color(), lines[0].get_color())
            for line in lines[4:6]:
                np.testing.assert_allclose(line.get_color(), lines[3].get_color())
            self.assertFalse(np.allclose(lines[0].get_color(), lines[3].get_color()))
            plot_output.plt.close(profile_figure)

            self.assertEqual(
                plot_output._component_curve_style("s", 150)["markevery"], 5
            )
            self.assertNotIn(
                "markevery", plot_output._component_curve_style("b", 150)
            )

            plots = plot_output.generate_plots(
                run, directory / "figures"
            )
            self.assertEqual(
                {path.name for path in plots},
                {
                    "density_profiles.pdf",
                    "velocity_dispersion_profiles.pdf",
                    "lagrangian_radii.pdf",
                    "central_density.pdf",
                    "central_velocity_dispersion.pdf",
                },
            )
            for path in plots:
                self.assertGreater(path.stat().st_size, 1000)

            arguments = plot_output.build_argument_parser().parse_args([str(directory)])
            self.assertEqual(tuple(arguments.formats), ("pdf",))

    def test_incomplete_optional_family_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            self.make_output(directory)
            (directory / "central_U_d.dat").unlink()
            with self.assertRaisesRegex(plot_output.OutputFormatError, "incomplete"):
                plot_output.load_output(directory)


if __name__ == "__main__":
    unittest.main()
