"""Readback, normalization and diagnostic checks for Heggie observers."""
import sys
from pathlib import Path
import unittest
import tempfile
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"script"))
import plot_heggie_comparison as plotting
from compare_heggie_resolution import compare


class HeggieTest(unittest.TestCase):
    def test_history_validation_and_resolution(self):
        data = np.ones((21, 13))
        data[:, 0] = np.arange(21)
        data[:, 1] = np.arange(21)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data.astype("<f8").tofile(root/"history.dat")
            (root/"status.txt").write_text("TIME_LIMIT\n")
            loaded, status = plotting.load_run(root)
            self.assertEqual(status, "TIME_LIMIT")
            self.assertEqual(max(compare(loaded, loaded, root/"comparison.pdf").values()), 0)
            data[1, 0] = 0
            data.astype("<f8").tofile(root/"history.dat")
            with self.assertRaisesRegex(ValueError, "unordered"):
                plotting.load_run(root)

    def test_observer_histories(self):
        root = Path(__file__).resolve().parents[1]/"output"
        for model in range(4):
            # Short C++ integration tests do not write production stop status.
            data = np.fromfile(root/f"check_heggie_{model}/history.dat", dtype="<f8").reshape(-1, 13)
            self.assertEqual(len(data), 21)
            self.assertTrue(np.all(np.isfinite(data)))
            self.assertTrue(np.all(np.diff(data[:, 0]) > 0))
            self.assertAlmostEqual(data[0, 5], 3*np.pi/16/np.sqrt(2), places=9)
            ratio = data[0, 3]/data[0, 2]
            self.assertAlmostEqual(ratio, 1e-12 if model == 0 else .01 if model == 1 else .06/.97)
            summary = plotting.diagnostics(data, "TEST")
            self.assertLess(summary["max_mass_relative_error"], 1e-12)
            self.assertLess(summary["max_step_energy_ledger_error"], 1e-12)

    def test_growth_rate_units(self):
        data = np.ones((101, 13))
        data[:, 0] = np.linspace(0, 1, 101)
        data[:, 2] = np.exp(2*data[:, 0])
        data[:, 3] = .01*np.exp(3*data[:, 0])
        data[:, 9] = 4
        data[:, 10] = 5
        _, _, xis, xib = plotting.segregation_curves(data)
        np.testing.assert_allclose(xis, 8, rtol=1e-12)
        np.testing.assert_allclose(xib, 15, rtol=1e-12)


if __name__ == "__main__":
    unittest.main()
