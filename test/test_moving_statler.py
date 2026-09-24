"""Moving-fluid public diagnostic schema and budget summary regression."""
import sys
from pathlib import Path
import unittest
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
import plot_moving_statler as moving


class MovingStatlerTest(unittest.TestCase):
    def test_identical_convergence_histories(self):
        times = np.geomspace(1, 100, 12)
        history = {key:np.ones_like(times) for key in ("rho0","nb","vms2","ratio","rc","rh")}
        history["time"] = times
        run = moving.atlas.RunData(history, np.array([]), {})
        with tempfile.TemporaryDirectory() as directory:
            result = moving.compare_runs(run, run, Path(directory))
            self.assertEqual(max(result.values()), 0)
            self.assertTrue((Path(directory)/"resolution_comparison.pdf").is_file())

    def test_smoke_output_and_conservation(self):
        directory = Path(__file__).resolve().parents[1] / "output/check_moving_statler"
        run = moving.atlas.load_run(directory)
        result = moving.summarize(directory, run)
        self.assertEqual(len(run.history["time"]), 21)
        self.assertGreater(result["cumulative_formed_binaries"], 0)
        self.assertLess(result["capture_energy_code_units"], 0)
        self.assertLess(result["mass_plus_outflow_max_relative_error"], 1e-12)
        self.assertLess(result["max_step_energy_ledger_error"], 1e-12)
        self.assertTrue(np.all(np.diff(run.snapshots["radius"], axis=1) > 0))
        self.assertLess(run.snapshots["radius"][0, 0], 0.01)


if __name__ == "__main__":
    unittest.main()
